"""Reuse identical native PWM scores for ANANSE binding and sampled MOAP rows.

The archive scans the same fixed regions repeatedly. Native GimmeMotifs uses
the same per-sequence function and genome background for any call with the
same median width, motifs and options. A miss executes the unchanged scanner.
"""
import hashlib
import importlib
import io
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

from r2lib import now, read_json, sha, write_json


def sha_payload(value):
    return hashlib.sha256(pickle.dumps(value, protocol=5)).hexdigest()


class NativeScanCache:
    def __init__(self, prefix, audit):
        from gimmemotifs.scanner import scan_regionfile_to_table

        self.native = scan_regionfile_to_table
        self.path = Path(prefix)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.audit = Path(audit)
        self.entries = {}
        self.history = []

    def identity(self, genome, width):
        from diskcache import Cache
        from gimmemotifs.config import CACHE_DIR

        gc_bins = [(0.0, 0.2), (0.8, 1.0)]
        gc_bins += [(round(b, 2), round(b + 0.05, 2)) for b in np.arange(0.2, 0.799, 0.05)]
        gc_bins.sort()
        native_key = f"d{genome}:{int(width)}:True:{str(gc_bins)}"
        with Cache(CACHE_DIR) as disk:
            record = disk.get(native_key)
        if record is None:
            return None
        background, bins = record
        return sha_payload((native_key, background.ids, background.seqs, bins))

    def record(self, row):
        self.history.append(row)
        write_json(self.audit, {"stage": "native scanner lookups", "data_dependent_ref_labels": False,
                               "events": self.history})

    def __call__(self, input_table, genome, scoring, pfmfile=None, ncpus=None, zscore=True,
                 gc=True, random_state=None, progress=None):
        from gimmemotifs.motif import read_motifs
        from gimmemotifs.scanner import as_fasta

        argv = (input_table, genome, scoring, pfmfile, ncpus, zscore, gc, random_state, progress)
        if scoring != "score" or not zscore or not gc or pfmfile is None or input_table.endswith("feather"):
            return self.native(*argv)
        motifs = read_motifs(pfmfile) if isinstance(pfmfile, (str, Path)) else list(pfmfile)
        motif_ids = [m.id for m in motifs]
        if not motif_ids or len(motif_ids) != len(set(motif_ids)):
            return self.native(*argv)
        checksum = sha_payload(sorted((m.id, m.hash) for m in motifs))
        frame = pd.read_table(input_table, index_col=0, comment="#")
        requested = list(frame.index)
        if not requested or len(requested) != len(set(requested)):
            return self.native(*argv)
        rng = np.random if random_state is None else random_state
        incoming = rng.get_state()
        # Consume the stock scanner's RNG for the width sample, in the same requested row order.
        sample = rng.choice(requested, size=1000, replace=False) if len(requested) >= 1000 else requested
        width = int(np.median([len(seq) for seq in as_fasta(sample, genome=genome).seqs]))
        background_hash = self.identity(genome, width)
        signature = (str(genome), scoring, bool(zscore), bool(gc), width, checksum)
        key = sha_payload(signature)
        table = self.path.with_name(f"{self.path.name}_{key}.pickle")
        manifest = table.with_suffix(".json")
        cached = self.entries.get(key)
        if cached is None and table.is_file() and manifest.is_file() and background_hash:
            meta = read_json(manifest)
            if meta["background_hash"] == background_hash and sha(table) == meta["feature_sha256"]:
                cached = (pd.read_pickle(table), meta)
                self.entries[key] = cached
        if cached is not None:
            reference, meta = cached
            if meta["background_hash"] == background_hash and set(requested) <= set(reference.index):
                scores = reference.loc[requested, motif_ids]
                assert scores.index.equals(frame.index)
                self.record({"when_UTC":now(),"operation":"exact_native_scores_reused","n_rows":len(requested),
                             "median_width_bp":width,"n_motifs":len(motif_ids),"background_hash":background_hash,
                             "native_feature_sha256":meta["feature_sha256"],"ncpus":ncpus,
                             "state_signal_included_in_cache":False,"native_width_sampling_rng_consumed":True})
                return scores
        # A cache miss executes the archive from its untouched incoming random state.
        rng.set_state(incoming)
        result = self.native(*argv)
        assert list(result.columns) == motif_ids and result.index.equals(frame.index)
        bg = self.identity(genome, width)
        if bg is not None and (cached is None or set(cached[0].index) <= set(result.index)):
            pending = table.with_suffix(".partial")
            result.to_pickle(pending, protocol=5)
            pending.replace(table)
            meta = {"native_signature_sha256":key,"background_hash":bg,"feature_sha256":sha(table),
                    "cached_n":len(result),"cache_path":str(table),"version_scope":self.path.name}
            write_json(manifest,meta)
            self.entries[key] = (result,meta)
        self.record({"when_UTC":now(),"operation":"original_native_scan_executed","n_rows":len(result),"median_width_bp":width,
                     "n_motifs":len(motif_ids),"background_hash":bg,"native_scoring":"GimmeMotifs stock score/zscore/gc",
                     "ncpus":ncpus, "cached":key in self.entries})
        return result


def use_task_scan_cache():
    binding = importlib.import_module("ananse.peakpredictor")
    activity = importlib.import_module("gimmemotifs.maelstrom.moap")
    cache = NativeScanCache(os.environ["R2_MOTIF_CACHE_PREFIX"], os.environ["R2_MOTIF_CACHE_AUDIT"])
    binding.scan_regionfile_to_table = cache
    activity.scan_regionfile_to_table = cache
