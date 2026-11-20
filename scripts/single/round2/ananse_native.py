"""Launch the archived CLI with explicit RNG and resource settings and record its native influence gate."""
import importlib
import os
import random
import runpy
import sys
from pathlib import Path

import numpy as np

from r2lib import BASE, FOLLOW, ROOT, write_tsv


def deny_reference_open(event, args):
    if event != "open" or not isinstance(args[0], (bytes, str, os.PathLike)):
        return
    path = os.path.abspath(os.fsdecode(args[0]))
    inaccessible = (BASE / "references", BASE / "evaluation", ROOT / "evaluation", FOLLOW / "evaluation")
    if any(path == str(p) or path.startswith(str(p) + "/") for p in inaccessible):
        raise PermissionError(f"Reference access blocked during prediction: {path}")


def main():
    sys.addaudithook(deny_reference_open)
    random.seed(20261118)
    np.random.seed(20261118)
    script, stage, *trailing = sys.argv[1:]
    sys.argv = [script, stage, *trailing]
    if stage == "binding" and os.environ.get("R2_MOTIF_CACHE_PREFIX"):
        from ananse_motif_cache import use_task_scan_cache
        use_task_scan_cache()
    if stage == "network":
        if os.environ.get("R2_ANANSE_PARTITIONED_NETWORK") == "1":
            from ananse_partitioned_network import install_compatible_partitioning
            install_compatible_partitioning()
        implementation = importlib.import_module("ananse.commands.network")
        stock_cluster = implementation.LocalCluster

        def cluster_with_limit(*args, **kwargs):
            kwargs["memory_limit"] = os.environ["R2_DASK_WORKER_MEMORY"]
            return stock_cluster(*args, **kwargs)

        implementation.LocalCluster = cluster_with_limit
    if stage == "influence":
        module = importlib.import_module("ananse.influence")
        Influence = module.Influence
        native_edge_filter = module._filter_network_edges
        filters = []

        def record_native_edge_cutoff(df, key, n_edges, whitelist=None):
            selected = native_edge_filter(df, key, n_edges, whitelist)
            field = module.GRN_COLUMNS.get(key,key)
            boundary = float(selected[field].min())
            tied_total_count = int((df[field] == boundary).sum())
            selected_count = int((selected[field] == boundary).sum())
            filters.append({"state":"source" if not filters else "target", "total_edges":len(df),"kept_edges":len(selected),
                            "requested_cutoff":n_edges,"boundary_prob":boundary,
                            "edges_tied_at_boundary":tied_total_count,"tied_edges_selected":selected_count,
                            "split_boundary_tie":0 < selected_count < tied_total_count})
            write_tsv(Path(os.environ["R2_EDGE_GATE"]),filters)
            return selected

        module._filter_network_edges = record_native_edge_cutoff

        native_run = Influence.run_target_score

        def record_native_gate(self):
            coverage = []
            for sym, degree in self.grn.out_degree():
                expr = self.expression_change.get(sym)
                coverage.append({"symbol": sym, "differential_network_out_degree": degree,
                                 "native_DE_present": expr is not None,
                                 "native_real_fold_change": "NA" if expr is None else expr.realfc,
                                 "native_DE_target_score": "NA" if expr is None else expr.score})
            write_tsv(Path(os.environ["R2_INFLUENCE_GATE"]), sorted(coverage, key=lambda r: r["symbol"]))
            return native_run(self)

        Influence.run_target_score = record_native_gate
    runpy.run_path(script, run_name="__main__")


if __name__ == "__main__":
    main()
