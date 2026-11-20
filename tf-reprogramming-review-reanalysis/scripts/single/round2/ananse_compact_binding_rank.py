"""Emit the archive's globally ranked binding table without a Dask string shuffle."""
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
from loguru import logger
from scipy.stats import rankdata
from sklearn.preprocessing import minmax_scale


def rank_chromosomes(directory, destination):
    files = sorted(Path(directory).glob("chr*.csv"))
    assert files, directory
    # Use the same Pandas C parser and float64 weight parsing as the stock Dask CSV reader.
    frames = [pd.read_csv(path, dtype={"tf_target":"string[pyarrow]", "weighted_binding":"float64"}) for path in files]
    full = pd.concat(frames, ignore_index=True, copy=False)
    del frames
    assert list(full.columns)==["tf_target","weighted_binding"], full.columns
    keys = full["tf_target"].array._data
    weight = full["weighted_binding"].to_numpy(copy=False)
    assert keys.null_count==0 and np.isfinite(weight).all()
    # Both transforms are those in the archive. Their min ties are independent of input order.
    value = minmax_scale(rankdata(weight,method="min")).astype(np.float32)
    order = pc.sort_indices(keys).to_numpy(zero_copy_only=False)
    assert len(order)==len(value)>0
    last = None
    with Path(destination).open("w",encoding="utf-8",newline="") as stream:
        for start in range(0,len(order),500_000):
            ix = order[start:start+500_000]
            names = pc.take(keys,pa.array(ix)).to_pylist()
            if last is not None:
                assert names[0]>last
            assert len(names)==len(set(names)), "Archive binding output contains duplicate regulator–target keys"
            last = names[-1]
            pd.DataFrame({"tf_target":names,"weighted_binding":value[ix]}).to_csv(stream, index=False, header=start==0)
    logger.info(f"Validated compact rank aggregation of {len(order)} distinct archive binding pairs from {len(files)} contigs")
    return len(order)
