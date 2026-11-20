"""Bound memory in the archived full network, leaving its per-edge formula intact."""
import inspect
import math
import os
import textwrap
from pathlib import Path

import dask.dataframe as dd
import numpy as np
import pandas as pd


def cartesian_rows(regulators, targets):
    from ananse import SEPARATOR

    product = regulators.merge(targets, how="outer")
    product["tf_target"] = product["tf"] + SEPARATOR + product["target"]
    return product[["tf", "tf_expression", "target_expression", "tf_target"]]


def create_expression_network(self, expression, tfs, column="tpm"):
    # These are the archive's full-rank transforms. Only their Cartesian product is split.
    expression[column] = np.log2(expression[column] + 1e-5, dtype=np.float32)
    tf_expr = expression[expression.index.isin(tfs)]
    tf_fname = self._save_temp_expression(tf_expr, "tf", column)
    target_fname = self._save_temp_expression(expression, "target", column)
    regulators = pd.read_parquet(tf_fname)
    targets = pd.read_parquet(target_fname)
    max_regulators = int(os.getenv("R2_ANANSE_CARTESIAN_REGULATORS_PER_PARTITION", "24"))
    count = max(1, math.ceil(len(regulators)/max_regulators))
    meta = cartesian_rows(regulators.head(1), targets.head(1)).iloc[:0]
    partitions = dd.from_pandas(regulators, npartitions=count, sort=False)
    return partitions.map_partitions(cartesian_rows, targets, meta=meta, align_dataframes=False).clear_divisions()


def materialize_network(result, output_cols, outfile):
    # to_csv consumes partitions in Dask order and avoids holding the entire joined frame in the driver.
    result[output_cols].to_csv(outfile, sep="\t", index=False, single_file=True)


def install_compatible_partitioning():
    import ananse.network as native

    original = native.Network.run_network
    source = textwrap.dedent(inspect.getsource(original))
    full_compute = "    result = result.compute()\n"
    dense_export = "        result[output_cols].to_csv(outfile, sep=\"\\t\", index=False)"
    memory_export = "        r2_materialize_network(result, output_cols, outfile)"
    assert source.count(full_compute) == source.count(dense_export) == 1
    source = source.replace(full_compute, "").replace(dense_export, memory_export)
    source = source.replace("        return result[output_cols]\n", "        return result[output_cols].compute()\n")
    namespace = {**original.__globals__,"r2_materialize_network":materialize_network}
    exec(compile(source,"<ANANSE archived run_network with partitioned export>","exec"),namespace)
    native.Network.run_network = namespace[original.__name__]
    native.Network.create_expression_network = create_expression_network

    archived_aggregate = native.Network.aggregate_binding
    stock_after_chromosomes = archived_aggregate
    if os.getenv("R2_ANANSE_COMPACT_RANK") == "1":
        from ananse_compact_binding_rank import rank_chromosomes

        original_source = textwrap.dedent(inspect.getsource(archived_aggregate))
        anchor = '    ddf = dd.read_csv(os.path.join(tmpdir, "*.csv")).set_index("tf_target")'
        assert original_source.count(anchor)==1
        bounded_tail = '''    tmpfile = os.path.join(tmpdir, "binding.csv")
    r2_rank_chromosomes(tmpdir, tmpfile)
    return dd.read_csv(tmpfile).set_index("tf_target", sorted=True)
'''
        source = original_source[:original_source.index(anchor)] + bounded_tail
        scope = {**archived_aggregate.__globals__, "r2_rank_chromosomes":rank_chromosomes}
        exec(compile(source,"<ANANSE native binding aggregation with compact global rank>","exec"),scope)
        stock_after_chromosomes = scope[archived_aggregate.__name__]

    def aggregate_binding(self,*args,**kwargs):
        cache = os.getenv("R2_ANANSE_RANKED_BINDING_CSV")
        if cache and Path(cache).is_file():
            native.logger.info(f"Reusing the verified archived ranked binding table: {cache}")
            return dd.read_csv(cache).set_index("tf_target", sorted=True)
        chromosomes = os.getenv("R2_ANANSE_ARCHIVE_CHROMOSOME_CACHE")
        if chromosomes and os.getenv("R2_ANANSE_COMPACT_RANK")=="1":
            target = Path(os.environ["TMPDIR"])/"cached_ranked_binding.csv"
            native.logger.info(f"Using the verified archived chromosome binding tables: {chromosomes}")
            rank_chromosomes(chromosomes,target)
            return dd.read_csv(str(target)).set_index("tf_target", sorted=True)
        return stock_after_chromosomes(self,*args,**kwargs)

    native.Network.aggregate_binding = aggregate_binding
