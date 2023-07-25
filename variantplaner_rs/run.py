import polars
from variantplaner_rs import VariantId

df = polars.DataFrame(
    data={
        "real_pos": list(range(10, 110, 10)),
        "ref": ["A", "C", "T", "G", "A", "C", "T", "G", "A", "C"],
        "alt": ["T", "G", "A", "C", "T", "G", "A", "C", "T", "G"],
    },
    schema={
        "real_pos": polars.UInt64,
        "ref": polars.Utf8,
        "alt": polars.Utf8,
    },
)

lf = df.lazy()

lf_id = lf.with_columns(
    id = polars.col("real_pos").variant_id.compute(polars.col("ref"), polars.col("alt"), polars.lit(120).cast(polars.UInt64))
)

print(lf_id.collect())
