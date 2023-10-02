import polars
from polars.type_aliases import IntoExpr
from polars.utils.udfs import _get_shared_lib_location

lib = _get_shared_lib_location(__file__)

@polars.api.register_expr_namespace("variant_id")
class VariantId:
    def __init__(self, expr: polars.Expr):
        self._expr = expr

    def compute(self, ref: polars.Expr, alt: polars.Expr, max_pos: int) -> polars.Expr:
        return self._expr._register_plugin(
            lib=lib,
            args=[ref, alt, max_pos],
            symbol="compute",
        )


__version__: str = "0.2.0"
__all__: list[str] = ["VariantId"]
