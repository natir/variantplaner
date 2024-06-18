import pathlib
import polars
from polars.type_aliases import IntoExpr
from polars.plugins import register_plugin_function

@polars.api.register_expr_namespace("variant_id")
class VariantId:
    def __init__(self, expr: polars.Expr):
        self._expr = expr

    def compute(self, ref: polars.Expr, alt: polars.Expr, max_pos: int) -> polars.Expr:
        return register_plugin_function(
            plugin_path=pathlib.Path(__file__).parent,
            function_name="compute",
            args=[self._expr, ref, alt, max_pos],
        )

    def partition(self) -> polars.Expr:
        return register_plugin_function(
            plugin_path=pathlib.Path(__file__).parent,
            function_name="partition",
            args=[self._expr],
        )


__version__: str = "0.3.0"
__all__: list[str] = ["VariantId"]
