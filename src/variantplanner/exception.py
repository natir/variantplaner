"""All custom exception could be generate by VariantPlanner."""

# std import
from __future__ import annotations

# 3rd party import


class NotAVCFError(Exception):
    """Exception raise if file read seems not be a vcf."""
