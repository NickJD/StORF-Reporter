"""Small compatibility helpers shared by command adapters."""

from __future__ import annotations

import argparse


def parse_bool(value: str | bool) -> bool:
    """Parse strict legacy CLI Boolean values without using ``eval``."""

    if isinstance(value, bool):
        return value
    lowered = value.casefold()
    if lowered in {"true", "1", "yes", "y"}:
        return True
    if lowered in {"false", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError("expected True or False")
