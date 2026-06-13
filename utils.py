from __future__ import annotations

from alignment import Alignment, alignment_to_strings


def format_alignment(alignment: Alignment) -> str:
    return "\n".join(alignment_to_strings(alignment))


def improvement_percent(initial_score: int, final_score: int) -> float:
    if initial_score == 0:
        return 0.0 if final_score == 0 else 100.0

    return ((final_score - initial_score) / abs(initial_score)) * 100.0
