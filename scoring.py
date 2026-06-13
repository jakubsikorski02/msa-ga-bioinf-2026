from __future__ import annotations

from alignment import Alignment, GAP


MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_OPEN_SCORE = -2
GAP_EXTEND_SCORE = -1


def score_alignment(alignment: Alignment) -> int:
    if not alignment:
        return 0

    sequence_count = len(alignment)
    width = len(alignment[0])
    score = 0

    for first_index in range(sequence_count - 1):
        for second_index in range(first_index + 1, sequence_count):
            previous_gap_side: int | None = None

            for column_index in range(width):
                first_symbol = alignment[first_index][column_index]
                second_symbol = alignment[second_index][column_index]

                if first_symbol == GAP and second_symbol == GAP:
                    previous_gap_side = None
                    continue

                if first_symbol != GAP and second_symbol != GAP:
                    score += MATCH_SCORE if first_symbol == second_symbol else MISMATCH_SCORE
                    previous_gap_side = None
                    continue

                current_gap_side = 1 if first_symbol == GAP else 2
                score += GAP_EXTEND_SCORE if previous_gap_side == current_gap_side else GAP_OPEN_SCORE
                previous_gap_side = current_gap_side

    return score
