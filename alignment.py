from __future__ import annotations

import random
from typing import Iterable


Alignment = list[list[str]]

DNA_ALPHABET = {"A", "C", "G", "T"}
GAP = "-"


def validate_dna_sequences(sequences: list[str]) -> None:
    if len(sequences) < 2:
        raise ValueError("At least two sequences are required.")

    for index, sequence in enumerate(sequences, start=1):
        if not sequence:
            raise ValueError(f"Sequence {index} is empty.")

        invalid_symbols = set(sequence.upper()) - DNA_ALPHABET
        if invalid_symbols:
            symbols = ", ".join(sorted(invalid_symbols))
            raise ValueError(f"Sequence {index} contains invalid DNA symbols: {symbols}.")


def normalize_sequences(sequences: Iterable[str]) -> list[str]:
    return [sequence.strip().upper() for sequence in sequences if sequence.strip()]


def alignment_to_strings(alignment: Alignment) -> list[str]:
    return ["".join(row) for row in alignment]


def strings_to_alignment(rows: list[str]) -> Alignment:
    if not rows:
        raise ValueError("Alignment cannot be empty.")

    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("All alignment rows must have the same length.")

    return [list(row) for row in rows]


def remove_full_gap_columns(alignment: Alignment) -> Alignment:
    if not alignment or not alignment[0]:
        return alignment

    kept_columns = [
        column_index
        for column_index in range(len(alignment[0]))
        if any(row[column_index] != GAP for row in alignment)
    ]

    return [[row[column_index] for column_index in kept_columns] for row in alignment]


def sequence_from_alignment_row(row: list[str]) -> str:
    return "".join(symbol for symbol in row if symbol != GAP)


def is_valid_alignment(alignment: Alignment, sequences: list[str]) -> bool:
    if len(alignment) != len(sequences):
        return False

    if not alignment:
        return False

    width = len(alignment[0])
    if width == 0 or any(len(row) != width for row in alignment):
        return False

    for row, sequence in zip(alignment, sequences):
        if sequence_from_alignment_row(row) != sequence:
            return False

    for column_index in range(width):
        if all(row[column_index] == GAP for row in alignment):
            return False

    return True


def consumed_counts(alignment: Alignment, end_column: int) -> list[int]:
    return [
        sum(1 for symbol in row[:end_column] if symbol != GAP)
        for row in alignment
    ]


def compatible_cut_positions(alignment: Alignment, target_counts: list[int]) -> list[int]:
    positions: list[int] = []
    width = len(alignment[0])

    for cut_position in range(width + 1):
        if consumed_counts(alignment, cut_position) == target_counts:
            positions.append(cut_position)

    return positions


def random_alignment(sequences: list[str], rng: random.Random) -> Alignment:
    positions = [0 for _ in sequences]
    rows = [[] for _ in sequences]

    while any(position < len(sequence) for position, sequence in zip(positions, sequences)):
        available = [
            index
            for index, (position, sequence) in enumerate(zip(positions, sequences))
            if position < len(sequence)
        ]
        selected = [
            index
            for index in available
            if rng.random() < 0.55
        ]

        if not selected:
            selected = [rng.choice(available)]

        for row_index, sequence in enumerate(sequences):
            if row_index in selected and positions[row_index] < len(sequence):
                rows[row_index].append(sequence[positions[row_index]])
                positions[row_index] += 1
            else:
                rows[row_index].append(GAP)

    return remove_full_gap_columns(rows)


def complete_alignment_from_prefix(
    sequences: list[str],
    prefix: Alignment,
    rng: random.Random,
) -> Alignment:
    positions = [
        sum(1 for symbol in row if symbol != GAP)
        for row in prefix
    ]
    rows = [row.copy() for row in prefix]

    while any(position < len(sequence) for position, sequence in zip(positions, sequences)):
        available = [
            index
            for index, (position, sequence) in enumerate(zip(positions, sequences))
            if position < len(sequence)
        ]
        selected = [
            index
            for index in available
            if rng.random() < 0.55
        ]

        if not selected:
            selected = [rng.choice(available)]

        for row_index, sequence in enumerate(sequences):
            if row_index in selected and positions[row_index] < len(sequence):
                rows[row_index].append(sequence[positions[row_index]])
                positions[row_index] += 1
            else:
                rows[row_index].append(GAP)

    return remove_full_gap_columns(rows)


def repair_alignment(alignment: Alignment, sequences: list[str], rng: random.Random) -> Alignment:
    repaired_rows: Alignment = []

    for row, sequence in zip(alignment, sequences):
        sequence_index = 0
        repaired_row: list[str] = []

        for symbol in row:
            if symbol == GAP:
                repaired_row.append(GAP)
            elif sequence_index < len(sequence):
                repaired_row.append(sequence[sequence_index])
                sequence_index += 1

        while sequence_index < len(sequence):
            repaired_row.append(sequence[sequence_index])
            sequence_index += 1

        repaired_rows.append(repaired_row)

    width = max(len(row) for row in repaired_rows)
    for row in repaired_rows:
        row.extend([GAP] * (width - len(row)))

    repaired_rows = remove_full_gap_columns(repaired_rows)

    if is_valid_alignment(repaired_rows, sequences):
        return repaired_rows

    return random_alignment(sequences, rng)
