from __future__ import annotations

from dataclasses import dataclass
import random
import time

from alignment import (
    Alignment,
    GAP,
    compatible_cut_positions,
    complete_alignment_from_prefix,
    consumed_counts,
    is_valid_alignment,
    random_alignment,
    remove_full_gap_columns,
    repair_alignment,
)
from scoring import score_alignment


@dataclass(frozen=True)
class GAParameters:
    population_size: int = 50
    generations: int = 150
    tournament_size: int = 3
    crossover_probability: float = 0.85
    mutation_probability: float = 0.25
    seed: int | None = 42


@dataclass(frozen=True)
class GAResult:
    best_alignment: Alignment
    best_score: int
    initial_best_score: int
    generations: int
    elapsed_seconds: float


def run_genetic_algorithm(sequences: list[str], parameters: GAParameters) -> GAResult:
    rng = random.Random(parameters.seed)
    start_time = time.perf_counter()
    population = [random_alignment(sequences, rng) for _ in range(parameters.population_size)]
    initial_best = max(score_alignment(individual) for individual in population)
    best_alignment = max(population, key=score_alignment)
    best_score = score_alignment(best_alignment)

    for _ in range(parameters.generations):
        next_population = [copy_alignment(best_alignment)]

        while len(next_population) < parameters.population_size:
            first_parent = tournament_selection(population, parameters.tournament_size, rng)
            second_parent = tournament_selection(population, parameters.tournament_size, rng)

            if rng.random() < parameters.crossover_probability:
                child = crossover(first_parent, second_parent, sequences, rng)
            else:
                child = copy_alignment(first_parent)

            if rng.random() < parameters.mutation_probability:
                child = shift_mutation(child, sequences, rng)

            child = repair_alignment(child, sequences, rng)
            next_population.append(child)

        population = next_population
        generation_best = max(population, key=score_alignment)
        generation_score = score_alignment(generation_best)

        if generation_score > best_score:
            best_alignment = copy_alignment(generation_best)
            best_score = generation_score

    elapsed_seconds = time.perf_counter() - start_time

    return GAResult(
        best_alignment=best_alignment,
        best_score=best_score,
        initial_best_score=initial_best,
        generations=parameters.generations,
        elapsed_seconds=elapsed_seconds,
    )


def tournament_selection(population: list[Alignment], tournament_size: int, rng: random.Random) -> Alignment:
    selected = rng.sample(population, k=min(tournament_size, len(population)))
    return copy_alignment(max(selected, key=score_alignment))


def crossover(
    first_parent: Alignment,
    second_parent: Alignment,
    sequences: list[str],
    rng: random.Random,
) -> Alignment:
    if not first_parent or not first_parent[0]:
        return random_alignment(sequences, rng)

    first_width = len(first_parent[0])
    cut_in_first = rng.randint(0, first_width)
    target_counts = consumed_counts(first_parent, cut_in_first)
    compatible_positions = compatible_cut_positions(second_parent, target_counts)

    prefix = [row[:cut_in_first] for row in first_parent]

    if compatible_positions:
        cut_in_second = rng.choice(compatible_positions)
        child = [
            first_row[:cut_in_first] + second_row[cut_in_second:]
            for first_row, second_row in zip(first_parent, second_parent)
        ]
        child = remove_full_gap_columns(child)

        if is_valid_alignment(child, sequences):
            return child

    return complete_alignment_from_prefix(sequences, prefix, rng)


def shift_mutation(alignment: Alignment, sequences: list[str], rng: random.Random) -> Alignment:
    mutated = copy_alignment(alignment)
    gap_blocks = find_gap_blocks(mutated)

    if not gap_blocks:
        return mutated

    row_index, start, end = rng.choice(gap_blocks)
    direction = rng.choice([-1, 1])
    row = mutated[row_index]

    if direction == -1 and start > 0 and row[start - 1] != GAP:
        row[start - 1], row[end] = row[end], row[start - 1]
    elif direction == 1 and end + 1 < len(row) and row[end + 1] != GAP:
        row[start], row[end + 1] = row[end + 1], row[start]

    mutated = remove_full_gap_columns(mutated)

    if is_valid_alignment(mutated, sequences):
        return mutated

    return alignment


def find_gap_blocks(alignment: Alignment) -> list[tuple[int, int, int]]:
    blocks: list[tuple[int, int, int]] = []

    for row_index, row in enumerate(alignment):
        column_index = 0

        while column_index < len(row):
            if row[column_index] != GAP:
                column_index += 1
                continue

            start = column_index
            while column_index + 1 < len(row) and row[column_index + 1] == GAP:
                column_index += 1

            blocks.append((row_index, start, column_index))
            column_index += 1

    return blocks


def copy_alignment(alignment: Alignment) -> Alignment:
    return [row.copy() for row in alignment]
