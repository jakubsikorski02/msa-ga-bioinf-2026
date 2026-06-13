from __future__ import annotations

import argparse

from alignment import normalize_sequences, strings_to_alignment, validate_dna_sequences
from ga import GAParameters, run_genetic_algorithm
from scoring import score_alignment
from utils import format_alignment, improvement_percent


DEFAULT_SEQUENCES = [
    "AGTCGTAG",
    "ATCGTCG",
    "GTAG",
    "GTAGAG",
]

PDF_ALIGNMENT = [
    "AGTCGTAG",
    "A-TCGTCG",
    "----GTAG",
    "-GTAG-AG",
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Multiple Sequence Alignment for DNA using a genetic algorithm.")
    parser.add_argument("--sequences", nargs="*", default=DEFAULT_SEQUENCES)
    parser.add_argument("--population-size", type=int, default=50)
    parser.add_argument("--generations", type=int, default=150)
    parser.add_argument("--tournament-size", type=int, default=3)
    parser.add_argument("--crossover-probability", type=float, default=0.85)
    parser.add_argument("--mutation-probability", type=float, default=0.25)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--pdf-alignment", action="store_true")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()

    if arguments.pdf_alignment:
        alignment = strings_to_alignment(PDF_ALIGNMENT)
        print("Alignment from PDF:")
        print(format_alignment(alignment))
        print()
        print(f"SP-score: {score_alignment(alignment)}")
        return

    sequences = normalize_sequences(arguments.sequences)
    validate_dna_sequences(sequences)

    parameters = GAParameters(
        population_size=arguments.population_size,
        generations=arguments.generations,
        tournament_size=arguments.tournament_size,
        crossover_probability=arguments.crossover_probability,
        mutation_probability=arguments.mutation_probability,
        seed=arguments.seed,
    )

    result = run_genetic_algorithm(sequences, parameters)

    print("Input sequences:")
    for sequence in sequences:
        print(sequence)

    print()
    print("Best alignment:")
    print(format_alignment(result.best_alignment))
    print()
    print(f"SP-score: {result.best_score}")
    print(f"Initial best SP-score: {result.initial_best_score}")
    print(f"Improvement: {improvement_percent(result.initial_best_score, result.best_score):.2f}%")
    print(f"Generations: {result.generations}")
    print(f"Runtime: {result.elapsed_seconds:.4f} s")
    print(f"Verification score: {score_alignment(result.best_alignment)}")


if __name__ == "__main__":
    main()
