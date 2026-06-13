from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, pstdev

from ga import GAParameters, run_genetic_algorithm
from utils import improvement_percent


@dataclass(frozen=True)
class Instance:
    name: str
    sequences: list[str]


@dataclass(frozen=True)
class ExperimentRow:
    parameter: str
    value: str
    average_initial_score: float
    average_final_score: float
    best_final_score: int
    standard_deviation: float
    average_improvement: float
    average_time: float


INSTANCES = [
    Instance("PDF", ["AGTCGTAG", "ATCGTCG", "GTAG", "GTAGAG"]),
    Instance("Small 1", ["ACGTACGT", "ACGTCGT", "ACGACGT", "CGTACG"]),
    Instance("Small 2", ["GATTACA", "GATACA", "GACTACA", "TTACTA"]),
]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Parameter experiments for DNA MSA genetic algorithm.")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--generations", type=int, default=80)
    parser.add_argument("--output", default="parameter_tests.csv")
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    rows = run_parameter_tests(arguments.repeats, arguments.generations)
    write_csv(rows, Path(arguments.output))
    print_results(rows)


def run_parameter_tests(repeats: int, generations: int) -> list[ExperimentRow]:
    tests = [
        ("population_size", [20, 40, 60]),
        ("mutation_probability", [0.10, 0.25, 0.40]),
        ("crossover_probability", [0.60, 0.85, 1.00]),
    ]
    rows: list[ExperimentRow] = []

    for parameter, values in tests:
        for value in values:
            final_scores: list[int] = []
            initial_scores: list[int] = []
            improvements: list[float] = []
            times: list[float] = []

            for instance_index, instance in enumerate(INSTANCES):
                for repeat_index in range(repeats):
                    seed = 1000 + instance_index * 100 + repeat_index
                    parameters = create_parameters(parameter, value, generations, seed)
                    result = run_genetic_algorithm(instance.sequences, parameters)
                    initial_scores.append(result.initial_best_score)
                    final_scores.append(result.best_score)
                    improvements.append(improvement_percent(result.initial_best_score, result.best_score))
                    times.append(result.elapsed_seconds)

            rows.append(
                ExperimentRow(
                    parameter=parameter,
                    value=str(value),
                    average_initial_score=mean(initial_scores),
                    average_final_score=mean(final_scores),
                    best_final_score=max(final_scores),
                    standard_deviation=pstdev(final_scores),
                    average_improvement=mean(improvements),
                    average_time=mean(times),
                )
            )

    return rows


def create_parameters(parameter: str, value: int | float, generations: int, seed: int) -> GAParameters:
    base = {
        "population_size": 50,
        "generations": generations,
        "tournament_size": 3,
        "crossover_probability": 0.85,
        "mutation_probability": 0.25,
        "seed": seed,
    }
    base[parameter] = value
    return GAParameters(**base)


def write_csv(rows: list[ExperimentRow], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow([
            "parameter",
            "value",
            "average_initial_score",
            "average_final_score",
            "best_final_score",
            "standard_deviation",
            "average_improvement_percent",
            "average_time_seconds",
        ])

        for row in rows:
            writer.writerow([
                row.parameter,
                row.value,
                f"{row.average_initial_score:.2f}",
                f"{row.average_final_score:.2f}",
                row.best_final_score,
                f"{row.standard_deviation:.2f}",
                f"{row.average_improvement:.2f}",
                f"{row.average_time:.4f}",
            ])


def print_results(rows: list[ExperimentRow]) -> None:
    for row in rows:
        print(
            f"{row.parameter}={row.value}: "
            f"avg_initial={row.average_initial_score:.2f}, "
            f"avg_final={row.average_final_score:.2f}, "
            f"best={row.best_final_score}, "
            f"std={row.standard_deviation:.2f}, "
            f"improvement={row.average_improvement:.2f}%, "
            f"time={row.average_time:.4f}s"
        )


if __name__ == "__main__":
    main()
