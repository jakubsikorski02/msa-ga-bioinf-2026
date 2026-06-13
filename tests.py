from __future__ import annotations

from scoring import score_alignment


def test_pdf_dna_example() -> None:
    alignment = [
        "AGTCGTAG",
        "A-TCGTCG",
        "----GTAG",
        "-GTAG-AG",
    ]

    score = score_alignment(alignment)

    print("Alignment:")
    for row in alignment:
        print(row)

    print()
    print(f"SP-score: {score}")


def main() -> None:
    test_pdf_dna_example()


if __name__ == "__main__":
    main()
