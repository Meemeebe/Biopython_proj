
VALID_DNA_CHARACTERS = set("ACGTUNRYKMSWBDHVX-")


def validate_sequence(sequence):
    """Validate a DNA sequence and return any invalid characters."""
    sequence = str(sequence).upper()

    invalid = set(sequence) - VALID_DNA_CHARACTERS
    return sorted(invalid)


def analyze_sequence(sequence):
    """Calculate basic statistics for a DNA sequence."""
    sequence = str(sequence).upper()

    invalid = validate_sequence(sequence)

    if invalid:
        raise ValueError(
            f"Invalid sequence characters: {', '.join(invalid)}"
        )

    if not sequence:
        raise ValueError("The DNA sequence is empty.")

    counts = {
        base: sequence.count(base)
        for base in "ATGC"
    }

    canonical_count = sum(counts.values())

    gc_content = (
        (counts["G"] + counts["C"]) / canonical_count * 100
        if canonical_count else None
    )

    return {
        "length": len(sequence),
        "counts": counts,
        "gc_content": gc_content,
        "canonical_count": canonical_count,
        "ambiguous_count": len(sequence) - canonical_count
    }
