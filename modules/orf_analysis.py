
from Bio.Seq import Seq

START_CODON = "ATG"
STOP_CODONS = {"TAA", "TAG", "TGA"}


def find_orfs(dna_sequence):
    """
    Find ORFs in all six reading frames.

    Coordinates returned are 1-based and inclusive.
    Reverse-strand coordinates refer to the original sequence.
    """
    dna = str(dna_sequence).upper()

    if not dna:
        return []

    if any(base not in "ACGT" for base in dna):
        raise ValueError(
            "ORF analysis requires a DNA sequence containing "
            "only A, C, G and T."
        )

    length = len(dna)
    strands = [
        ("+", dna),
        ("-", str(Seq(dna).reverse_complement()))
    ]

    results = []

    for strand, sequence in strands:
        for frame in range(3):
            starts = []

            for pos in range(frame, len(sequence) - 2, 3):
                codon = sequence[pos:pos + 3]

                if codon == START_CODON:
                    starts.append(pos)

                elif codon in STOP_CODONS:
                    for start in starts:
                        # Include the stop codon in the ORF coordinates
                        end = pos + 3

                        if strand == "+":
                            original_start = start + 1
                            original_end = end
                        else:
                            original_start = length - end + 1
                            original_end = length - start

                        coding_sequence = sequence[start:end]
                        protein = str(
                            Seq(coding_sequence).translate(to_stop=True)
                        )

                        results.append({
                            "strand": strand,
                            "frame": frame + 1,
                            "start": original_start,
                            "end": original_end,
                            "length": end - start,
                            "protein": protein,
                            "dna": coding_sequence
                        })

                    # The stop codon ends all ORFs started
                    # in this reading frame.
                    starts = []

    results.sort(key=lambda item: (item["start"], item["strand"]))
    return results
