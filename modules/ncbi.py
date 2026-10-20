
from Bio import Entrez, SeqIO
from io import StringIO


def fetch_sequence(accession, email):
    Entrez.email = email

    with Entrez.efetch(
        db="nuccore",
        id=accession,
        rettype="gb",
        retmode="text"
    ) as handle:
        content = handle.read()

    record = SeqIO.read(StringIO(content), "genbank")

    return record
