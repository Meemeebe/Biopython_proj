# Gene Sequence Analysis & Annotation Platform

A Python-based bioinformatics web application for exploring DNA
sequences using **Biopython**, **Streamlit**, and **Matplotlib**. The
platform provides an interactive interface for basic sequence
statistics, GC-content visualization, open reading frame (ORF)
detection, GenBank feature inspection, pairwise sequence alignment, and
optional sequence retrieval from NCBI.

## Features

-   **FASTA and GenBank upload** --- load a single sequence record from
    a supported file.
-   **Sequence statistics** --- display sequence ID, description,
    length, nucleotide counts, ambiguous-base count, and GC content.
-   **GC-content visualization** --- calculate and plot GC content using
    a configurable sliding window.
-   **ORF analysis** --- search six reading frames for candidate ORFs
    beginning with `ATG` and ending at an in-frame stop codon (`TAA`,
    `TAG`, or `TGA`), then translate them into protein sequences.
-   **GenBank annotation viewer** --- display supported `gene`, `CDS`,
    `tRNA`, and `rRNA` features, including coordinates and qualifiers
    when available.
-   **Pairwise alignment** --- compare two DNA sequences using
    Biopython's `PairwiseAligner`.
-   **NCBI Entrez retrieval** --- retrieve a nucleotide record by
    accession number, subject to internet access and NCBI usage
    guidelines.
-   **Downloadable report** --- download a text report containing basic
    sequence information.

> **Note:** This is an educational bioinformatics project. ORF detection
> and sequence similarity are exploratory analyses and should not be
> interpreted as clinical or diagnostic results.

## Technology Stack

-   Python 3
-   [Biopython](https://biopython.org/)
-   [Streamlit](https://streamlit.io/)
-   Matplotlib
-   Pandas

## Project Structure

``` text
Biopython_proj/
├── app.py
├── requirements.txt
├── README.md
├── modules/
│   ├── __init__.py
│   ├── sequence_analysis.py
│   ├── visualization.py
│   ├── orf_analysis.py
│   ├── annotation.py
│   ├── alignment.py
│   └── ncbi.py
└── sample_data/
    ├── sample.fasta
    └── sample.gb
```

The exact files in `sample_data/` may vary depending on which sample
files you have created.

## Installation

### 1. Clone the repository

``` bash
git clone https://github.com/meemee/Biopython_proj.git
cd Biopython_proj
```

If your repository has a different URL, replace the URL above with your
repository's HTTPS URL.

### 2. Create and activate a virtual environment (recommended)

**Windows PowerShell:**

``` powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
```

**macOS/Linux:**

``` bash
python3 -m venv .venv
source .venv/bin/activate
```

### 3. Install dependencies

Make sure `requirements.txt` contains:

``` text
biopython
streamlit
matplotlib
pandas
```

Then run:

``` bash
python -m pip install -r requirements.txt
```

### 4. Start the application

From the project root directory, run:

``` bash
streamlit run app.py
```

Streamlit will display a local URL in the terminal, usually
`http://localhost:8501`. Open that address in your browser.

## How to Use

1.  Start the application.
2.  Upload a **single-record FASTA** (`.fasta`, `.fa`, `.fna`) or
    **GenBank** (`.gb`, `.gbk`) file.
3.  Review the sequence ID, description, length, nucleotide composition,
    and GC content.
4.  Choose a sliding-window size to visualize GC content along the
    sequence.
5.  Select **Find ORFs** to inspect candidate ORFs and their translated
    proteins.
6.  Upload a GenBank record to inspect its annotated gene and
    coding-sequence features.
7.  Enter two DNA sequences in the alignment section and select **Align
    Sequences**.
8.  To retrieve a record online, enter an NCBI accession number and a
    valid contact email in the NCBI section.
9.  Download the basic sequence analysis report.

## Input Formats

### FASTA

FASTA files contain a header beginning with `>` followed by the
nucleotide sequence. Example:

``` text
>example_sequence Example DNA sequence
ATGAAAGCTTTTGGATAA
```

### GenBank

GenBank files contain sequence data and may include annotations such as
genes, coding sequences (CDS), product names, organism details, and
feature coordinates. FASTA files generally do not include these feature
annotations.

Use a valid GenBank file for testing the annotation viewer.

## Analysis Notes and Limitations

-   The application currently expects one sequence record per uploaded
    file.
-   GC content is calculated from unambiguous `A`, `T`, `G`, and `C`
    bases; ambiguous characters are excluded from the GC denominator.
-   GC windows containing ambiguous bases are skipped in the
    sliding-window plot.
-   ORF detection reports candidate open reading frames in all six
    reading frames. A predicted ORF does not necessarily represent a
    real gene.
-   Pairwise alignment scores depend on the scoring parameters
    configured in the alignment module; a score is not the same as a
    percent identity.
-   NCBI retrieval requires internet access and is subject to NCBI
    service policies and rate limits.
-   The downloadable report currently includes basic sequence statistics
    and the DNA sequence. It may be extended to include ORF, annotation,
    and alignment results.

## Troubleshooting

-   **`ModuleNotFoundError`**: activate your virtual environment and
    install dependencies with
    `python -m pip install -r requirements.txt`.
-   **Streamlit command not found**: run
    `python -m streamlit run app.py`.
-   **Invalid or unsupported sequence**: confirm that the uploaded file
    is a valid FASTA or GenBank file and contains one sequence record.
-   **No ORFs found**: the input may not contain a complete
    start-to-stop ORF in the searched reading frames.
-   **NCBI request fails**: check the accession number, internet
    connection, and contact email, and try again later if the service is
    unavailable.

## Future Enhancements

-   Export complete reports including ORFs, annotations, and alignment
    results.
-   Add interactive genome-feature maps.
-   Support multiple-record FASTA and GenBank files.
-   Add richer alignment metrics and sequence comparison visualizations.
-   Improve error handling and automated tests.
-   Add an option to export extracted annotations as CSV.

## Learning Outcomes

This project demonstrates practical use of Python and Biopython for
sequence parsing, nucleotide composition analysis, GC-content
visualization, ORF detection and translation, GenBank feature
extraction, pairwise sequence alignment, and retrieval of nucleotide
records from NCBI.

## License

No license has been specified yet. Add a license file if you intend to
permit others to reuse, modify, or distribute this project.
