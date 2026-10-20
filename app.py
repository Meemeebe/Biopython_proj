
import streamlit as st
from Bio import SeqIO, Entrez
from io import StringIO

from modules.sequence_analysis import (
    analyze_sequence,
    validate_sequence,
)
from modules.visualization import create_gc_plot
from modules.orf_analysis import find_orfs
from modules.annotation import extract_annotations
from modules.alignment import align_sequences
from modules.ncbi import fetch_sequence


# --------------------------------------------------
# PAGE CONFIGURATION
# --------------------------------------------------

st.set_page_config(
    page_title="Gene Sequence Analysis & Annotation Platform",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 Gene Sequence Analysis & Annotation Platform")

st.write(
    "Analyze DNA sequences, detect ORFs, visualize GC content, "
    "inspect annotations, and compare sequences."
)

st.divider()


# --------------------------------------------------
# MODULE 1: FILE UPLOAD AND SEQUENCE STATISTICS
# --------------------------------------------------

st.header("1. Upload and Analyze a Sequence")

uploaded_file = st.file_uploader(
    "Upload a FASTA or GenBank file",
    type=["fasta", "fa", "fna", "gb", "gbk"],
)

if uploaded_file is not None:

    try:
        content = uploaded_file.getvalue().decode("utf-8")
        filename = uploaded_file.name.lower()

        if filename.endswith((".fasta", ".fa", ".fna")):
            file_format = "fasta"
        else:
            file_format = "genbank"

        record = SeqIO.read(
            StringIO(content),
            file_format,
        )

        sequence = str(record.seq).upper()

        # Validate the sequence.
        invalid = validate_sequence(sequence)

        if invalid:
            st.error(
                "Invalid sequence characters: "
                + ", ".join(invalid)
            )
            st.stop()

        if not sequence:
            st.error("The uploaded sequence is empty.")
            st.stop()

        stats = analyze_sequence(sequence)

        st.success("Sequence loaded successfully!")

        # Display the main statistics.
        col1, col2, col3 = st.columns(3)

        col1.metric("Sequence ID", record.id)
        col2.metric("Length", f"{stats['length']} bp")

        if stats["gc_content"] is not None:
            col3.metric(
                "GC Content",
                f"{stats['gc_content']:.2f}%",
            )
        else:
            col3.metric("GC Content", "N/A")

        st.subheader("Sequence Description")
        st.write(record.description)

        st.subheader("Nucleotide Composition")

        cols = st.columns(4)

        for col, base in zip(cols, "ATGC"):
            percentage = (
                stats["counts"][base] / stats["length"] * 100
            )

            col.metric(
                f"Base {base}",
                stats["counts"][base],
                f"{percentage:.2f}% of total sequence",
            )

        st.write(
            f"Ambiguous bases: {stats['ambiguous_count']}"
        )

        with st.expander("View DNA Sequence"):
            st.code(sequence)

        st.divider()


        # --------------------------------------------------
        # MODULE 2: SLIDING-WINDOW GC CONTENT
        # --------------------------------------------------

        st.header("2. GC Content Visualization")

        window_size = st.number_input(
            "Select sliding-window size (bp)",
            min_value=1,
            max_value=len(sequence),
            value=min(100, len(sequence)),
            step=1,
        )

        gc_fig = create_gc_plot(
            sequence,
            int(window_size),
        )

        if gc_fig is not None:
            st.pyplot(gc_fig)
            import matplotlib.pyplot as plt
            plt.close(gc_fig)
        else:
            st.warning(
                "No valid windows are available for plotting. "
                "Check the sequence and window size."
            )

        st.divider()


        # --------------------------------------------------
        # MODULE 3: ORF DETECTION AND PROTEIN TRANSLATION
        # --------------------------------------------------

        st.header("3. Open Reading Frame (ORF) Analysis")

        if st.button("Find ORFs", key="find_orfs"):

            try:
                orfs = find_orfs(sequence)

                if not orfs:
                    st.info(
                        "No complete ORFs were found in the "
                        "six reading frames."
                    )
                else:
                    st.success(
                        f"Found {len(orfs)} candidate ORFs."
                    )

                    st.dataframe(
                        [
                            {
                                "Strand": orf["strand"],
                                "Frame": orf["frame"],
                                "Start": orf["start"],
                                "End": orf["end"],
                                "Length (bp)": orf["length"],
                                "Protein": orf["protein"],
                            }
                            for orf in orfs
                        ],
                        use_container_width=True,
                    )

                    st.subheader("ORF Details")

                    for index, orf in enumerate(orfs, start=1):
                        with st.expander(f"ORF {index}"):

                            st.write(
                                f"Strand: {orf['strand']}"
                            )
                            st.write(
                                f"Reading frame: {orf['frame']}"
                            )
                            st.write(
                                f"Start: {orf['start']}"
                            )
                            st.write(
                                f"End: {orf['end']}"
                            )
                            st.write(
                                f"Length: {orf['length']} bp"
                            )

                            st.write("DNA sequence")
                            st.code(orf["dna"])

                            st.write("Translated protein")
                            st.code(orf["protein"])

            except ValueError as error:
                st.warning(str(error))

        st.divider()


        # --------------------------------------------------
        # MODULE 4: GENBANK ANNOTATION
        # --------------------------------------------------

        st.header("4. Gene Annotation")

        if file_format == "genbank":

            annotations = extract_annotations(record)

            st.write("Sequence name:", record.name)
            st.write(
                "Organism:",
                record.annotations.get("organism", "Unknown"),
            )

            if annotations:
                st.dataframe(
                    annotations,
                    use_container_width=True,
                )
            else:
                st.info(
                    "No supported gene, CDS, tRNA or rRNA "
                    "features were found."
                )

        else:
            st.info(
                "Upload a GenBank file to view existing "
                "gene and CDS annotations. FASTA files "
                "usually contain sequence data only."
            )

        st.divider()


        # --------------------------------------------------
        # MODULE 5: PAIRWISE SEQUENCE ALIGNMENT
        # --------------------------------------------------

        st.header("5. Pairwise Sequence Alignment")

        sequence1 = st.text_area(
            "First DNA sequence",
            placeholder="Enter a DNA sequence",
            key="alignment_sequence1",
        )

        sequence2 = st.text_area(
            "Second DNA sequence",
            placeholder="Enter another DNA sequence",
            key="alignment_sequence2",
        )

        if st.button("Align Sequences"):

            seq1 = "".join(sequence1.upper().split())
            seq2 = "".join(sequence2.upper().split())

            if not seq1 or not seq2:
                st.warning("Enter both DNA sequences.")

            elif (
                set(seq1) - set("ACGT")
                or set(seq2) - set("ACGT")
            ):
                st.warning(
                    "Alignment currently accepts only A, C, G and T."
                )

            else:
                try:
                    alignment = align_sequences(seq1, seq2)

                    st.metric(
                        "Alignment Score",
                        alignment.score,
                    )

                    st.code(str(alignment))

                except Exception as error:
                    st.error(
                        f"Alignment failed: {error}"
                    )

        st.divider()


        # --------------------------------------------------
        # MODULE 6: DOWNLOAD REPORT
        # --------------------------------------------------

        st.header("6. Download Sequence Report")

        gc_text = (
            f"{stats['gc_content']:.2f}%"
            if stats["gc_content"] is not None
            else "N/A"
        )

        report = f"""
GENE SEQUENCE ANALYSIS & ANNOTATION PLATFORM
=============================================

Sequence ID: {record.id}
Description: {record.description}
File format: {file_format.upper()}

SEQUENCE STATISTICS
-------------------
Length: {stats['length']} bp
GC content: {gc_text}

Nucleotide counts:
A: {stats['counts']['A']}
T: {stats['counts']['T']}
G: {stats['counts']['G']}
C: {stats['counts']['C']}

Ambiguous bases: {stats['ambiguous_count']}

DNA SEQUENCE
------------
{sequence}
"""

        st.download_button(
            label="Download Analysis Report",
            data=report,
            file_name="sequence_analysis_report.txt",
            mime="text/plain",
        )

    except Exception as error:
        st.error(
            "Unable to process the uploaded file. "
            "Please check that it is a valid, single-record "
            "FASTA or GenBank file."
        )
        st.caption(f"Details: {error}")


# --------------------------------------------------
# MODULE 7: NCBI SEQUENCE RETRIEVAL
# --------------------------------------------------

st.divider()
st.header("7. Fetch a Sequence from NCBI")

st.write(
    "Retrieve a nucleotide record from NCBI using an accession number."
)

accession = st.text_input(
    "NCBI accession number",
    placeholder="Example: NM_007294",
)

email = st.text_input(
    "Email address for NCBI Entrez",
    placeholder="your_email@example.com",
)

if st.button("Fetch from NCBI"):

    if not accession.strip() or not email.strip():
        st.warning(
            "Enter both an accession number and an email address."
        )

    else:
        try:
            with st.spinner("Fetching sequence from NCBI..."):

                ncbi_record = fetch_sequence(
                    accession.strip(),
                    email.strip(),
                )

            st.success("NCBI record retrieved successfully!")

            st.write("Record ID:", ncbi_record.id)
            st.write("Description:", ncbi_record.description)
            st.write(
                "Organism:",
                ncbi_record.annotations.get("organism", "Unknown"),
            )
            st.write(
                "Sequence length:",
                f"{len(ncbi_record.seq)} bp",
            )

            with st.expander("View Retrieved Sequence"):
                st.code(str(ncbi_record.seq))

            with st.expander("View GenBank Annotations"):
                annotations = extract_annotations(ncbi_record)

                if annotations:
                    st.dataframe(
                        annotations,
                        use_container_width=True,
                    )
                else:
                    st.info("No supported features were found.")

            st.download_button(
                "Download Retrieved GenBank Record",
                data=ncbi_record.format("genbank"),
                file_name=f"{ncbi_record.id}.gb",
                mime="text/plain",
            )

        except Exception as error:
            st.error(
                "Could not retrieve the record. Check the accession "
                "number, internet connection and email address."
            )
            st.caption(f"Details: {error}")


st.divider()
st.caption(
    "Gene Sequence Analysis & Annotation Platform | "
    "Python, Biopython, Streamlit and Matplotlib"
)
