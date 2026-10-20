
def extract_annotations(record):
    annotations = []

    for feature in record.features:
        if feature.type not in {"gene", "CDS", "tRNA", "rRNA"}:
            continue

        qualifiers = feature.qualifiers

        annotations.append({
            "Feature": feature.type,
            "Start": int(feature.location.start) + 1,
            "End": int(feature.location.end),
            "Strand": feature.location.strand,
            "Gene": ", ".join(qualifiers.get("gene", ["N/A"])),
            "Product": ", ".join(
                qualifiers.get("product", ["N/A"])
            )
        })

    return annotations
