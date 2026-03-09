import pandas as pd

# 1. The 55 gene symbols from the user's file
symbols = [
    "CHEK2", "H2AX", "RPS3", "CBX5", "CBX3", "PARP1", "XRCC5", "PRKDC", "BCLAF1", "ATAD3A",
    "ERCC6L2", "LMNA", "MACROD1", "MAPT", "MTA1", "NPM1", "NUCKS1", "OTUB1", "PARK7", "PPP1CA",
    "RAD52", "RPA1", "SFN", "UBA1", "UBE2N", "UBE2V2", "USP7", "VCP", "XRCC6", "MEN1",
    "HTATSF1", "TOP2A", "XRCC1", "ERCC3", "TAF1", "SFPQ", "MTREX", "LIG3", "HMGA2", "XPC",
    "TP53BP1", "ATR", "SNW1", "MDC1", "CHD4", "NIPBL", "POGZ", "RBBP6", "SMC5", "DDB2",
    "CDC5L", "NABP2", "CHCHD6", "INIP", "ZMYND8"
]

# 2. Classification logic based on Repairtoire images & expert knowledge
# Priorities: Image evidence > Known Pathway > General Function
classification_map = {
    # --- Core Pathways (Image Evidence) ---
    "XRCC5": {"Class": "NHEJ", "Detail": "NHEJ (Core); Signaling"},
    "XRCC6": {"Class": "NHEJ", "Detail": "NHEJ (Core); Signaling"},
    "PRKDC": {"Class": "NHEJ", "Detail": "NHEJ (Core); Signaling"},
    "PARP1": {"Class": "BER", "Detail": "BER (Sensor)"},
    "XRCC1": {"Class": "BER", "Detail": "BER (Scaffold)"},
    "LIG3": {"Class": "BER", "Detail": "BER (Ligase)"},
    "XPC": {"Class": "NER", "Detail": "NER (Global Genome Recognition)"},
    "DDB2": {"Class": "NER", "Detail": "NER (UV Damage Recognition)"},
    "ERCC3": {"Class": "NER", "Detail": "NER (Helicase XPB)"},
    "RAD52": {"Class": "HR", "Detail": "HR (Strand Annealing)"},
    "RPA1": {"Class": "Multiple (HR/NER/Signal)", "Detail": "Replication/Repair Hub (ssDNA Binding)"},
    "H2AX": {"Class": "Signaling", "Detail": "DDR Signaling (Histone H2AFX)"},
    "TP53BP1": {"Class": "Signaling/NHEJ", "Detail": "DDR Signaling; Pathway Choice"},
    "ATR": {"Class": "Signaling", "Detail": "DDR Signaling (Replication Stress)"},
    "MDC1": {"Class": "Signaling", "Detail": "DDR Signaling (Mediator)"},
    "CHEK2": {"Class": "Signaling", "Detail": "Checkpoint Kinase"},
    "UBE2V2": {"Class": "DDT/TLS", "Detail": "Damage Tolerance (MMS2)"},

    # --- Regulators & Other Functional Groups ---
    "UBE2N": {"Class": "DDT/TLS", "Detail": "Damage Tolerance (UBC13 partner of MMS2)"},
    "CHD4": {"Class": "Chromatin Remodeling", "Detail": "NuRD Complex"},
    "MTA1": {"Class": "Chromatin Remodeling", "Detail": "NuRD Complex"},
    "CBX3": {"Class": "Chromatin", "Detail": "Heterochromatin Protein (HP1 gamma)"},
    "CBX5": {"Class": "Chromatin", "Detail": "Heterochromatin Protein (HP1 alpha)"},
    "HMGA2": {"Class": "Chromatin", "Detail": "Architectural Transcription Factor"},
    "TOP2A": {"Class": "Topological Stress", "Detail": "Topoisomerase II alpha"},
    "VCP": {"Class": "Protein Turnover", "Detail": "Ubiquitin-dependent Extraction"},
    "USP7": {"Class": "Protein Turnover", "Detail": "Deubiquitinase (p53/MDM2 regulation)"},
    "UBA1": {"Class": "Protein Turnover", "Detail": "E1 Ubiquitin-activating Enzyme"},
    "OTUB1": {"Class": "Protein Turnover", "Detail": "Deubiquitinase (DNA damage response)"},
    "NPM1": {"Class": "Chaperone/Regulation", "Detail": "Nucleophosmin (Histone Chaperone)"},
    "BCLAF1": {"Class": "Transcription/Splicing", "Detail": "BCL2-associated Transcription Factor"},
    "SFPQ": {"Class": "Splicing/Repair", "Detail": "Splicing Factor (Proline/Glutamine Rich)"},
    "TAF1": {"Class": "Transcription", "Detail": "TFIID Subunit (H4K16ac Reader)"},
    "SMC5": {"Class": "Structural Maintenance", "Detail": "SMC5/6 Complex (HR/Replication)"},
    "INIP": {"Class": "Structural Maintenance", "Detail": "Soss Complex subunit"},

    # --- Less Specific / Other ---
    "RPS3": {"Class": "Other/Repair-Associated", "Detail": "Ribosomal Protein (DNA glycosylase activity reported)"},
    "ATAD3A": {"Class": "Mitochondrial", "Detail": "Mitochondrial DNA organization"},
    "ERCC6L2": {"Class": "Repair-Associated", "Detail": "Snf2-family helicase"},
    "LMNA": {"Class": "Nuclear Envelope", "Detail": "Lamin A/C (Genome Stability)"},
    "MACROD1": {"Class": "ADP-Ribosylation", "Detail": "Removes ADP-ribose"},
    "MAPT": {"Class": "Cytoskeleton", "Detail": "Microtubule-associated protein"},
    "NUCKS1": {"Class": "Chromatin", "Detail": "Nuclear Casein Kinase Substrate"},
    "PARK7": {"Class": "Oxidative Stress", "Detail": "DJ-1 (Antioxidant)"},
    "PPP1CA": {"Class": "Signaling", "Detail": "Protein Phosphatase 1"},
    "SFN": {"Class": "Signaling", "Detail": "14-3-3 Sigma (Cell Cycle Checkpoint)"},
    "MEN1": {"Class": "Transcription", "Detail": "Menin (Genome Stability)"},
    "HTATSF1": {"Class": "Transcription/Splicing", "Detail": "HIV Tat Specific Factor 1"},
    "MTREX": {"Class": "RNA Processing", "Detail": "Exosome complex"},
    "SNW1": {"Class": "Splicing/Transcription", "Detail": "Ski-interacting protein"},
    "NIPBL": {"Class": "Cohesin Loading", "Detail": "Cohesin loading factor"},
    "POGZ": {"Class": "Chromatin", "Detail": "Zinc finger protein (Mitosis/DDR)"},
    "RBBP6": {"Class": "Regulation", "Detail": "P53-interacting protein"},
    "CDC5L": {"Class": "Splicing/Repair", "Detail": "Pre-mRNA splicing factor"},
    "NABP2": {"Class": "SSB Binding", "Detail": "Soss Complex (SSB1)"},
    "CHCHD6": {"Class": "Mitochondrial", "Detail": "Mitochondrial protein"},
    "ZMYND8": {"Class": "Chromatin", "Detail": "Transcriptional Repressor (DDR recruited)"}
}

# 3. Create data list
data = []
for sym in symbols:
    info = classification_map.get(sym, {"Class": "Unclassified/Other", "Detail": "Check literature"})
    data.append({
        "Symbol": sym,
        "Primary_Classification": info["Class"],
        "Functional_Detail": info["Detail"]
    })

# 4. Create DataFrame
df = pd.DataFrame(data)

# 5. Save
csv_filename = "lac_ddr_final_classification.csv"
df.to_csv(csv_filename, index=False)

print(f"CSV file created: {csv_filename}")
print(df.head(10))