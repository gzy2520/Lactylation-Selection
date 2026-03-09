import pandas as pd
import ast
import re
from pathlib import Path



EXACT_MAP = {
    # --- [1] 同源重组 (HR) ---
    "同源重组（HR）": {
        "GOBP_HOMOLOGOUS_RECOMBINATION",
        "KEGG_HOMOLOGOUS_RECOMBINATION",
        "REACTOME_HOMOLOGY_DIRECTED_REPAIR",
        "GOBP_RECOMBINATIONAL_REPAIR",
        "REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE",
        "BIOCARTA_ATRBRCA_PATHWAY",
        "GOBP_INTERSTRAND_CROSS_LINK_REPAIR",
        "BAE_BRCA1_TARGETS_UP",
        "BAE_BRCA1_TARGETS_DN",
        "PUJANA_BRCA1_PCC_NETWORK",
        "PUJANA_BRCA2_PCC_NETWORK"
    },

    # --- [2] 非同源末端连接 (NHEJ) ---
    "非同源末端连接（NHEJ）": {
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_CLASSICAL_NONHOMOLOGOUS_END_JOINING",
        "KEGG_NON_HOMOLOGOUS_END_JOINING",
        "GOCC_DNA_DEPENDENT_PROTEIN_KINASE_COMPLEX",
        "GOCC_NONHOMOLOGOUS_END_JOINING_COMPLEX",
        "GNF2_XRCC5",
        "MORF_XRCC5"
    },

    # --- [3] 碱基切除修复 (BER) ---
    "碱基切除修复（BER）": {
        "GOBP_BASE_EXCISION_REPAIR",
        "KEGG_BASE_EXCISION_REPAIR",
        "REACTOME_BASE_EXCISION_REPAIR",
        "WP_BASE_EXCISION_REPAIR",
        "GOBP_DNA_ALKYLATION_REPAIR",
        "GCM_APEX1",
        "GNF2_APEX1"
    },

    # --- [4] 核苷酸切除修复 (NER) ---
    "核苷酸切除修复（NER）": {
        "GOBP_NUCLEOTIDE_EXCISION_REPAIR",
        "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
        "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
        "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
        "GOBP_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR",
        "DACOSTA_UV_RESPONSE_VIA_ERCC3_DN",
        "ENK_UV_RESPONSE_KERATINOCYTE_DN",
        "WP_NUCLEOTIDE_EXCISION_REPAIR_IN_XERODERMA_PIGMENTOSUM",
        "GCM_ERCC4",
        "MORF_DDB1"
    },

    # --- [5] 错配修复 (MMR) ---
    "错配修复（MMR）": {
        "GOBP_MISMATCH_REPAIR",
        "KEGG_MISMATCH_REPAIR",
        "REACTOME_MISMATCH_REPAIR",
        "WP_DNA_MISMATCH_REPAIR",
        "WATANABE_COLON_CANCER_MSI_VS_MSS_DN"
    },

    # --- [6] 跨损伤合成 (TLS) ---
    "跨损伤合成（TLS）": {
        "GOBP_TRANSLESION_SYNTHESIS",
        "REACTOME_TRANSLESION_SYNTHESIS",
        "WP_TRANSLESION_SYNTHESIS",
        "GOBP_DNA_DAMAGE_BYPASS"
    },

    # --- [7] 直接修复 (DRR) ---
    "直接修复（DRR）": {
        "GOBP_DNA_DIRECT_REPAIR",
        "REACTOME_DIRECT_REVERSAL_OF_DAMAGE",
        "GOBP_REVERSAL_OF_ALKYLATION_DAMAGE"
    },

    # --- [8] 损伤应答检查点 (CP) & 通用DDR ---
    "损伤应答检查点（CP）": {
        "GOBP_DNA_INTEGRITY_CHECKPOINT_SIGNALING",
        "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
        "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS",
        "BIOCARTA_ATM_PATHWAY",
        "WP_DNA_IRDAMAGE_AND_CELLULAR_RESPONSE_VIA_ATR",
        "WP_DNA_IRDOUBLE_STRAND_BREAKS_AND_CELLULAR_RESPONSE_VIA_ATM",

        # p53 相关
        "KEGG_P53_SIGNALING_PATHWAY",
        "INGA_TP53_TARGETS",
        "PEREZ_TP53_TARGETS",

        # 凋亡
        "GOBP_APOPTOTIC_PROCESS",
        "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",

        # 通用
        "GOBP_DNA_DAMAGE_RESPONSE",
        "GOBP_DNA_REPAIR",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR",
        "REACTOME_DNA_REPAIR",
        "KAUFFMANN_DNA_REPAIR_GENES",
        "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
        "GOBP_REGULATION_OF_DNA_REPAIR",
        "GOCC_SITE_OF_DNA_DAMAGE",
        "GOCC_SITE_OF_DOUBLE_STRAND_BREAK",
        "GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE"
    }
}

# 2. 模糊匹配
FUZZY_CONFIG = {
    "同源重组（HR）": ["HOMOLOGOUS RECOMBINATION", "BRCA1", "BRCA2", "RAD51"],
    "非同源末端连接（NHEJ）": ["NON-HOMOLOGOUS", "NHEJ", "XRCC5", "DNA-PK", "LIG4"],
    "碱基切除修复（BER）": ["BASE EXCISION", "BER", "PARP", "APEX", "XRCC1"],
    "核苷酸切除修复（NER）": ["NUCLEOTIDE EXCISION", "NER", "XPC", "ERCC", "XPA"],
    "损伤应答检查点（CP）": ["CHECKPOINT", "ATM", "ATR", "TP53", "P53", "DNA DAMAGE RESPONSE"]
}

SPECIFIC_CATS = {
    "同源重组（HR）", "非同源末端连接（NHEJ）", "碱基切除修复（BER）",
    "核苷酸切除修复（NER）", "错配修复（MMR）", "跨损伤合成（TLS）", "直接修复（DRR）"
}
GENERAL_CAT = "损伤应答检查点（CP）"


def build_lookup_tables():
    lookup_map = {}
    for cat, pathways in EXACT_MAP.items():
        for p in pathways:
            lookup_map[p.upper()] = cat

    regex_map = {}
    for cat, keywords in FUZZY_CONFIG.items():
        pattern_str = "|".join(map(re.escape, keywords))
        regex_map[cat] = re.compile(pattern_str, re.IGNORECASE)
    return lookup_map, regex_map


PATHWAY_LOOKUP, REGEX_MAP = build_lookup_tables()

def parse_pathways(val):
    if pd.isna(val): return []
    try:
        return [x.strip().upper() for x in ast.literal_eval(val)]
    except:
        s = str(val).strip().upper()
        return [s] if s and s != 'NAN' else []


def classify_gene_pathways(pathways):
    counts = {k: 0 for k in EXACT_MAP.keys()}
    counts["未分类"] = 0
    hit_categories = set()

    for p in pathways:
        cat_found = None
        # 1. 查表
        if p in PATHWAY_LOOKUP:
            cat_found = PATHWAY_LOOKUP[p]
        # 2. 正则
        else:
            for cat, regex in REGEX_MAP.items():
                if regex.search(p):
                    cat_found = cat
                    break

        if cat_found:
            counts[cat_found] += 1
            hit_categories.add(cat_found)
        else:
            counts["未分类"] += 1

    # 3. 优先级清洗
    is_specific = not hit_categories.isdisjoint(SPECIFIC_CATS)
    if is_specific and GENERAL_CAT in hit_categories:
        counts[GENERAL_CAT] = 0

    return counts


def assign_final_label(counts):
    valid_counts = {k: v for k, v in counts.items() if k != "未分类" and v > 0}
    if not valid_counts: return "未归类"

    sorted_cats = sorted(valid_counts.keys(), key=lambda k: valid_counts[k], reverse=True)
    if len(sorted_cats) == 1: return sorted_cats[0]

    short_names = {
        "同源重组（HR）": "HR", "非同源末端连接（NHEJ）": "NHEJ",
        "碱基切除修复（BER）": "BER", "核苷酸切除修复（NER）": "NER",
        "错配修复（MMR）": "MMR", "损伤应答检查点（CP）": "CP",
        "跨损伤合成（TLS）": "TLS", "直接修复（DRR）": "DRR"
    }
    labels = [short_names.get(c, c) for c in sorted_cats]
    return f"混合({'+'.join(labels)})"


def run_classification(input_csv, output_csv):
    print("正在读取数据...")
    try:
        df = pd.read_csv(input_csv, encoding="utf-8-sig")
    except:
        df = pd.read_csv(input_csv, encoding="gbk")

    print(f"读取 {len(df)} 行数据，正在处理...")

    df["parsed"] = df["STANDARD_NAME"].apply(parse_pathways)
    counts_df = df["parsed"].apply(classify_gene_pathways).apply(pd.Series)

    # 保存数值列
    for col in counts_df.columns:
        df[f"重叠数量_{col}"] = counts_df[col]

    # 生成最终标签 (重新根据 row 计算，确保逻辑一致)
    def get_label(row):
        c = {cat: row.get(f"重叠数量_{cat}", 0) for cat in EXACT_MAP.keys()}
        return assign_final_label(c)

    df["GENE_BELONG_CATEGORY"] = df.apply(get_label, axis=1)

    if "parsed" in df.columns: del df["parsed"]

    cols = ["Symbol", "GENE_BELONG_CATEGORY"] + \
           [c for c in df.columns if c.startswith("重叠数量_")] + \
           ["STANDARD_NAME"]

    final_cols = [c for c in cols if c in df.columns]
    df[final_cols].to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"完成！结果已保存至: {output_csv}")


if __name__ == "__main__":
    INPUT = Path("/Users/gzy2520/Desktop/lac/symbol_and_pathway.csv")
    OUTPUT = Path("/Users/gzy2520/Desktop/lac/classified_symbol.csv")

    run_classification(INPUT, OUTPUT)