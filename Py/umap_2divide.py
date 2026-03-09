import pandas as pd
import ast
import re
import numpy as np
import umap

# ================= 配置区 =================
# 输入文件（请修改为你自己的路径）
INPUT_FILE = "/Users/gzy2520/Desktop/lac/filter/symbol_and_pathway.csv"
OUTPUT_FOR_R = "/Users/gzy2520/Desktop/lac/filter/umap_2divide_data.csv"

# --- 1. 精细分类字典 (用于计算 UMAP 坐标) ---
# 我们把“通用DSB”单独拎出来，算作 DSB 的一部分
EXACT_MAP = {
    # === [DSB 阵营] ===
    "HR": {
        "GOBP_HOMOLOGOUS_RECOMBINATION", "KEGG_HOMOLOGOUS_RECOMBINATION",
        "REACTOME_HOMOLOGY_DIRECTED_REPAIR", "GOBP_RECOMBINATIONAL_REPAIR",
        "REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE", "BIOCARTA_ATRBRCA_PATHWAY",
        "BAE_BRCA1_TARGETS_UP", "BAE_BRCA1_TARGETS_DN", "PUJANA_BRCA1_PCC_NETWORK",
        "REACTOME_HOMOLOGOUS_RECOMBINATION_REPAIR_OF_REPLICATION_INDEPENDENT_DOUBLE_STRAND_BREAKS"
    },
    "NHEJ": {
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
        "KEGG_NON_HOMOLOGOUS_END_JOINING", "GOCC_DNA_DEPENDENT_PROTEIN_KINASE_COMPLEX",
        "GOCC_NONHOMOLOGOUS_END_JOINING_COMPLEX", "GNF2_XRCC5", "MORF_XRCC5",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_CLASSICAL_NONHOMOLOGOUS_END_JOINING"
    },
    "DSB_General": {  # 新增：通用双链断裂术语
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR", "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR",
        "REACTOME_DOUBLE_STRAND_BREAK_REPAIR", "WP_DNA_IRDOUBLE_STRAND_BREAKS_AND_CELLULAR_RESPONSE_VIA_ATM",
        "GOCC_SITE_OF_DOUBLE_STRAND_BREAK", "GOBP_PROTEIN_LOCALIZATION_TO_SITE_OF_DOUBLE_STRAND_BREAK"
    },

    # === [SSB 阵营] ===
    "BER": {
        "GOBP_BASE_EXCISION_REPAIR", "KEGG_BASE_EXCISION_REPAIR",
        "REACTOME_BASE_EXCISION_REPAIR", "WP_BASE_EXCISION_REPAIR", "GCM_APEX1"
    },
    "NER": {
        "GOBP_NUCLEOTIDE_EXCISION_REPAIR", "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
        "REACTOME_NUCLEOTIDE_EXCISION_REPAIR", "GOBP_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR",
        "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER", "WP_NUCLEOTIDE_EXCISION_REPAIR",
        "DACOSTA_UV_RESPONSE_VIA_ERCC3_DN", "GCM_ERCC4"
    },
    "MMR": {
        "GOBP_MISMATCH_REPAIR", "KEGG_MISMATCH_REPAIR", "REACTOME_MISMATCH_REPAIR", "WP_DNA_MISMATCH_REPAIR"
    },
    "TLS": {
        "GOBP_TRANSLESION_SYNTHESIS", "REACTOME_TRANSLESION_SYNTHESIS", "GOBP_DNA_DAMAGE_BYPASS"
    },
    "DRR": {
        "GOBP_DNA_DIRECT_REPAIR", "REACTOME_DIRECT_REVERSAL_OF_DAMAGE"
    },

    # === [Others/General 阵营] ===
    "Others": {
        "GOBP_DNA_INTEGRITY_CHECKPOINT_SIGNALING", "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
        "DNA_DAMAGE_CHECKPOINT", "BIOCARTA_ATM_PATHWAY", "KEGG_P53_SIGNALING_PATHWAY",
        "INGA_TP53_TARGETS", "PEREZ_TP53_TARGETS", "GOBP_DNA_DAMAGE_RESPONSE",
        "GOBP_DNA_REPAIR", "GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE"
    }
}

# --- 2. 模糊匹配 (用于兜底) ---
FUZZY_CONFIG = {
    "HR": ["HOMOLOGOUS RECOMBINATION", "BRCA1", "BRCA2", "RAD51"],
    "NHEJ": ["NON-HOMOLOGOUS", "NHEJ", "XRCC5", "DNA-PK"],
    "DSB_General": ["DOUBLE STRAND BREAK", "DSB REPAIR"],  # 优先匹配双链关键词
    "BER": ["BASE EXCISION", "BER", "PARP", "APEX", "XRCC1"],
    "NER": ["NUCLEOTIDE EXCISION", "NER", "XPC", "ERCC"],
    "MMR": ["MISMATCH REPAIR", "MMR", "MSH2", "MLH1"],
    "TLS": ["TRANSLESION", "TLS", "REV1"],
    "DRR": ["DIRECT REPAIR", "MGMT"],
    "Others": ["CHECKPOINT", "ATM", "ATR", "TP53", "P53", "DNA DAMAGE RESPONSE"]
}


def process_ssb_dsb():
    print(f"读取: {INPUT_FILE}")
    try:
        df = pd.read_csv(INPUT_FILE, encoding="utf-8-sig")
    except:
        df = pd.read_csv(INPUT_FILE, encoding="gbk")

    # 1. 解析 Pathway List
    EXACT_MAP_UPPER = {k: {x.upper() for x in v} for k, v in EXACT_MAP.items()}
    REGEX_MAP = {k: re.compile("|".join(map(re.escape, v)), re.IGNORECASE) for k, v in FUZZY_CONFIG.items()}

    def parse_pathways(val):
        try:
            return [x.strip().upper() for x in ast.literal_eval(val)]
        except:
            s = str(val).strip().upper()
            return [s] if s and s != 'NAN' else []

    df["PATHWAY_LIST"] = df["STANDARD_NAME"].apply(parse_pathways)

    # 2. 计算 9 个细分分类的权重
    def get_counts(pathways):
        counts = {k: 0 for k in EXACT_MAP.keys()}
        for p in pathways:
            cat_found = None
            for cat, p_set in EXACT_MAP_UPPER.items():
                if p in p_set:
                    cat_found = cat
                    break
            if not cat_found:
                for cat, regex in REGEX_MAP.items():
                    if regex.search(p):
                        cat_found = cat
                        break
            if cat_found:
                counts[cat_found] += 1
        return counts

    print("计算细分通路权重...")
    count_data = df["PATHWAY_LIST"].apply(get_counts).tolist()
    count_df = pd.DataFrame(count_data)

    # 3. 聚合为三大阵营 (用于画图颜色)
    # DSB阵营: HR + NHEJ + 通用DSB
    count_df["Macro_DSB"] = count_df["HR"] + count_df["NHEJ"] + count_df["DSB_General"]
    # SSB阵营: BER + NER + MMR + TLS + DRR
    count_df["Macro_SSB"] = count_df["BER"] + count_df["NER"] + count_df["MMR"] + count_df["TLS"] + count_df["DRR"]
    # Others阵营: Others
    count_df["Others"] = count_df["Others"]

    # 4. 计算 UMAP 坐标 (使用所有细分列，保证聚类精度)
    # 我们使用 9 个细分列进行 UMAP，这样 HR 和 NHEJ 依然会分开，而不是混成一团
    umap_features = count_df[list(EXACT_MAP.keys())]
    row_sums = umap_features.sum(axis=1)
    prop_df = umap_features.div(row_sums.replace(0, 1), axis=0)

    # 填充无数据的行
    mask_zero = row_sums == 0
    feature_matrix = prop_df.values
    if mask_zero.any():
        feature_matrix[mask_zero] = np.random.uniform(0, 0.001, (mask_zero.sum(), feature_matrix.shape[1]))

    print("生成 UMAP 坐标...")
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=25, n_jobs=1)
    embedding = reducer.fit_transform(feature_matrix)

    # 5. 准备输出数据
    # 计算用于画饼图的“宏观比例”
    macro_sums = count_df[["Macro_DSB", "Macro_SSB", "Others"]].sum(axis=1)
    macro_prop = count_df[["Macro_DSB", "Macro_SSB", "Others"]].div(macro_sums.replace(0, 1), axis=0)

    df_out = pd.concat([df, macro_prop], axis=1)
    df_out["UMAP_1"] = embedding[:, 0]
    df_out["UMAP_2"] = embedding[:, 1]
    df_out["Total_Hits"] = macro_sums

    # 导出
    cols = ["Symbol", "UMAP_1", "UMAP_2", "Total_Hits", "Macro_DSB", "Macro_SSB", "Others"]
    df_out[cols].to_csv(OUTPUT_FOR_R, index=False, encoding="utf-8-sig")
    print(f"完成！文件已保存至: {OUTPUT_FOR_R}")


if __name__ == "__main__":
    process_ssb_dsb()