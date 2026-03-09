import pandas as pd
import ast
import re
import numpy as np
import umap

INPUT_FILE = "/Users/gzy2520/Desktop/lac/filter/symbol_and_pathway.csv"
OUTPUT_FOR_R = "/Users/gzy2520/Desktop/lac/filter/umap_8type_data.csv"

EXACT_MAP = {
    "同源重组（HR）": {
        "GOBP_HOMOLOGOUS_RECOMBINATION", "KEGG_HOMOLOGOUS_RECOMBINATION",
        "REACTOME_HOMOLOGY_DIRECTED_REPAIR", "GOBP_RECOMBINATIONAL_REPAIR",
        "REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE", "BIOCARTA_ATRBRCA_PATHWAY",
        "GOBP_INTERSTRAND_CROSS_LINK_REPAIR", "BAE_BRCA1_TARGETS_UP",
        "BAE_BRCA1_TARGETS_DN", "PUJANA_BRCA1_PCC_NETWORK", "PUJANA_BRCA2_PCC_NETWORK"
    },
    "非同源末端连接（NHEJ）": {
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_CLASSICAL_NONHOMOLOGOUS_END_JOINING",
        "KEGG_NON_HOMOLOGOUS_END_JOINING", "GOCC_DNA_DEPENDENT_PROTEIN_KINASE_COMPLEX",
        "GOCC_NONHOMOLOGOUS_END_JOINING_COMPLEX", "GNF2_XRCC5", "MORF_XRCC5"
    },
    "碱基切除修复（BER）": {
        "GOBP_BASE_EXCISION_REPAIR", "KEGG_BASE_EXCISION_REPAIR",
        "REACTOME_BASE_EXCISION_REPAIR", "WP_BASE_EXCISION_REPAIR",
        "GOBP_DNA_ALKYLATION_REPAIR", "GCM_APEX1", "GNF2_APEX1"
    },
    "核苷酸切除修复（NER）": {
        "GOBP_NUCLEOTIDE_EXCISION_REPAIR", "KEGG_NUCLEOTIDE_EXCISION_REPAIR",
        "REACTOME_NUCLEOTIDE_EXCISION_REPAIR", "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
        "GOBP_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR", "DACOSTA_UV_RESPONSE_VIA_ERCC3_DN",
        "ENK_UV_RESPONSE_KERATINOCYTE_DN", "WP_NUCLEOTIDE_EXCISION_REPAIR_IN_XERODERMA_PIGMENTOSUM",
        "GCM_ERCC4", "MORF_DDB1"
    },
    "错配修复（MMR）": {
        "GOBP_MISMATCH_REPAIR", "KEGG_MISMATCH_REPAIR",
        "REACTOME_MISMATCH_REPAIR", "WP_DNA_MISMATCH_REPAIR",
        "WATANABE_COLON_CANCER_MSI_VS_MSS_DN"
    },
    "跨损伤合成（TLS）": {
        "GOBP_TRANSLESION_SYNTHESIS", "REACTOME_TRANSLESION_SYNTHESIS",
        "WP_TRANSLESION_SYNTHESIS", "GOBP_DNA_DAMAGE_BYPASS"
    },
    "直接修复（DRR）": {
        "GOBP_DNA_DIRECT_REPAIR", "REACTOME_DIRECT_REVERSAL_OF_DAMAGE",
        "GOBP_REVERSAL_OF_ALKYLATION_DAMAGE"
    },
    "损伤应答检查点（CP）": {
        "GOBP_DNA_INTEGRITY_CHECKPOINT_SIGNALING", "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
        "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS", "BIOCARTA_ATM_PATHWAY",
        "WP_DNA_IRDAMAGE_AND_CELLULAR_RESPONSE_VIA_ATR",
        "WP_DNA_IRDOUBLE_STRAND_BREAKS_AND_CELLULAR_RESPONSE_VIA_ATM",
        "KEGG_P53_SIGNALING_PATHWAY", "INGA_TP53_TARGETS", "PEREZ_TP53_TARGETS",
        "GOBP_APOPTOTIC_PROCESS", "GOBP_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY",
        "GOBP_DNA_DAMAGE_RESPONSE", "GOBP_DNA_REPAIR", "GOBP_DOUBLE_STRAND_BREAK_REPAIR",
        "REACTOME_DNA_REPAIR", "KAUFFMANN_DNA_REPAIR_GENES",
        "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR", "GOBP_REGULATION_OF_DNA_REPAIR",
        "GOCC_SITE_OF_DNA_DAMAGE", "GOCC_SITE_OF_DOUBLE_STRAND_BREAK",
        "GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE"
    }
}

FUZZY_CONFIG = {
    "同源重组（HR）": ["HOMOLOGOUS RECOMBINATION", "BRCA1", "BRCA2", "RAD51"],
    "非同源末端连接（NHEJ）": ["NON-HOMOLOGOUS", "NHEJ", "XRCC5", "DNA-PK", "LIG4"],
    "碱基切除修复（BER）": ["BASE EXCISION", "BER", "PARP", "APEX", "XRCC1"],
    "核苷酸切除修复（NER）": ["NUCLEOTIDE EXCISION", "NER", "XPC", "ERCC", "XPA"],
    "错配修复（MMR）": ["MISMATCH REPAIR", "MMR", "MSH2", "MSH6", "MLH1"],
    "跨损伤合成（TLS）": ["TRANSLESION", "TLS", "DNA POLYMERASE ZETA", "REV1"],
    "直接修复（DRR）": ["DIRECT REPAIR", "MGMT", "ALKBH"],
    "损伤应答检查点（CP）": ["CHECKPOINT", "ATM", "ATR", "TP53", "P53", "DNA DAMAGE RESPONSE"]
}

# 中英文映射
NAME_MAPPING = {
    "同源重组（HR）": "HR", "非同源末端连接（NHEJ）": "NHEJ",
    "碱基切除修复（BER）": "BER", "核苷酸切除修复（NER）": "NER",
    "错配修复（MMR）": "MMR", "跨损伤合成（TLS）": "TLS",
    "直接修复（DRR）": "DRR", "损伤应答检查点（CP）": "CP"
}


def process_data_8type():
    print(f"读取: {INPUT_FILE}")
    try:
        df = pd.read_csv(INPUT_FILE, encoding="utf-8-sig")
    except:
        df = pd.read_csv(INPUT_FILE, encoding="gbk")

    # 1. 预编译正则和集合
    EXACT_MAP_UPPER = {k: {x.upper() for x in v} for k, v in EXACT_MAP.items()}
    REGEX_MAP = {k: re.compile("|".join(map(re.escape, v)), re.IGNORECASE) for k, v in FUZZY_CONFIG.items()}

    def parse_pathways(val):
        try:
            return [x.strip().upper() for x in ast.literal_eval(val)]
        except:
            s = str(val).strip().upper()
            return [s] if s and s != 'NAN' else []

    df["PATHWAY_LIST"] = df["STANDARD_NAME"].apply(parse_pathways)

    # 2. 计数核心函数
    def get_counts(pathways):
        # 初始化 8 个分类的计数
        counts = {k: 0 for k in EXACT_MAP.keys()}

        for p in pathways:
            cat_found = None
            # A. 精确匹配
            for cat, p_set in EXACT_MAP_UPPER.items():
                if p in p_set:
                    cat_found = cat
                    break
            # B. 模糊匹配 (如果精确没找到)
            if not cat_found:
                for cat, regex in REGEX_MAP.items():
                    if regex.search(p):
                        cat_found = cat
                        break

            if cat_found:
                counts[cat_found] += 1

        return counts

    print("正在根据 8 个分类计算...")
    count_data = df["PATHWAY_LIST"].apply(get_counts).tolist()

    # 转为 DataFrame 并重命名列为英文简写
    count_df = pd.DataFrame(count_data)
    count_df.rename(columns=NAME_MAPPING, inplace=True)
    categories = list(NAME_MAPPING.values())  # ['HR', 'NHEJ', 'BER', 'NER', 'MMR', 'TLS', 'DRR', 'CP']

    # 3. 计算比例
    row_sums = count_df.sum(axis=1)
    prop_df = count_df.div(row_sums.replace(0, 1), axis=0)  # 防止除以0

    # 4. 生成 UMAP
    print("正在生成 UMAP ")
    feature_matrix = prop_df.values
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=25, n_jobs=1)
    embedding = reducer.fit_transform(feature_matrix)

    # 5. 合并结果
    df_out = pd.concat([df, prop_df], axis=1)
    df_out["UMAP_1"] = embedding[:, 0]
    df_out["UMAP_2"] = embedding[:, 1]
    df_out["Total_Hits"] = row_sums

    # 标记主类 (用于标签)
    def label_group(row):
        if row["Total_Hits"] == 0: return "None"
        # 值大于0的列
        active = [col for col in categories if row[col] > 0]
        if len(active) == 1:
            return active[0]
        else:
            return "Hybrid"

    df_out["Group"] = df_out.apply(label_group, axis=1)

    # 导出
    cols = ["Symbol", "UMAP_1", "UMAP_2", "Group", "Total_Hits"] + categories
    df_out[cols].to_csv(OUTPUT_FOR_R, index=False, encoding="utf-8-sig")
    print(f"处理完成！文件已保存至: {OUTPUT_FOR_R}")


if __name__ == "__main__":
    process_data_8type()