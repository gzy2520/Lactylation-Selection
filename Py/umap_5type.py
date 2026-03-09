import pandas as pd
import ast
import numpy as np
import umap

# ================= 配置区 =================
# 输入文件路径
INPUT_FILE = "/Users/gzy2520/Desktop/lac/filter/symbol_and_pathway.csv"
# 输出给 R 语言用的文件路径
OUTPUT_FOR_R = "/Users/gzy2520/Desktop/lac/filter/umap_5type_data.csv"

# 你的精确匹配字典
EXACT_MAP = {
    "NER": {  # 简化列名，方便 R 处理
        "DACOSTA_UV_RESPONSE_VIA_ERCC3_COMMON_DN", "DACOSTA_UV_RESPONSE_VIA_ERCC3_DN",
        "DACOSTA_UV_RESPONSE_VIA_ERCC3_TTD_UP", "DACOSTA_UV_RESPONSE_VIA_ERCC3_UP",
        "DAZARD_RESPONSE_TO_UV_NHEK_DN", "DAZARD_RESPONSE_TO_UV_NHEK_UP",
        "DAZARD_RESPONSE_TO_UV_SCC_DN", "DAZARD_UV_RESPONSE_CLUSTER_G5",
        "DAZARD_UV_RESPONSE_CLUSTER_G6", "ENK_UV_RESPONSE_EPIDERMIS_DN",
        "ENK_UV_RESPONSE_EPIDERMIS_UP", "ENK_UV_RESPONSE_KERATINOCYTE_DN",
        "ENK_UV_RESPONSE_KERATINOCYTE_UP", "GCM_ERCC4", "GENTILE_UV_HIGH_DOSE_DN",
        "GENTILE_UV_LOW_DOSE_DN", "GENTILE_UV_RESPONSE_CLUSTER_D4",
        "GENTILE_UV_RESPONSE_CLUSTER_D5", "GENTILE_UV_RESPONSE_CLUSTER_D6",
        "GENTILE_UV_RESPONSE_CLUSTER_D8", "GOBP_NUCLEOTIDE_EXCISION_REPAIR",
        "GOBP_PYRIMIDINE_DIMER_REPAIR", "GOBP_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR",
        "GOBP_UV_DAMAGE_EXCISION_REPAIR", "HP_DEFECTIVE_DNA_REPAIR_AFTER_ULTRAVIOLET_RADIATION_DAMAGE",
        "KEGG_NUCLEOTIDE_EXCISION_REPAIR", "MURAKAMI_UV_RESPONSE_6HR_UP",
        "REACTOME_DNA_DAMAGE_RECOGNITION_IN_GG_NER",
        "REACTOME_FORMATION_OF_TRANSCRIPTION_COUPLED_NER_TC_NER_REPAIR_COMPLEX",
        "REACTOME_GAP_FILLING_DNA_REPAIR_SYNTHESIS_AND_LIGATION_IN_GG_NER",
        "REACTOME_GLOBAL_GENOME_NUCLEOTIDE_EXCISION_REPAIR_GG_NER",
        "REACTOME_NUCLEOTIDE_EXCISION_REPAIR",
        "REACTOME_REPAIR_SYNTHESIS_FOR_GAP_FILLING_BY_DNA_POL_IN_TC_NER",
        "REACTOME_TRANSCRIPTION_COUPLED_NUCLEOTIDE_EXCISION_REPAIR_TC_NER",
        "WP_NUCLEOTIDE_EXCISION_REPAIR", "WP_NUCLEOTIDE_EXCISION_REPAIR_IN_XERODERMA_PIGMENTOSUM"
    },
    "HR": {
        "BAE_BRCA1_TARGETS_DN", "BAE_BRCA1_TARGETS_UP", "BIOCARTA_ATRBRCA_PATHWAY",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_SINGLE_STRAND_ANNEALING",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_SYNTHESIS_DEPENDENT_STRAND_ANNEALING",
        "GOBP_HOMOLOGOUS_RECOMBINATION", "GOBP_MITOTIC_RECOMBINATION",
        "GOBP_NEGATIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
        "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
        "GOBP_POSITIVE_REGULATION_OF_DNA_RECOMBINATION", "GOBP_RECOMBINATIONAL_REPAIR",
        "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_HOMOLOGOUS_RECOMBINATION",
        "GOBP_REPLICATION_BORN_DOUBLE_STRAND_BREAK_REPAIR_VIA_SISTER_CHROMATID_EXCHANGE",
        "HONRADO_BREAST_CANCER_BRCA1_VS_BRCA2", "KEGG_HOMOLOGOUS_RECOMBINATION",
        "MACLACHLAN_BRCA1_TARGETS_UP", "REACTOME_HOMOLOGOUS_DNA_PAIRING_AND_STRAND_EXCHANGE",
        "REACTOME_HOMOLOGOUS_RECOMBINATION_REPAIR_OF_REPLICATION_INDEPENDENT_DOUBLE_STRAND_BREAKS",
        "REACTOME_HOMOLOGY_DIRECTED_REPAIR", "VANTVEER_BREAST_CANCER_BRCA1_UP"
    },
    "NHEJ": {
        "COLLIS_PRKDC_REGULATORS", "COLLIS_PRKDC_SUBSTRATES", "GNF2_XRCC5",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_ALTERNATIVE_NONHOMOLOGOUS_END_JOINING",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_CLASSICAL_NONHOMOLOGOUS_END_JOINING",
        "GOBP_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
        "GOBP_NEGATIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
        "GOBP_POSITIVE_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
        "GOBP_REGULATION_OF_DOUBLE_STRAND_BREAK_REPAIR_VIA_NONHOMOLOGOUS_END_JOINING",
        "GOCC_DNA_DEPENDENT_PROTEIN_KINASE_COMPLEX",
        "GOCC_DNA_DEPENDENT_PROTEIN_KINASE_DNA_LIGASE_4_COMPLEX",
        "GOCC_NONHOMOLOGOUS_END_JOINING_COMPLEX", "KEGG_NON_HOMOLOGOUS_END_JOINING",
        "MORF_XRCC5"
    },
    "BER": {
        "GCM_APEX1", "GNF2_APEX1", "GOBP_BASE_EXCISION_REPAIR",
        "GOBP_BASE_EXCISION_REPAIR_GAP_FILLING", "GOBP_DNA_ALKYLATION_REPAIR",
        "GOMF_DNA_APURINIC_OR_APYRIMIDINIC_SITE_ENDONUCLEASE_ACTIVITY",
        "GOMF_DNA_N_GLYCOSYLASE_ACTIVITY", "KEGG_BASE_EXCISION_REPAIR",
        "REACTOME_BASE_EXCISION_REPAIR",
        "REACTOME_PCNA_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR",
        "REACTOME_POLB_DEPENDENT_LONG_PATCH_BASE_EXCISION_REPAIR",
        "REACTOME_RECOGNITION_AND_ASSOCIATION_OF_DNA_GLYCOSYLASE_WITH_SITE_CONTAINING_AN_AFFECTED_PURINE",
        "WP_BASE_EXCISION_REPAIR"
    },
    "Others": {
        "BIOCARTA_ATM_PATHWAY", "BIOCARTA_G1_PATHWAY", "BIOCARTA_G2_PATHWAY",
        "BIOCARTA_P53HYPOXIA_PATHWAY", "DNA_DAMAGE_CHECKPOINT",
        "GOBP_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
        "GOBP_DNA_INTEGRITY_CHECKPOINT_SIGNALING",
        "GOBP_DNA_REPLICATION_CHECKPOINT_SIGNALING",
        "GOBP_MITOTIC_DNA_INTEGRITY_CHECKPOINT_SIGNALING",
        "GOBP_MITOTIC_INTRA_S_DNA_DAMAGE_CHECKPOINT_SIGNALING",
        "GOBP_POSITIVE_REGULATION_OF_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
        "GOBP_REGULATION_OF_DNA_DAMAGE_CHECKPOINT",
        "GOBP_REGULATION_OF_DNA_DAMAGE_RESPONSE_SIGNAL_TRANSDUCTION_BY_P53_CLASS_MEDIATOR",
        "GOBP_SIGNAL_TRANSDUCTION_IN_RESPONSE_TO_DNA_DAMAGE", "INGA_TP53_TARGETS",
        "KANNAN_TP53_TARGETS_UP", "KEGG_CELL_CYCLE", "KEGG_P53_SIGNALING_PATHWAY",
        "PEREZ_TP53_TARGETS", "REACTOME_G1_S_DNA_DAMAGE_CHECKPOINTS",
        "REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT",
        "REACTOME_P53_INDEPENDENT_G1_S_DNA_DAMAGE_CHECKPOINT",
        "SCHAVOLT_TARGETS_OF_TP53_AND_TP63",
        "SCIAN_CELL_CYCLE_TARGETS_OF_TP53_AND_TP73_UP",
        "TANG_SENESCENCE_TP53_TARGETS_DN", "WHITFIELD_CELL_CYCLE_G2",
        "WHITFIELD_CELL_CYCLE_G2_M", "WHITFIELD_CELL_CYCLE_M_G1",
        "WP_DNA_DAMAGE_RESPONSE", "WP_DNA_IRDAMAGE_AND_CELLULAR_RESPONSE_VIA_ATR",
        "WP_DNA_IRDOUBLE_STRAND_BREAKS_AND_CELLULAR_RESPONSE_VIA_ATM",
        "ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_24HR",
        "ZHOU_CELL_CYCLE_GENES_IN_IR_RESPONSE_6HR"
    }
}

FUZZY_KEYWORDS = {
    "HR": ["HOMOLOGOUS RECOMBINATION", "HR", "BRCA", "RAD51", "RAD52", "DSB REPAIR"],
    "BER": ["BASE EXCISION REPAIR", "BER", "PARP", "APEX", "XRCC1"],
    "NHEJ": ["NON-HOMOLOGOUS END JOINING", "NHEJ", "XRCC5", "XRCC6", "PRKDC", "DNA-PK"],
    "Others": ["CHECKPOINT", "ATM", "ATR", "CHEK", "TP53", "P53", "CELL CYCLE ARREST"],
    "NER": ["NUCLEOTIDE EXCISION REPAIR", "NER", "ERCC", "XPC", "DDB2", "UV REPAIR"]
}


def process_data():
    print(f"读取: {INPUT_FILE}")
    try:
        df = pd.read_csv(INPUT_FILE, encoding="utf-8-sig")
    except:
        df = pd.read_csv(INPUT_FILE, encoding="gbk")

    # 1. 预处理
    EXACT_MAP_UPPER = {k: {x.upper() for x in v} for k, v in EXACT_MAP.items()}

    def parse_pathways(val):
        try:
            return [x.strip().upper() for x in ast.literal_eval(val)]
        except:
            s = str(val).strip().upper()
            return [s] if s and s != 'NAN' else []

    df["PATHWAY_LIST"] = df["STANDARD_NAME"].apply(parse_pathways)

    # 2. 计数
    def get_counts(pathways):
        # 注意：这里 keys 的顺序很重要，要和下面 DataFrame 建立时一致
        counts = {k: 0 for k in EXACT_MAP.keys()}
        for p in pathways:
            matched = False
            for cat, p_set in EXACT_MAP_UPPER.items():
                if p in p_set:
                    counts[cat] += 1
                    matched = True
                    break
            if not matched:
                for cat, keywords in FUZZY_KEYWORDS.items():
                    if any(kw in p for kw in keywords):
                        counts[cat] += 1
                        matched = True
                        break
        return counts

    # 计算每个基因在 5 个分类的计数
    count_data = df["PATHWAY_LIST"].apply(get_counts).tolist()
    categories = list(EXACT_MAP.keys())  # ['NER', 'HR', 'NHEJ', 'BER', 'Others']

    # 转换为 DataFrame
    count_df = pd.DataFrame(count_data, columns=categories)

    # 3. 计算比例 (用于画饼图)
    # 增加一个 epsilon 防止全0行除法报错，或者稍后处理
    row_sums = count_df.sum(axis=1)

    # 我们只保留有至少命中1个通路的基因，或者你要保留全部？
    # 建议保留全部，未命中的在图上显示为灰色点

    # 计算比例 (Proportion)
    # 如果总和是0，比例全设为0
    prop_df = count_df.div(row_sums.replace(0, 1), axis=0)

    # 将比例列合并回主表
    df_out = pd.concat([df, prop_df], axis=1)
    df_out["Total_Hits"] = row_sums

    # 4. 生成 UMAP 坐标
    # 使用比例矩阵作为特征，这样功能相似的基因（比如都是 50% HR + 50% NHEJ）会聚在一起
    print("正在计算 UMAP 坐标...")
    feature_matrix = prop_df.values

    # 处理全0行（没有通路的基因），加极微小噪音防止报错，或者让它们聚在一起
    mask_zero = row_sums == 0
    if mask_zero.any():
        feature_matrix[mask_zero] = np.random.uniform(0, 0.001, (mask_zero.sum(), 5))

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=25, n_jobs=1)
    embedding = reducer.fit_transform(feature_matrix)

    df_out["UMAP_1"] = embedding[:, 0]
    df_out["UMAP_2"] = embedding[:, 1]

    # 5. 标记分类标签 (用于在 R 里筛选标签)
    def label_text(row):
        if row["Total_Hits"] == 0: return "None"
        # 找出比例大于0的列
        active = [col for col in categories if row[col] > 0]
        if len(active) == 1:
            return active[0]
        else:
            return "Hybrid"

    df_out["Group"] = df_out.apply(label_text, axis=1)

    # 6. 导出
    # 我们只需要 Symbol, 5个分类的比例, UMAP坐标, Group, Total_Hits
    cols_to_export = ["Symbol", "UMAP_1", "UMAP_2", "Group", "Total_Hits"] + categories
    df_out[cols_to_export].to_csv(OUTPUT_FOR_R, index=False)
    print(f"处理完成！文件已保存至: {OUTPUT_FOR_R}")
    print("现在请打开 RStudio 运行下一步的代码。")


if __name__ == "__main__":
    process_data()