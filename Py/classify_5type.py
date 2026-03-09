import pandas as pd
import ast

# ==========================================
# 1. 配置区域
# ==========================================

# --- A. 精确匹配字典 (优先级最高) ---
# 这里直接使用你整理好的高质量 Set
EXACT_MAP = {
    "核苷酸切除修复（NER）": {
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
    "同源重组修复（HR）": {
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
    "非同源末端连接（NHEJ）": {
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
    "碱基切除修复（BER）": {
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
    "其他（Others）": {
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

# --- B. 模糊匹配关键词 (优先级低，用于兜底) ---
FUZZY_KEYWORDS = {
    "同源重组修复（HR）": ["HOMOLOGOUS RECOMBINATION", "HR", "BRCA", "RAD51", "RAD52", "DSB REPAIR"],
    "碱基切除修复（BER）": ["BASE EXCISION REPAIR", "BER", "PARP", "APEX", "XRCC1"],
    "非同源末端连接（NHEJ）": ["NON-HOMOLOGOUS END JOINING", "NHEJ", "XRCC5", "XRCC6", "PRKDC", "DNA-PK"],
    "其他（Others）": ["CHECKPOINT", "ATM", "ATR", "CHEK", "TP53", "P53", "CELL CYCLE ARREST"],
    "核苷酸切除修复（NER）": ["NUCLEOTIDE EXCISION REPAIR", "NER", "ERCC", "XPC", "DDB2", "UV REPAIR"]
}


# ==========================================
# 主处理函数
# ==========================================

def run_hybrid_classification(input_csv, output_csv):
    print(f"正在读取文件: {input_csv}")
    try:
        df = pd.read_csv(input_csv, encoding="utf-8-sig")
    except:
        df = pd.read_csv(input_csv, encoding="gbk")

    # 1. 预处理：将 Exact Map 转为全大写 Set，确保匹配无误
    EXACT_MAP_UPPER = {k: {x.upper() for x in v} for k, v in EXACT_MAP.items()}

    # 记录哪些通路是靠关键词“救”回来的（用于人工核查）
    rescued_pathways = []

    def parse_pathways(val):
        """解析字符串列表"""
        try:
            return [x.strip().upper() for x in ast.literal_eval(val)]
        except:
            s = str(val).strip().upper()
            return [s] if s and s != 'NAN' else []

    df["PATHWAY_LIST"] = df["STANDARD_NAME"].apply(parse_pathways)

    # 2. 核心分类逻辑
    def count_hybrid(pathways):
        counts = {k: 0 for k in EXACT_MAP.keys()}
        counts["未分类"] = 0

        for p in pathways:
            matched = False

            # A. 第一轮：精确匹配 (优先)
            for cat, p_set in EXACT_MAP_UPPER.items():
                if p in p_set:
                    counts[cat] += 1
                    matched = True
                    break

            # B. 第二轮：模糊匹配 (兜底)
            if not matched:
                for cat, keywords in FUZZY_KEYWORDS.items():
                    # 检查该通路是否包含任意关键词
                    if any(kw in p for kw in keywords):
                        counts[cat] += 1
                        matched = True
                        # 记录一下，这是靠关键词匹配上的
                        rescued_pathways.append(f"[{cat}] {p}")
                        break  # 假设一个通路只归一类，避免重复计数

            # C. 还没匹配上 -> 未分类
            if not matched:
                counts["未分类"] += 1

        return counts

    print("正在执行混合匹配（精确优先 + 关键词兜底）...")
    df["COUNTS"] = df["PATHWAY_LIST"].apply(count_hybrid)

    # 3. 展开结果列
    for cat in EXACT_MAP.keys():
        df[f"重叠数量_{cat}"] = df["COUNTS"].apply(lambda x: x[cat])
    df["重叠数量_未分类"] = df["COUNTS"].apply(lambda x: x["未分类"])

    # 4. 判定归属
    def assign_category(c_dict):
        # 排除未分类，看谁最高
        valid = {k: v for k, v in c_dict.items() if k != "未分类"}
        max_v = max(valid.values()) if valid else 0

        if max_v == 0: return "未归类"

        tops = [k for k, v in valid.items() if v == max_v]
        if len(tops) == 1:
            return tops[0]
        else:
            # 简写混合类
            mapping = {"同源重组修复（HR）": "HR", "碱基切除修复（BER）": "BER",
                       "非同源末端连接（NHEJ）": "NHEJ", "其他（Others）": "Others",
                       "核苷酸切除修复（NER）": "NER"}
            names = [mapping.get(t, t) for t in tops]
            return f"混合类（{'、'.join(names)}）"

    df["GENE_BELONG_CATEGORY"] = df["COUNTS"].apply(assign_category)

    # 5. 保存与报告
    cols = ["Symbol", "GENE_BELONG_CATEGORY"] + \
           [f"重叠数量_{k}" for k in EXACT_MAP.keys()] + \
           ["重叠数量_未分类", "STANDARD_NAME"]

    df[cols].to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"结果已保存至: {output_csv}")

    # 打印统计
    print("\n=== 最终分类统计 ===")
    print(df["GENE_BELONG_CATEGORY"].value_counts())

    # 打印兜底情况
    unique_rescued = list(set(rescued_pathways))
    if unique_rescued:
        print("\n" + "=" * 60)
        print(f"【提示】有 {len(unique_rescued)} 个通路不在精确列表中，但通过关键词被成功召回：")
        print("(建议检查这些通路，如果确实相关，以后可以加到精确列表中)")
        for rp in unique_rescued[:10]:  # 只打印前10个
            print(f"  -> {rp}")
        if len(unique_rescued) > 10: print(f"  ... 以及其他 {len(unique_rescued) - 10} 个")
        print("=" * 60 + "\n")
    else:
        print("\n所有匹配均来自精确列表，未触发模糊匹配。")


# ==========================================
# 执行
# ==========================================
if __name__ == "__main__":
    INPUT_FILE = "/Users/gzy2520/Desktop/lac/symbol_and_pathway.csv"
    OUTPUT_FILE = "/Users/gzy2520/Desktop/lac/classified_symbol.csv"

    run_hybrid_classification(INPUT_FILE, OUTPUT_FILE)