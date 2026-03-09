import pandas as pd
import ast

KEYWORDS = {
    "DSB_Repair": [  # 双链断裂：同源重组 + NHEJ + Fanconi
        "HOMOLOGOUS", "RECOMBINATION", "NHEJ", "NON_HOMOLOGOUS", "NONHOMOLOGOUS",
        "DOUBLE_STRAND", "DSB", "BRCA", "RAD51", "FANCONI", "XRCC5", "DNA_PK",
        "PRKDC", "KU70", "KU80", "LIGASE_4", "END_JOINING", "INTERSTRAND"
    ],
    "SSB_Excision": [  # 单链/切除：BER + NER + MMR + UV
        "BASE_EXCISION", "BER", "NUCLEOTIDE_EXCISION", "NER", "MISMATCH", "MMR",
        "UV_RESPONSE", "RESPONSE_TO_UV", "ULTRAVIOLET", "SINGLE_STRAND", "EXCISION",
        "ERCC", "XPC", "DDB1", "APEX", "PARP", "XRCC1", "MSH2", "MSH6", "MLH1", "PMS2"
    ],
    "Replication": [  # DNA复制与压力
        "DNA_REPLICATION", "REPLISOME", "REPLICATION_FORK", "SYNTHESIS_OF_DNA",
        "STRAND_ELONGATION", "DNA_SYNTHESIS", "MCM_COMPLEX", "PCNA", "GINS"
    ],
    "Chromatin_Epi": [  # 染色质与表观
        "CHROMATIN", "HISTONE", "METHYLATION", "EPIGENETIC", "NUCLEOSOME",
        "HDAC", "SWI_SNF", "RPA", "SMC", "COHESIN", "CONDENSIN", "TOPOLOGICAL"
    ],
    "Checkpoint_p53": [  # 检查点与信号
        "CHECKPOINT", "ATM", "ATR", "CHK1", "CHK2", "TP53", "P53", "APOPTOSIS",
        "G1_PHASE", "G2_PHASE", "G2_M", "MITOTIC", "SPINDLE", "CELL_CYCLE"
    ],
    "Oncogenic_Other": [  # 剩下的癌症信号
        "CANCER", "TUMOR", "MYC", "KRAS", "IMMUNE", "VIRAL"
    ]
}

# 优先级列表：如果一个基因同时命中多个，按此顺序归类
PRIORITY = [
    "DSB_Repair",
    "SSB_Excision",
    "Replication",
    "Chromatin_Epi",
    "Checkpoint_p53",
    "Oncogenic_Other"
]


def run_multi_classification(input_csv, output_csv):
    try:
        df = pd.read_csv(input_csv, encoding="utf-8-sig")
    except:
        df = pd.read_csv(input_csv, encoding="gbk")

    print(f"读取数据: {len(df)} 行")

    def parse_pathways(val):
        try:
            return [x.strip().upper() for x in ast.literal_eval(val)]
        except:
            s = str(val).strip().upper()
            return [s] if s and s != 'NAN' else []

    # --- 核心功能 1: 判定主分类 ---
    def get_main_category(pathways):
        # 收集该基因命中了哪些分类
        hits = set()
        combined_text = " ".join(pathways)

        for cat, kw_list in KEYWORDS.items():
            if any(kw in combined_text for kw in kw_list):
                hits.add(cat)

        # 按优先级返回
        for cat in PRIORITY:
            if cat in hits:
                return cat
        return "Other"

    # --- 核心功能 2: 计算UMAP所需的数值特征 ---
    def calculate_scores(pathways):
        combined_text = " ".join(pathways)
        scores = {}
        for cat, kw_list in KEYWORDS.items():
            # 计算关键词出现的总次数，作为该维度的得分
            count = sum(combined_text.count(kw) for kw in kw_list)
            scores[f"Score_{cat}"] = count
        return pd.Series(scores)

    print("正在执行多分类及特征计算...")

    # 1. 转换通路列
    df["pathway_list"] = df["STANDARD_NAME"].apply(parse_pathways)

    # 2. 获取主分类
    df["Category"] = df["pathway_list"].apply(get_main_category)

    # 3. 计算分数矩阵 (用于画图)
    scores_df = df["pathway_list"].apply(calculate_scores)
    result_df = pd.concat([df, scores_df], axis=1)

    # 清理临时列
    if "pathway_list" in result_df.columns:
        del result_df["pathway_list"]

    print("\n分类统计结果:")
    print(result_df["Category"].value_counts())

    result_df.to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"\n结果已保存至: {output_csv}")


if __name__ == "__main__":
    input_file = "/Users/gzy2520/Desktop/lac/symbol_and_pathway.csv"
    output_file = "/Users/gzy2520/Desktop/lac/divided_symbol.csv"

    run_multi_classification(input_file, output_file)