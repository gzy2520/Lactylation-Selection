from pathlib import Path
import pandas as pd

def parse_full_msigdb_tsv(tsv_file_path):
    pathway_map = {}  # key: 基因Symbol, value: 所有对应通路
    current_pathway = None
    current_genes = []

    # 逐行读取TSV文件
    with open(tsv_file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue  # 跳过空行

            # 按制表符分割
            parts = [p.strip() for p in line.split('\t') if p.strip()]
            if len(parts) < 2:
                continue  # 跳过无效行（列数不足）

            # 1. 识别新通路的开始（STANDARD_NAME行）
            if parts[0] == 'STANDARD_NAME':
                # 先保存上一个通路的基因映射（避免遗漏）
                if current_pathway and current_genes:
                    # 基因去重+统一大写（与CSV基因格式对齐）
                    unique_genes = list(set([g.upper().strip() for g in current_genes if g.strip()]))
                    for gene in unique_genes:
                        if gene not in pathway_map:
                            pathway_map[gene] = []
                        # 通路去重（避免同一基因重复关联同一通路）
                        if current_pathway not in pathway_map[gene]:
                            pathway_map[gene].append(current_pathway)
                # 更新当前通路名，重置基因列表
                current_pathway = parts[1]
                current_genes = []

            # 2. 识别当前通路的基因列表
            elif parts[0] == 'GENE_SYMBOLS' and len(parts) > 1:
                # 拆分基因
                genes = [g.strip() for g in parts[1].split(',') if g.strip()]
                current_genes.extend(genes)

        # 处理最后一个通路
        if current_pathway and current_genes:
            unique_genes = list(set([g.upper().strip() for g in current_genes if g.strip()]))
            for gene in unique_genes:
                if gene not in pathway_map:
                    pathway_map[gene] = []
                if current_pathway not in pathway_map[gene]:
                    pathway_map[gene].append(current_pathway)

    print(f"共识别 {len(pathway_map)} 个基因的通路关联")
    return pathway_map


def clean_csv_genes(csv_file_path):
    df = pd.read_csv(csv_file_path)
    # 清理
    df_clean = df.dropna(subset=["Symbol"])
    df_clean["Symbol"] = df_clean["Symbol"].str.upper().str.strip()  # 统一格式
    return df_clean


def match_genes_to_full_pathways(csv_clean, pathway_map):

    # 为每个基因匹配所有通路
    def get_all_matched_pathways(sym):
        sym_upper = sym.upper().strip()
        # 无匹配时返回标记，有匹配时返回所有通路
        matched = pathway_map.get(sym_upper, ["无匹配通路"])
        return sorted(matched) if matched != ["无匹配通路"] else matched

    # 应用匹配逻辑
    csv_clean["STANDARD_NAME"] = csv_clean["Symbol"].apply(get_all_matched_pathways)

    # 统计匹配情况
    match_stats = {
        "总有效基因数": len(csv_clean),
        "匹配成功数（≥1个通路）": len(csv_clean[csv_clean["STANDARD_NAME"].apply(lambda x: x[0] != "无匹配通路")]),
        "无匹配数": len(csv_clean[csv_clean["STANDARD_NAME"].apply(lambda x: x[0] == "无匹配通路")]),
        "无匹配基因列表": csv_clean[csv_clean["STANDARD_NAME"].apply(lambda x: x[0] == "无匹配通路")]["Symbol"].tolist()
    }
    return csv_clean, match_stats


if __name__ == "__main__":
    TSV_FILE_PATH = Path("/Users/gzy2520/Desktop/lac/genesets.tsv")
    CSV_FILE_PATH = Path("/Users/gzy2520/Desktop/lac/lac_ddr_unique.csv")
    OUTPUT_FILE_PATH = Path("/Users/gzy2520/Desktop/lac/symbol_and_pathway.csv")

    # 1. 读取并解析完整TSV文件
    print(f"正在读取TSV文件：{TSV_FILE_PATH}")
    full_pathway_map = parse_full_msigdb_tsv(TSV_FILE_PATH)

    # 2. 清理CSV基因（去空值、统一格式）
    print(f"\n正在读取CSV基因文件：{CSV_FILE_PATH}")
    clean_csv = clean_csv_genes(CSV_FILE_PATH)

    # 3. 匹配所有通路
    print("\n正在匹配基因与通路...")
    final_result, stats = match_genes_to_full_pathways(clean_csv, full_pathway_map)

    # 4. 输出统计结果
    print("\n" + "=" * 50)
    print("=== 完整通路匹配结果统计 ===")
    for key, value in stats.items():
        print(f"{key}：{value}")
    print("=" * 50)

    # 5. 保存最终结果（包含所有通路）
    final_result.to_csv(OUTPUT_FILE_PATH, index=False, encoding="utf-8-sig")
    print(f"\n最终结果已保存至：{OUTPUT_FILE_PATH}")
    print("\n前3行结果预览（Symbol + 对应所有通路）：")
    print(final_result[["Symbol", "STANDARD_NAME"]].head(3))

