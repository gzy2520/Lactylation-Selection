from pathlib import Path
import pandas as pd


def extract_all_genes_from_msigdb_tsv(tsv_file_path):
    """
    解析 TSV 文件，提取所有基因及其对应的通路
    返回字典: {基因Symbol: [通路1, 通路2, ...]}
    """
    pathway_map = {}  # key: 基因Symbol, value: 通路列表
    current_pathway = None

    print(f"正在解析文件: {tsv_file_path}")

    with open(tsv_file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            # 按制表符分割
            parts = [p.strip() for p in line.split('\t')]
            if len(parts) < 2:
                continue

            # 1. 获取通路名称
            if parts[0] == 'STANDARD_NAME':
                current_pathway = parts[1]

            # 2. 获取该通路下的基因并建立映射
            elif parts[0] == 'GENE_SYMBOLS' and current_pathway:
                # 基因之间通常用逗号分隔
                genes = [g.strip().upper() for g in parts[1].split(',') if g.strip()]

                for gene in genes:
                    if gene not in pathway_map:
                        pathway_map[gene] = set()  # 使用集合自动去重
                    pathway_map[gene].add(current_pathway)

    print(f"解析完成，共发现 {len(pathway_map)} 个唯一基因")
    return pathway_map


def save_to_csv(pathway_map, output_path):
    """
    将提取的结果转换为 DataFrame 并保存
    """
    data = []
    # 排序基因使结果更整齐
    for gene in sorted(pathway_map.keys()):
        # 将通路集合转为逗号分隔的字符串
        pathways_str = "|".join(sorted(pathway_map[gene]))
        data.append({
            "Gene_Symbol": gene,
            "Pathways": pathways_str,
            "Pathway_Count": len(pathway_map[gene])
        })

    df = pd.DataFrame(data)
    df.to_csv(output_path, index=False, encoding="utf-8-sig")
    print(f"结果已成功保存至: {output_path}")
    return df


if __name__ == "__main__":
    # --- 配置路径 ---
    TSV_FILE_PATH = Path("/Users/gzy2520/Desktop/lac/genesets.tsv")
    OUTPUT_FILE_PATH = Path("/Users/gzy2520/Desktop/lac/all_genes_pathways_summary.csv")

    # 1. 提取所有数据
    all_gene_data = extract_all_genes_from_msigdb_tsv(TSV_FILE_PATH)

    # 2. 保存结果
    if all_gene_data:
        df_result = save_to_csv(all_gene_data, OUTPUT_FILE_PATH)

        # 3. 预览结果
        print("\n结果预览 (前10行):")
        print(df_result.head(10))

        # 统计一些有趣的数据
        max_pathways = df_result["Pathway_Count"].max()
        most_active_gene = df_result[df_result["Pathway_Count"] == max_pathways]["Gene_Symbol"].iloc[0]
        print(f"\n统计信息:")
        print(f"- 参与通路最多的基因是: {most_active_gene} (参与了 {max_pathways} 个通路)")
    else:
        print("未能在文件中找到有效数据，请检查 TSV 格式。")