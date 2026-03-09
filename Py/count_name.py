import pandas as pd
import ast
from collections import Counter


def extract_and_count_pathways(input_csv, output_csv=None):
    try:
        # 1. 读取CSV文件
        try:
            df = pd.read_csv(input_csv, encoding='utf-8-sig')
        except UnicodeDecodeError:
            df = pd.read_csv(input_csv, encoding='gbk')

        print(f"成功读取文件，共 {len(df)} 行数据。正在统计...")

        # 2. 使用 Counter 来存储通路及其出现次数
        pathway_counter = Counter()

        # 3. 遍历 STANDARD_NAME 列
        for item in df['STANDARD_NAME']:
            # 检查是否为空值
            if pd.isna(item):
                continue

            try:
                # 将字符串 "['A', 'B']" 转换为列表 ['A', 'B']
                pathway_list = ast.literal_eval(item)

                if isinstance(pathway_list, list):
                    # 如果是列表，批量更新计数器
                    pathway_counter.update(pathway_list)
                else:
                    # 如果解析出来的不是列表，作为单个条目计数
                    pathway_counter[str(pathway_list)] += 1

            except (ValueError, SyntaxError):
                # 如果格式不是标准的Python列表字符串，直接作为原始字符串计数
                pathway_counter[str(item)] += 1

        # 4. 将计数器转换为 DataFrame以便展示和保存
        # items() 返回 [(name, count), (name, count)...]
        result_df = pd.DataFrame(pathway_counter.items(), columns=['Pathway_Name', 'Count'])

        # 5. 排序：先按 Count 降序(从大到小)，如果相同则按 Pathway_Name 升序
        result_df = result_df.sort_values(by=['Count', 'Pathway_Name'], ascending=[False, True])

        # 6. 如果指定了输出路径，保存结果
        if output_csv:
            result_df.to_csv(output_csv, index=False, encoding='utf-8-sig')
            print(f"统计结果已保存至: {output_csv}")

        return result_df

    except Exception as e:
        print(f"发生错误: {e}")
        return pd.DataFrame()


if __name__ == "__main__":
    input_file = "/Users/gzy2520/Desktop/lac/symbol_and_pathway.csv"
    output_file = "/Users/gzy2520/Desktop/lac/pathway_counts.csv"

    # 运行函数
    df_result = extract_and_count_pathways(input_file, output_file)

    # 打印前 50 个
    if not df_result.empty:
        print("\n====== 通路频率统计 (Top 50) ======")
        # 设置pandas打印选项，防止显示不全
        pd.set_option('display.max_rows', 100)
        pd.set_option('display.max_colwidth', None)
        print(df_result.head(50).to_string(index=False))
        print(f"\n总计检测到 {len(df_result)} 种不同的通路名称。")