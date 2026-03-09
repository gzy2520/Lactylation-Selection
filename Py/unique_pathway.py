import pandas as pd
import ast


def extract_unique_pathways(csv_file_path):
    try:
        # 读取CSV文件
        df = pd.read_csv(csv_file_path, encoding='utf-8-sig')

        # 使用 set 集合来存储，自动去除重复项
        unique_names = set()

        # 遍历 STANDARD_NAME 列
        for item in df['STANDARD_NAME']:
            # 检查是否为空值
            if pd.isna(item):
                continue

            try:
                # ast.literal_eval 将字符串 "['A', 'B']" 转换为列表 ['A', 'B']
                pathway_list = ast.literal_eval(item)

                # 将列表中的每个通路名称添加到集合中
                if isinstance(pathway_list, list):
                    for name in pathway_list:
                        unique_names.add(name)
                else:
                    # 如果解析出来的不是列表，直接添加（以防万一格式不同）
                    unique_names.add(str(pathway_list))

            except (ValueError, SyntaxError):
                # 如果格式不是标准的Python列表字符串，直接添加
                unique_names.add(item)

        # 将集合转换为排序后的列表
        sorted_pathways = sorted(list(unique_names))

        return sorted_pathways

    except Exception as e:
        print(f"发生错误: {e}")
        return []



pathways = extract_unique_pathways("/Users/gzy2520/Desktop/lac/symbol_and_pathway.csv")
for p in pathways:
    print(p)