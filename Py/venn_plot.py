import pandas as pd
import matplotlib.pyplot as plt
from venn import venn
import math
from matplotlib.lines import Line2D

# ================= 配置区域 =================
# 请确保此处路径正确
INPUT_FILE = '/Users/gzy2520/Desktop/lac/filter/venn_input.csv'
OUTPUT_FILE = 'venn_final_corrected.png'

# 颜色定义
COLORS_MAP = {
    "HR": "#E64B35",
    "NHEJ": "#4DBBD5",
    "BER": "#00A087",
    "NER": "#3C5488",
    "CP": "#F39B7F"
}

# 设置一个足够大的数，确保显示所有基因不折叠
MAX_GENES_SHOW = 1000


# ================= 核心逻辑 =================

def get_labels_and_data(df):
    """读取数据并确保列顺序一致"""
    # 强制排序，保证每次运行顺序一致
    labels = sorted([c for c in df.columns if c in COLORS_MAP])
    return labels


def calibrate_coordinates(labels):
    """
    【关键步骤】几何校准
    画一个只有数字 ID 的虚拟图，获取每个逻辑区域的准确 (x, y) 坐标。
    """
    print("正在进行几何坐标校准...")

    # 构造虚拟数据：让区域 K (二进制值) 拥有 K 个元素
    # 这样图上显示的数字就是它的 ID
    dummy_sets = {label: set() for label in labels}

    # 遍历所有可能的组合 (1 到 31)
    for i in range(1, 32):
        # 生成 i 个唯一标识符
        dummy_items = {f"dummy_{i}_{x}" for x in range(i)}

        # 将这些标识符放入对应的集合中
        for idx, label in enumerate(labels):
            # 检查第 idx 位是否为 1
            if (i >> idx) & 1:
                dummy_sets[label].update(dummy_items)

    # 画在临时画布上
    fig_cal, ax_cal = plt.subplots()
    venn(dummy_sets, ax=ax_cal, legend_loc=None)

    # 提取坐标映射：ID -> (x, y)
    id_to_coord = {}
    for t in ax_cal.texts:
        try:
            val = int(t.get_text())  # 获取图上的数字（即我们的 ID）
            id_to_coord[val] = t.get_position()
        except ValueError:
            continue

    plt.close(fig_cal)  # 关闭临时图
    return id_to_coord


def calculate_real_intersections(df, labels):
    """计算真实数据的交集内容"""
    real_data = {}

    # 这里的 i 必须与校准时的逻辑完全一致
    for i in range(1, 32):
        condition = pd.Series([True] * len(df))

        for idx, label in enumerate(labels):
            if (i >> idx) & 1:
                condition = condition & (df[label] == 1)
            else:
                condition = condition & (df[label] == 0)

        symbols = df[condition]['Symbol'].tolist()
        if symbols:  # 只记录非空集合
            real_data[i] = symbols

    return real_data


def find_nearest_id(x, y, id_to_coord):
    """找到距离当前坐标最近的校准点 ID"""
    min_dist = float('inf')
    best_id = -1

    for pid, (px, py) in id_to_coord.items():
        dist = math.sqrt((x - px) ** 2 + (y - py) ** 2)
        if dist < min_dist:
            min_dist = dist
            best_id = pid

    return best_id


def main():
    # 1. 读取数据
    try:
        df = pd.read_csv(INPUT_FILE)
    except FileNotFoundError:
        print(f"找不到文件: {INPUT_FILE}")
        return

    labels = get_labels_and_data(df)
    colors_list = [COLORS_MAP[l] for l in labels]

    # 2. 获取坐标映射 (Calibration)
    coord_map = calibrate_coordinates(labels)

    # 3. 计算真实交集内容
    real_intersections = calculate_real_intersections(df, labels)

    # 构造真实绘图用的集合字典
    real_sets = {}
    for label in labels:
        real_sets[label] = set(df[df[label] == 1]['Symbol'])

    # 4. 正式绘图
    fig, ax = plt.subplots(figsize=(14, 14))  # 画布稍微加大一点

    venn(real_sets, ax=ax, cmap=colors_list, alpha=0.5, legend_loc=None)

    plt.title("DNA Repair Pathway Classification", fontsize=20, weight='bold', pad=40)

    # 5. 智能标注
    annotations_count = 0

    # 计算所有点的中心，用于发散引线
    if coord_map:
        center_x = sum(p[0] for p in coord_map.values()) / len(coord_map)
        center_y = sum(p[1] for p in coord_map.values()) / len(coord_map)
    else:
        center_x, center_y = 0.5, 0.45

    for t_obj in ax.texts:
        x, y = t_obj.get_position()

        # 找到这个位置对应的逻辑 ID
        region_id = find_nearest_id(x, y, coord_map)

        # 获取该 ID 对应的真实基因
        genes = real_intersections.get(region_id, [])
        count = len(genes)

        # --- 核心修改：处理交集内的数字显示 ---
        if count == 0:
            t_obj.set_visible(False)
            continue
        else:
            t_obj.set_visible(True)  # 保持数字可见
            t_obj.set_fontweight('bold')  # 加粗数字
            t_obj.set_fontsize(12)  # 字体稍大
            # t_obj.set_text(str(count))  # 确保显示的是数量

        annotations_count += 1

        # 显示全部 Symbol
        label_text = "\n".join(genes)

        # 计算引线位置 (向外发散算法)
        dx = x - center_x
        dy = y - center_y
        dist = math.sqrt(dx ** 2 + dy ** 2)

        # 拉伸系数：离中心越近拉得越远，避免挤在中间
        stretch = 1.7 + (0.5 / (dist + 0.1))

        # 特殊处理：如果点就在正中心附近
        if dist < 0.05:
            new_x, new_y = x + 0.4, y + 0.4
        else:
            new_x = center_x + dx * stretch
            new_y = center_y + dy * stretch

        # 微调防重叠 (上下交错)
        if annotations_count % 2 == 0:
            new_y += 0.06
        else:
            new_y -= 0.06

        # 绘制引线标注
        ax.annotate(
            text=label_text,
            xy=(x, y),
            xytext=(new_x, new_y),
            fontsize=9,
            fontweight='bold',
            color='black',
            ha='center', va='center',
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.95),
            arrowprops=dict(
                arrowstyle='->',
                color='#333333',
                linewidth=1.2,
                connectionstyle="arc3,rad=0.1"  # 弧形引线
            )
        )

    # 6. 添加图例 (放在底部横排)
    legend_elements = [Line2D([0], [0], marker='o', color='w', label=label,
                              markerfacecolor=COLORS_MAP[label], markersize=15)
                       for label in labels]

    ax.legend(handles=legend_elements,
              loc="upper center",
              bbox_to_anchor=(0.5, 0.02),  # 放在底部
              ncol=len(labels),  # 横向排列
              title="Pathway",
              fontsize=12,
              frameon=False)

    # 保存图片
    plt.savefig(OUTPUT_FILE, dpi=300, bbox_inches='tight', pad_inches=0.2)
    print(f"绘图完成，已保存为: {OUTPUT_FILE}")
    plt.show()


if __name__ == "__main__":
    main()