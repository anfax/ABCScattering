#!/usr/bin/env python3
import os
import re
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import concurrent.futures

# 参数设置
input_pattern = "Aoutput*.txt"
header = "#  Coll.                   a    v    j    k    a'   v'   j'   k'     Re(S)     Im(S)     |S|^2"

# 构建或清空目录 s_matrix
def build_dir(d):
    if os.path.isdir(d):
        for f in os.listdir(d):
            fpath = os.path.join(d, f)
            if os.path.isfile(fpath):
                os.remove(fpath)
    else:
        os.makedirs(d, exist_ok=True)

# 从文件中提取能量列表，匹配 "E(eV) =" 行，提取 "=" 后第一个数字
def parse_energy_list(lines):
    energy_list = []
    for line in lines:
        if "E(eV) =" in line:
            m = re.search(r"E\(eV\)\s*=\s*([^\s]+)", line)
            if m:
                energy_list.append(m.group(1).strip())
    return energy_list

# 从文件中提取通道数据：读取首个"Channel"行后连续非空、不为分隔行的行，
# 每行提取 a, v, j, k（列表按顺序存储），遇到 a=="2" 时停止处理
def parse_channels(lines):
    channels = []
    flag = False
    for line in lines:
        if not flag:
            if line[:7] == "Channel":
                flag = True
            continue
        if line.strip() == "":
            break
        if len(line) > 1 and line[1] == "-":
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            _ = int(parts[0])
        except:
            continue
        channels.append({
            'a': parts[1],
            'v': parts[2],
            'j': parts[3],
            'k': parts[4]
        })
    return channels

# 对某一对通道（n, m）进行数据匹配处理：
# 在整个文件中通过正则匹配连续出现的 ch_a 与 ch_b，
# 按行将能量列表与匹配到的数据行进行“粘贴”，保存到数据文件中。
def process_channel_pair(n, m, channels, all_lines, energy_list, s_dir):
    ch_a = "    ".join([channels[n]['a'], channels[n]['v'], channels[n]['j'], channels[n]['k']])
    ch_b = "    ".join([channels[m]['a'], channels[m]['v'], channels[m]['j'], channels[m]['k']])
    pattern = re.compile(re.escape(ch_a) + r"\s+" + re.escape(ch_b))
    buf_lines = [line.rstrip() for line in all_lines if pattern.search(line)]
    if len(buf_lines) == 0:
        return None
    n_lines = min(len(energy_list), len(buf_lines))
    paste_lines = []
    for i in range(n_lines):
        paste_lines.append(energy_list[i] + "\t" + buf_lines[i])
    filename = os.path.join(s_dir,
        f"a={channels[n]['a']}-{channels[m]['a']}_v={channels[n]['v']}-{channels[m]['v']}_j={channels[n]['j']}-{channels[m]['j']}_k={channels[n]['k']}-{channels[m]['k']}.dat")
    with open(filename, "w") as f:
        f.write(header + "\n" + "\n".join(paste_lines))
    print(f"生成文件: {filename}")
    return filename

def main():
    files = glob.glob(input_pattern)
    if not files:
        print("未找到匹配 Aoutput*.txt 的文件。")
        sys.exit(1)
    infile = files[0]
    with open(infile, "r") as f:
        lines = f.readlines()

    energy_list = parse_energy_list(lines)
    if not energy_list:
        print("无法从输入文件中提取能量列表")
        sys.exit(1)

    channels = parse_channels(lines)
    if not channels:
        print("无法从输入文件中提取通道数据")
        sys.exit(1)
    valid_channels = []
    for ch in channels:
        if ch['a'] == "2":
            break
        valid_channels.append(ch)
    channels = valid_channels
    n_max = len(channels)
    s_dir = os.path.join(os.getcwd(), "s_matrix")
    build_dir(s_dir)

    generated_files = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for n in range(n_max):
            for m in range(n, n_max):
                futures.append(executor.submit(process_channel_pair, n, m, channels, lines, energy_list, s_dir))
        for future in concurrent.futures.as_completed(futures):
            res = future.result()
            if res:
                generated_files.append(res)
    if not generated_files:
        print("未生成任何数据文件。")
        sys.exit(1)

    # 绘图：对所有数据文件读取第一列（碰撞能量）与最后一列（|S|^2）的数据
    plt.figure(figsize=(10, 6))
    for f in generated_files:
        try:
            data = np.loadtxt(f, comments='#', delimiter="\t")
            if data.ndim == 1 or data.shape[0] < 2:
                continue
            coll = data[:, 0].astype(float)
            S_squared = data[:, -1].astype(float)
            plt.plot(coll, S_squared, label=os.path.basename(f))
        except Exception as e:
            print(f"绘制 {f} 时出错: {e}")
    plt.xlabel("Collision Energy")
    plt.ylabel("|S|^2")
    plt.title("S-matrix Plot")
    plt.legend(fontsize='small', loc='best')
    plt.tight_layout()
    plotfile = "s_matrix_plot.png"
    plt.savefig(plotfile)
    print(f"图像已保存为: {plotfile}")
    plt.show()

if __name__ == '__main__':
    main()