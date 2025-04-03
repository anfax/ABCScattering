import os
import re
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

def parse_scat_file(filepath):
    """解析单个散射文件，返回结构化数据"""
    data = {
        'energy': None,
        'transitions': defaultdict(list),
        'v_total': defaultdict(float)
    }
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
        
    # 提取能量值（如果文件包含）
    energy_lines = [l for l in lines if "E(eV)" in l]
    if energy_lines:
        energy_line = energy_lines[0]
        match = re.search(r"[-+]?\d*\.?\d+E?[-+]?\d+", energy_line)
        if match:
            data['energy'] = float(match.group())
    
    # 解析数据表
    for line in lines:
        cols = re.split(r'\s+', line.strip())
        # 跳过空行和注释行
        if len(cols) < 7 or cols[0] == '#':
            continue
            
        # 处理跃迁数据行 (格式: a v j k a' v' j' |S|^2)
        if cols[0].isdigit() and '.' in cols[-1]:
            try:
                transition = {
                    'a': int(cols[0]),
                    'v': int(cols[1]),
                    'j': int(cols[2]),
                    'k': int(cols[3]),
                    'a_prime': int(cols[4]),
                    'v_prime': int(cols[5]),
                    'j_prime': int(cols[6]),
                    '|S|^2': float(cols[-1])
                }
                data['transitions'][transition['v_prime']].append(transition)
                # 自动计算v'总概率
                data['v_total'][transition['v_prime']] += transition['|S|^2']
            except Exception as e:
                print(f"解析行错误: {line.strip()} -> {str(e)}")
                
    return data

def visualize(data, output_dir):
    """生成三种可视化图表"""
    os.makedirs(output_dir, exist_ok=True)
    base_name = f"E_{data['energy']:.3f}eV" if data['energy'] else "Unknown_Energy"

    # 1. 振动能级总分布
    if data['v_total']:
        v_primes = sorted(data['v_total'].keys())
        totals = [data['v_total'][v] for v in v_primes]
        plt.figure(figsize=(8, 4))
        plt.bar(v_primes, totals, color='teal')
        plt.title(f"Vibrational Distribution @ {base_name}")
        plt.xlabel("Final v'")
        plt.ylabel("Total |S|²")
        plt.yscale('log')
        plt.savefig(os.path.join(output_dir, f"{base_name}_vibrational.png"))
        plt.close()

    # 2. 各v'的转动分布
    if data['transitions']:
        for v_prime in data['transitions']:
            transitions = data['transitions'][v_prime]
            j_primes = sorted({t['j_prime'] for t in transitions})
            s2_values = [sum(t['|S|^2'] for t in transitions if t['j_prime'] == j) 
                        for j in j_primes]
            
            plt.figure(figsize=(6, 4))
            plt.bar(j_primes, s2_values, color='orange')
            plt.title(f"Rotational Distribution (v'={v_prime}) @ {base_name}")
            plt.xlabel("Final j'")
            plt.ylabel("|S|²")
            plt.savefig(os.path.join(output_dir, f"{base_name}_v{v_prime}_rotational.png"))
            plt.close()

    # 3. 热图（仅当存在跃迁数据时生成）
    if data['transitions']:
        try:
            max_v = max(data['transitions'].keys())
            all_j = [t['j_prime'] for v in data['transitions'] for t in data['transitions'][v]]
            max_j = max(all_j) if all_j else 0
            
            heatmap = np.zeros((max_v + 1, max_j + 1))
            for v in data['transitions']:
                for t in data['transitions'][v]:
                    heatmap[v, t['j_prime']] += t['|S|^2']
            
            plt.figure(figsize=(8, 6))
            plt.imshow(heatmap, cmap='viridis', origin='lower', aspect='auto')
            plt.colorbar(label='|S|²')
            plt.title(f"Heatmap @ {base_name}")
            plt.xlabel("j'")
            plt.ylabel("v'")
            plt.savefig(os.path.join(output_dir, f"{base_name}_heatmap.png"))
            plt.close()
        except Exception as e:
            print(f"生成热图失败: {str(e)}")

if __name__ == "__main__":
    input_dir = "./"
    output_dir = "fig"
    
    # 批量处理所有.dat文件
    for filename in os.listdir(input_dir):
        if filename.endswith(".dat"):
            filepath = os.path.join(input_dir, filename)
            print(f"Processing: {filepath}")
            try:
                data = parse_scat_file(filepath)
                visualize(data, output_dir)
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}")