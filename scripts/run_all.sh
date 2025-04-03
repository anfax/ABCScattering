#!/bin/bash
set -u
set -e

# 自动检测调度系统资源 (gres)
if command -v scontrol >/dev/null 2>&1; then
    scheduler="slurm"
    echo "Detected SLURM scheduler"
    # 自动获取 gres 值
    GRES=$(scontrol show node | grep -m1 'Gres=' | sed -n 's/.*Gres:\([^,[:space:]]*\).*/\1/p')
    if [ -z "$GRES" ]; then
         GRES=""
    fi
    echo "Auto-detected gres: $GRES"
else
    scheduler="pbs"
    echo "Detected PBS scheduler"
fi

# 如果使用 SLURM，则更新 slurm_setup.sh 中的 gres 变量
if [ "$scheduler" == "slurm" ]; then
    if [ "$GRES" == "" ]; then
        # 如果自动检测到 gres 为 all，则删除 gres 设置行
        sed -i '/^#SBATCH --gres=/d' "$(dirname "$0")/slurm_setup.sh"
    else
        sed -i 's/^gres=".*"/gres="'"$GRES"'"/' "$(dirname "$0")/slurm_setup.sh"
    fi
fi

echo "---------------------------------------------------------"
echo "Step 1: 运行 abc_workflow.sh 构建目录和输入文件"
bash "$(dirname "$0")/abc_workflow.sh"
echo "---------------------------------------------------------"

if [ "$scheduler" == "slurm" ]; then
    echo "Step 2: 生成 SLURM 作业脚本"
    bash "$(dirname "$0")/slurm_setup.sh"
    echo "Step 3: 提交 SLURM 作业"
    bash "$(dirname "$0")/sbatch_all.sh"
else
    echo "Step 2: 生成 PBS 作业脚本"
    bash "$(dirname "$0")/pbs_setup.sh"
    echo "Step 3: 提交 PBS 作业"
    bash "$(dirname "$0")/qsub_all.sh"
fi
echo "---------------------------------------------------------"

echo "作业已经提交，请等待作业完成后再进行后处理。"
read -p "完成后请按 Enter 键继续运行后处理步骤..." dummy

echo "Step 4: 处理 S 矩阵数据 (get_s.sh)"
bash "$(dirname "$0")/get_s.sh"
echo "---------------------------------------------------------"

echo "Step 5: 提取通道信息 (get_channels.sh)"
bash "$(dirname "$0")/get_channels.sh"
echo "---------------------------------------------------------"

echo "所有步骤已完成。"
