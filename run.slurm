#! /bin/bash
#SBATCH --job-name filter
#SBATCH -o filter-%j.out # 把输出结果STDOUT保存在哪一个文件
#SBATCH -e filter-%j.err # 把报错结果STDERR保存在哪一个文件
#SBATCH -n 960 #指定每个进程使用核数，不指定默认为1
#SBATCH -t 12-00:00 # 运行总时间，天数-小时数-分钟， D-HH:MM
#SBATCH -p batch # 提交到哪一个分区i
#SBATCH -w node10

source activate design
cd /lustre/home/yangsen/workspace/ProteinFilter
bash run.sh