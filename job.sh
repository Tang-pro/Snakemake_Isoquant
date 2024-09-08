#在LSF提交snakemake任务
snakemake --snakefile isoquant_smk.py --cluster "bsub -J merge -n 4 -q gpu" -j 100
#将工作流程可视化
module load graphviz/2.40.1
snakemake --snakefile isoquant_smk.py --dag | dot -Tpng > dag.png
