# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 绘制双性状 MiXeR 图（Venn / 散点）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
#docker run --rm \
#  -v "/path/to/project/GWAS/MiXeR:/home" \
#  -w /home \
#  ghcr.io/precimed/gsa-mixer:latest \
#  python /tools/mixer/precimed/mixer_figures.py two \
#  --json-fit /home/results/figures/single_right_knee_vs_single_right_hip.fit.json \
#  --json-test /home/results/figures/single_right_knee_vs_single_right_hip.test.json \
#  --out /home/results/figures/single_right_knee_vs_single_right_hip \
#  --trait1 single_right_knee \
#  --trait2 single_right_hip \
#  --statistic mean std \
#  --ext png


docker run --rm \
  -v "/path/to/project/GWAS/MiXeR:/home" \
  -w /home \
  ghcr.io/precimed/gsa-mixer:latest \
  python /tools/mixer/precimed/mixer_figures.py two \
  --json-fit /home/results/figures/double_right_vs_single_right_knee.fit.json \
  --json-test /home/results/figures/double_right_vs_single_right_knee.test.json \
  --out /home/results/figures/double_right_vs_single_right_knee \
  --trait1 double_right \
  --trait2 single_right_knee \
  --statistic mean std \
  --ext png

docker run --rm \
  -v "/path/to/project/GWAS/MiXeR:/home" \
  -w /home \
  ghcr.io/precimed/gsa-mixer:latest \
  python /tools/mixer/precimed/mixer_figures.py two \
  --json-fit /home/results/figures/double_right_vs_single_right_hip.fit.json \
  --json-test /home/results/figures/double_right_vs_single_right_hip.test.json \
  --out /home/results/figures/double_right_vs_single_right_hip \
  --trait1 double_right \
  --trait2 single_right_hip \
  --statistic mean std \
  --ext png  