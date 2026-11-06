# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 合并 20 次重复结果为单一参数集
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
# Combine double_right fit1
docker run --rm \
  -v "/path/to/project/GWAS/MiXeR:/home" \
  -w /home \
  ghcr.io/precimed/gsa-mixer:latest \
  python /tools/mixer/precimed/mixer_figures.py combine \
  --json /home/results/fit/double_right.fit.rep@.json \
  --out /home/results/figures/double_right.fit

# Combine double_right–knee fit2
docker run --rm \
  -v "/path/to/project/GWAS/MiXeR:/home" \
  -w /home \
  ghcr.io/precimed/gsa-mixer:latest \
  python /tools/mixer/precimed/mixer_figures.py combine \
  --json /home/results/fit/double_right_vs_single_right_knee.fit.rep@.json \
  --out /home/results/figures/double_right_vs_single_right_knee.fit

# Combine double_right–hip fit2
docker run --rm \
  -v "/path/to/project/GWAS/MiXeR:/home" \
  -w /home \
  ghcr.io/precimed/gsa-mixer:latest \
  python /tools/mixer/precimed/mixer_figures.py combine \
  --json /home/results/fit/double_right_vs_single_right_hip.fit.rep@.json \
  --out /home/results/figures/double_right_vs_single_right_hip.fit

# Combine double_right test1
docker run --rm \
  -v "/path/to/project/GWAS/MiXeR:/home" \
  -w /home \
  ghcr.io/precimed/gsa-mixer:latest \
  python /tools/mixer/precimed/mixer_figures.py combine \
  --json /home/results/test/double_right.test.rep@.json \
  --out /home/results/figures/double_right.test

# Combine double_right–knee test2
docker run --rm \
  -v "/path/to/project/GWAS/MiXeR:/home" \
  -w /home \
  ghcr.io/precimed/gsa-mixer:latest \
  python /tools/mixer/precimed/mixer_figures.py combine \
  --json /home/results/test/double_right_vs_single_right_knee.test.rep@.json \
  --out /home/results/figures/double_right_vs_single_right_knee.test

# Combine double_right–hip test2
docker run --rm \
  -v "/path/to/project/GWAS/MiXeR:/home" \
  -w /home \
  ghcr.io/precimed/gsa-mixer:latest \
  python /tools/mixer/precimed/mixer_figures.py combine \
  --json /home/results/test/double_right_vs_single_right_hip.test.rep@.json \
  --out /home/results/figures/double_right_vs_single_right_hip.test

