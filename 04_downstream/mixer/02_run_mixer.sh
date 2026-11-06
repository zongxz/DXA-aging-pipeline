#!/bin/bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 20 次重复的 fit1 / fit2 / test1 / test2（docker 运行 gsa-mixer）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
set -euo pipefail

MIXER_DIR="/path/to/project/GWAS/MiXeR"

IMAGE="ghcr.io/precimed/gsa-mixer:latest"

mkdir -p "$MIXER_DIR/results/fit"
mkdir -p "$MIXER_DIR/results/test"
mkdir -p "$MIXER_DIR/results/logs"

for REP in $(seq 1 20); do

    echo "======================================================"
    echo "Starting MiXeR replicate $REP"
    echo "======================================================"

    # ==================================================
    # 1. double_right fit1
    # ==================================================

    docker run --rm \
      -v "$MIXER_DIR:/home" \
      -w /home \
      "$IMAGE" \
      python /tools/mixer/precimed/mixer.py fit1 \
      --lib /tools/mixer/src/build/lib/libbgmg.so \
      --bim-file /home/ref/1000G.EUR.QC.@.bim \
      --ld-file /home/ref/1000G.EUR.QC.@.run4.ld \
      --trait1-file /home/sumstats/double_right.sumstats.gz \
      --extract /home/ref/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps \
      --chr2use 1-22 \
      --threads 24 \
      --out /home/results/fit/double_right.fit.rep${REP} \
      2>&1 | tee "$MIXER_DIR/results/logs/double_right.fit.rep${REP}.log"


    # ==================================================
    # 2. double_right vs single_right_hip fit2
    # ==================================================

    docker run --rm \
      -v "$MIXER_DIR:/home" \
      -w /home \
      "$IMAGE" \
      python /tools/mixer/precimed/mixer.py fit2 \
      --lib /tools/mixer/src/build/lib/libbgmg.so \
      --bim-file /home/ref/1000G.EUR.QC.@.bim \
      --ld-file /home/ref/1000G.EUR.QC.@.run4.ld \
      --trait1-file /home/sumstats/double_right.sumstats.gz \
      --trait2-file /home/sumstats/single_right_hip.sumstats.gz \
      --trait1-params-file /home/results/fit/double_right.fit.rep${REP}.json \
      --trait2-params-file /home/results/fit/single_right_hip.fit.rep${REP}.json \
      --extract /home/ref/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps \
      --chr2use 1-22 \
      --threads 24 \
      --out /home/results/fit/double_right_vs_single_right_hip.fit.rep${REP} \
      2>&1 | tee "$MIXER_DIR/results/logs/double_right_vs_single_right_hip.fit.rep${REP}.log"


    # ==================================================
    # 3. double_right test1
    # ==================================================

    docker run --rm \
      -v "$MIXER_DIR:/home" \
      -w /home \
      "$IMAGE" \
      python /tools/mixer/precimed/mixer.py test1 \
      --lib /tools/mixer/src/build/lib/libbgmg.so \
      --bim-file /home/ref/1000G.EUR.QC.@.bim \
      --ld-file /home/ref/1000G.EUR.QC.@.run4.ld \
      --trait1-file /home/sumstats/double_right.sumstats.gz \
      --load-params-file /home/results/fit/double_right.fit.rep${REP}.json \
      --chr2use 1-22 \
      --threads 24 \
      --out /home/results/test/double_right.test.rep${REP} \
      2>&1 | tee "$MIXER_DIR/results/logs/double_right.test.rep${REP}.log"


    # ==================================================
    # 4. double_right vs single_right_hip test2
    # ==================================================

    docker run --rm \
      -v "$MIXER_DIR:/home" \
      -w /home \
      "$IMAGE" \
      python /tools/mixer/precimed/mixer.py test2 \
      --lib /tools/mixer/src/build/lib/libbgmg.so \
      --bim-file /home/ref/1000G.EUR.QC.@.bim \
      --ld-file /home/ref/1000G.EUR.QC.@.run4.ld \
      --trait1-file /home/sumstats/double_right.sumstats.gz \
      --trait2-file /home/sumstats/single_right_hip.sumstats.gz \
      --load-params-file /home/results/fit/double_right_vs_single_right_hip.fit.rep${REP}.json \
      --chr2use 1-22 \
      --threads 24 \
      --out /home/results/test/double_right_vs_single_right_hip.test.rep${REP} \
      2>&1 | tee "$MIXER_DIR/results/logs/double_right_vs_single_right_hip.test.rep${REP}.log"


    # ==================================================
    # 5. double_right vs single_right_knee fit2
    # ==================================================

    docker run --rm \
      -v "$MIXER_DIR:/home" \
      -w /home \
      "$IMAGE" \
      python /tools/mixer/precimed/mixer.py fit2 \
      --lib /tools/mixer/src/build/lib/libbgmg.so \
      --bim-file /home/ref/1000G.EUR.QC.@.bim \
      --ld-file /home/ref/1000G.EUR.QC.@.run4.ld \
      --trait1-file /home/sumstats/double_right.sumstats.gz \
      --trait2-file /home/sumstats/single_right_knee.sumstats.gz \
      --trait1-params-file /home/results/fit/double_right.fit.rep${REP}.json \
      --trait2-params-file /home/results/fit/single_right_knee.fit.rep${REP}.json \
      --extract /home/ref/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${REP}.snps \
      --chr2use 1-22 \
      --threads 24 \
      --out /home/results/fit/double_right_vs_single_right_knee.fit.rep${REP} \
      2>&1 | tee "$MIXER_DIR/results/logs/double_right_vs_single_right_knee.fit.rep${REP}.log"


    # ==================================================
    # 6. double_right vs single_right_knee test2
    # ==================================================

    docker run --rm \
      -v "$MIXER_DIR:/home" \
      -w /home \
      "$IMAGE" \
      python /tools/mixer/precimed/mixer.py test2 \
      --lib /tools/mixer/src/build/lib/libbgmg.so \
      --bim-file /home/ref/1000G.EUR.QC.@.bim \
      --ld-file /home/ref/1000G.EUR.QC.@.run4.ld \
      --trait1-file /home/sumstats/double_right.sumstats.gz \
      --trait2-file /home/sumstats/single_right_knee.sumstats.gz \
      --load-params-file /home/results/fit/double_right_vs_single_right_knee.fit.rep${REP}.json \
      --chr2use 1-22 \
      --threads 24 \
      --out /home/results/test/double_right_vs_single_right_knee.test.rep${REP} \
      2>&1 | tee "$MIXER_DIR/results/logs/double_right_vs_single_right_knee.test.rep${REP}.log"


    echo "======================================================"
    echo "Replicate $REP completed successfully"
    echo "======================================================"

done