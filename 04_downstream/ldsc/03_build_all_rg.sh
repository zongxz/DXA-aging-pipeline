# ---------------------------------------------------------------------------
# DXA-aging-pipeline / 04_downstream
# 从多个 *_rg.log 汇总出 all.rg 矩阵（供 LAVA 计算样本重叠）
# 注意: 文中 /path/to/... 为占位符，运行前请按本机实际路径替换（见 config/）。
# ---------------------------------------------------------------------------
FILES=($(ls *_rg.log))        # assuming the output format is [phenotype]_rg.log
N=$(echo ${#FILES[@]})          # and that all combinations of penotypes have been analysed

for I in ${FILES[@]}; do
        PHEN=$(echo $I | sed 's/_rg\.log//')
        
        # subset log files to relevant output
        tail -n$(($N+4)) $I | head -$((N+1)) > $PHEN.rg     # (adapt as necessary)
        
        # add to single data set
        if [[ $I == ${FILES[0]} ]]; then
            cat $PHEN.rg > all.rg       # only including the header for the first phenotypes
        else
            cat $PHEN.rg | sed '1d' >> all.rg
        fi
done