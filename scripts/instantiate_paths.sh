#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# DXA-aging-pipeline / scripts
# 把仓库里的 /path/to/... 占位符替换为本机实际路径（首次运行前执行一次）
# 用法: bash scripts/instantiate_paths.sh --project /data/proj --epi /data/epi \
#                                        --ref /data/ref --tools /opt/bio
# 注意: 只改写代码文件（.R/.sh/.py/.ipynb），不改 README 与 docs。
# ---------------------------------------------------------------------------
set -euo pipefail

PROJECT=""
EPI=""
REF=""
TOOLS=""
SHARED=""
EXTERNAL=""
RAWGENO=""
DRY=0

usage() {
    sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//'
    echo
    echo "  --project   Python 侧项目根      (默认 /path/to/project)"
    echo "  --epi       R 侧项目根           (默认 /path/to/epi)"
    echo "  --ref       共享参考数据根        (默认 /path/to/ref)"
    echo "  --tools     第三方工具根          (默认 /path/to/tools)"
    echo "  --shared    上游共享表根          (默认 /path/to/shared_project)"
    echo "  --external  早期项目代码根        (默认 /path/to/external_project)"
    echo "  --rawgeno   原始基因型根          (默认 /path/to/raw_genotype)"
    echo "  --dry-run   只打印将要替换的条数"
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --project)  PROJECT="$2";  shift 2 ;;
        --epi)      EPI="$2";      shift 2 ;;
        --ref)      REF="$2";      shift 2 ;;
        --tools)    TOOLS="$2";    shift 2 ;;
        --shared)   SHARED="$2";   shift 2 ;;
        --external) EXTERNAL="$2"; shift 2 ;;
        --rawgeno)  RAWGENO="$2";  shift 2 ;;
        --dry-run)  DRY=1;         shift ;;
        -h|--help)  usage ;;
        *) echo "未知参数: $1"; usage ;;
    esac
done

[[ -n "$PROJECT$EPI$REF$TOOLS$SHARED$EXTERNAL$RAWGENO" ]] || usage

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# 占位符 -> 实际路径（只替换已提供的）
PAIRS=()
[[ -n "$SHARED"   ]] && PAIRS+=("/path/to/shared_project|$SHARED")
[[ -n "$EXTERNAL" ]] && PAIRS+=("/path/to/external_project|$EXTERNAL")
[[ -n "$RAWGENO"  ]] && PAIRS+=("/path/to/raw_genotype|$RAWGENO")
[[ -n "$PROJECT"  ]] && PAIRS+=("/path/to/project|$PROJECT")
[[ -n "$EPI"      ]] && PAIRS+=("/path/to/epi|$EPI")
[[ -n "$REF"      ]] && PAIRS+=("/path/to/ref|$REF")
[[ -n "$TOOLS"    ]] && PAIRS+=("/path/to/tools|$TOOLS")

echo "仓库根目录: $ROOT"
echo "待替换规则: ${#PAIRS[@]} 条"

FILES=$(find "$ROOT" -type f \( -name "*.R" -o -name "*.sh" -o -name "*.py" -o -name "*.ipynb" \) \
        ! -path "*/scripts/instantiate_paths.sh" | sort)

total=0
for f in $FILES; do
    n=0
    for p in "${PAIRS[@]}"; do
        old="${p%%|*}"; new="${p##*|}"
        c=$(grep -oF -- "$old" "$f" 2>/dev/null | wc -l || true)
        n=$((n + c))
    done
    if [[ "$n" -eq 0 ]]; then
        continue
    fi
    total=$((total + n))
    if [[ "$DRY" -eq 1 ]]; then
        echo "  [DRY] $n 处  ${f#$ROOT/}"
    else
        cp -n "$f" "$f.bak" 2>/dev/null || true
        for p in "${PAIRS[@]}"; do
            old="${p%%|*}"; new="${p##*|}"
            # 用 | 作分隔符，路径中若含 & 或 / 也能安全替换
            sed -i.tmp "s|${old}|${new}|g" "$f"
            rm -f "$f.tmp"
        done
        echo "  [OK ] $n 处  ${f#$ROOT/}"
    fi
done

echo
if [[ "$DRY" -eq 1 ]]; then
    echo "共 $total 处待替换（dry-run，未改动任何文件）"
else
    echo "共替换 $total 处；原文件备份为 *.bak（已存在则不覆盖）"
    echo "如需回滚: find . -name '*.bak' -exec sh -c 'mv \"\$1\" \"\${1%.bak}\"' _ {} \\;"
fi
