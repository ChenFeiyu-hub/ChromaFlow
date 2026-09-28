#!/usr/bin/env bash
# ============================================================================
# lib/common.sh — 公共库：路径解析 / 配置加载 / 日志
# ============================================================================
# 本文件属于「软件本体」，用户不应修改。
# 它负责把「代码位置」「项目配置」「输入/输出目录」三者解耦：
#   - CHROMAFLOW_HOME : 代码所在目录（由 bin/chromaflow 自动探测）
#   - 项目配置文件     : 用户为每个数据集准备一份，声明物种与 I/O 路径
#   - 所有输出路径     : 一律从 OUTPUT_DIR 派生，绝不写回代码目录
# ============================================================================

# ─────────────── 日志与工具检查 ───────────────
log_section() {
    echo ""
    echo "============================================"
    echo " $1"
    echo " $(date)"
    echo "============================================"
    echo ""
}

log_info()  { echo "[INFO]  $*"; }
log_warn()  { echo "[WARN]  $*" >&2; }
log_error() { echo "[ERROR] $*" >&2; }
die()       { log_error "$*"; exit 1; }

# 输出目录覆盖：若设置了 CHROMAFLOW_OUTDIR（CLI 全局 --outdir），各步用它作为
# 本步「输出」目录，替换默认值；输入仍从 config 的 RESULTS_DIR 解析。
# 用法: OUT="$(step_outdir "<默认输出目录>")"
step_outdir() { echo "${CHROMAFLOW_OUTDIR:-$1}"; }

# 从任意"基因列表/结果表"文件提取基因ID到每行一个：
#   - 有表头且含 gene_id/gene/genes/gene_name/id 等(不分大小写) → 取该列并跳过表头
#     （兼容 lookup 的 query_binding.tsv、candidate_targets.tsv 等多列带表头文件）
#   - 无可识别表头 → 视为无表头，取第一列（含首行）
# 用法: normalize_gene_list <in> <out>
normalize_gene_list() {
    python3 - "$1" "$2" <<'PYEOF'
import sys
inp, out = sys.argv[1], sys.argv[2]
lines = [l.rstrip("\n") for l in open(inp) if l.strip() and not l.lstrip().startswith("#")]
if not lines:
    open(out,"w").write(""); sys.exit(0)
delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else None)
def split(l): return l.split(delim) if delim else l.split()
first = split(lines[0])
aliases = {"gene_id","geneid","gene","genes","gene_name","genename","id","locus","locusname","name"}
hdr_idx = next((i for i,h in enumerate(first) if h.strip().strip('"').lower() in aliases), None)
genes=[]
body = lines[1:] if hdr_idx is not None else lines
col = hdr_idx if hdr_idx is not None else 0
for l in body:
    p = split(l)
    if col < len(p):
        g = p[col].strip().strip('"')
        if g: genes.append(g)
seen=set(); uniq=[g for g in genes if not (g in seen or seen.add(g))]
open(out,"w").write("\n".join(uniq)+("\n" if uniq else ""))
PYEOF
}

check_tool() {
    command -v "$1" &>/dev/null || { log_error "$1 is not installed or not on PATH"; return 1; }
}

# 可选功能开关：SKIP_ORGANELLE / SKIP_AGAT / SKIP_KEGG ... (config 或环境变量)
# 用法: if is_skipped organelle; then ... fi
is_skipped() {
    local var="SKIP_$(echo "$1" | tr '[:lower:]' '[:upper:]')" val
    val="$(eval "echo \"\${${var}:-}\"")"
    [[ "${val}" == "true" || "${val}" == "1" || "${val}" == "yes" ]]
}

# ─────────────── 压缩识别与按需解压 ───────────────
# 按「文件头魔数」判断，而不是看扩展名：
#   · 已解压的 .gff3/.fa 直接原样使用，不做无谓拷贝
#   · 命名不规范的文件（叫 .gff3 实为 gzip、Phytozome 下载常见）也能正确识别
#   · 支持 gzip/bgzip、bzip2、xz、zstd、zip
detect_compression() {
    local f="${1:-}" magic
    [[ -s "${f}" ]] || { echo "none"; return 1; }
    magic=$(head -c 6 "${f}" 2>/dev/null | od -An -v -tx1 | tr -d ' \n')
    case "${magic}" in
        1f8b*)       echo gzip  ;;   # gzip / bgzip
        425a68*)     echo bzip2 ;;
        fd377a585a*) echo xz    ;;
        28b52ffd*)   echo zstd  ;;
        504b0304*)   echo zip   ;;
        *)           echo none  ;;
    esac
}

_decompress_cmd() {
    case "$1" in
        gzip)  echo "gunzip -c"  ;;
        bzip2) echo "bzip2 -dc"  ;;
        xz)    echo "xz -dc"     ;;
        zstd)  echo "zstd -dc"   ;;
        zip)   echo "unzip -p"   ;;
        *)     echo "cat"        ;;
    esac
}

# ensure_plain_file <src> <dest_dir> [目标文件名]
#   未压缩 → 回显源路径本身（零拷贝）
#   已压缩 → 解压到 dest_dir 并回显新路径（已存在且比源新则直接复用）
#   FORCE_DECOMPRESS=true 可强制重解压
# 注意：本函数只把路径写到 stdout，提示信息一律走 stderr，便于 $(...) 取值。
ensure_plain_file() {
    local src="${1:-}" dest_dir="${2:-}" want="${3:-}" comp base out cmd tool
    [[ -s "${src}" ]] || return 1
    comp="$(detect_compression "${src}")"
    if [[ "${comp}" == "none" ]]; then
        echo "${src}"; return 0
    fi

    tool="$(_decompress_cmd "${comp}" | awk '{print $1}')"
    command -v "${tool}" &>/dev/null || {
        log_error "File uses ${comp} compression; decompression tool ${tool} is missing: ${src}" >&2
        return 1
    }

    if [[ -z "${want}" ]]; then
        base="$(basename "${src}")"
        case "${base}" in
            *.gz|*.bz2|*.xz|*.zst|*.zip) want="${base%.*}" ;;
            # 扩展名没体现压缩（如 .gff3 实为 gzip）→ 插入 .plain 避免覆盖源文件
            *.*) want="${base%.*}.plain.${base##*.}" ;;
            *)   want="${base}.plain" ;;
        esac
    fi
    mkdir -p "${dest_dir}" 2>/dev/null || true
    out="${dest_dir}/${want}"
    [[ "${out}" == "${src}" ]] && out="${out}.plain"

    if [[ -s "${out}" && "${out}" -nt "${src}" && "${FORCE_DECOMPRESS:-false}" != "true" ]]; then
        log_info "Reusing decompressed file: ${out}" >&2
        echo "${out}"; return 0
    fi
    log_info "Detected ${comp} compression; decompressing: $(basename "${src}") → ${out}" >&2
    cmd="$(_decompress_cmd "${comp}")"
    ${cmd} "${src}" > "${out}" || { log_error "Decompression failed: ${src}" >&2; return 1; }
    echo "${out}"
}

# ─────────────── 实测标量值（与 reference.manifest 并列）───────────────
# reference.manifest 记「文件在哪」；derived.env 记「量出来的数」
#   (MACS_GSIZE / READ_LENGTH / GENOME_GC ...)
# 一律不回写用户的 chromaflow.config：config 是「你的声明」，derived.env 是
# 「脚本的观测」，两者分离才不会互相覆盖，也方便随 outputs 一起丢弃重算。
# 优先级: config 显式赋值 > derived.env > 现场测量
derived_values_path() { echo "${OUTPUT_DIR}/.chromaflow/derived.env"; }

# derived_set <KEY> <VALUE>   由测量它的步骤调用
derived_set() {
    local key="$1" val="$2" df
    [[ -n "${key}" ]] || return 0
    df="$(derived_values_path)"
    mkdir -p "$(dirname "${df}")" 2>/dev/null || return 0
    if [[ -f "${df}" ]]; then grep -v "^${key}=" "${df}" > "${df}.tmp" 2>/dev/null || : > "${df}.tmp"
    else : > "${df}.tmp"; fi
    printf '%s=%s\n' "${key}" "${val}" >> "${df}.tmp"
    mv "${df}.tmp" "${df}"
}

derived_get() {
    local df; df="$(derived_values_path)"
    [[ -f "${df}" ]] || return 1
    local line; line="$(grep -m1 "^$1=" "${df}" 2>/dev/null)" || return 1
    echo "${line#*=}"
}

# 只填「当前为空」的变量（config 显式声明始终优先）
load_derived_values() {
    local df k v cur; df="$(derived_values_path)"
    [[ -f "${df}" ]] || return 0
    while IFS='=' read -r k v; do
        [[ -n "${k}" && "${k}" != \#* ]] || continue
        [[ "${k}" =~ ^[A-Z][A-Z0-9_]*$ ]] || continue
        cur="$(eval "echo \"\${${k}:-}\"")"
        [[ -n "${cur}" ]] && continue
        declare -g "${k}=${v}"
        CHROMAFLOW_FROM_DERIVED="${CHROMAFLOW_FROM_DERIVED:-} ${k}"
    done < "${df}"
}

# 某个值是不是从 derived.env 来的（doctor 用来标注来源）
value_origin() {
    local key="$1"
    [[ " ${CHROMAFLOW_FROM_DERIVED:-} " == *" ${key} "* ]] && { echo "measured"; return; }
    [[ -n "$(eval "echo \"\${${key}:-}\"")" ]] && { echo "config"; return; }
    echo "unset"
}

# 文件签名：路径|大小|参数，用于判断实测值是否还对得上当前输入
file_signature() {
    local f="${1:-}" extra="${2:-}"
    [[ -s "${f}" ]] || { echo ""; return 1; }
    echo "$(basename "${f}")|$(stat -c%s "${f}" 2>/dev/null)|${extra}"
}

# ─────────────── 实测：有效基因组大小 / GC / 读长 ───────────────
# 有效基因组大小 (MACS 的 -g)：MACS 用它算 Poisson 背景 λ = N·d/gsize。
# 首选 khmer 的 unique-kmers.py（deepTools/nf-core 同款做法：统计长度为读长的
# 唯一 k-mer 数，即「理论上可唯一比对到的碱基数」）；khmer 不可用时退到
# 非 N 碱基数（偏大 ~5-10%，作为上界仍可用）。
# measure_effective_genome_size <fasta> <k> → "数值<空格>方法"，失败返回非 0
# （回显方法名而不是设全局变量：本函数总在 $(...) 子 shell 里被调用）
measure_effective_genome_size() {
    local fa="${1:-}" k="${2:-150}" val=""
    [[ -s "${fa}" ]] || return 1
    if command -v unique-kmers.py &>/dev/null; then
        val=$(unique-kmers.py -k "${k}" "${fa}" 2>&1 \
              | grep -oP 'Total estimated number of unique .*-mers:\s*\K[0-9]+' | tail -1)
        if [[ -n "${val}" ]]; then
            echo "${val} unique-kmers(k=${k})"; return 0
        fi
    fi
    # 兜底：非 N 碱基计数（一次流式扫描，无额外依赖）
    val=$(awk '!/^>/ { s=$0; gsub(/[NnXx]/,"",s); n += length(s) } END { print n+0 }' "${fa}")
    [[ "${val:-0}" -gt 0 ]] || return 1
    echo "${val} non-N-bases(upper-bound)"; return 0
}

# 基因组 GC%（一次扫描，顺带给 04 的背景模型做参照）
measure_genome_gc() {
    local fa="${1:-}"
    [[ -s "${fa}" ]] || return 1
    awk '!/^>/ {
            s = toupper($0)
            gc += gsub(/[GC]/, "", s)
            at += gsub(/[AT]/, "", s)
         }
         END { if (gc+at > 0) printf "%.1f", 100*gc/(gc+at); else print "" }' "${fa}"
}

# 从 fastq 抽样实测读长（取抽样中的最大值：trimming 只会变短不会变长）
# detect_read_length <fastq> [抽样 read 数，默认 2000]
detect_read_length() {
    local fq="${1:-}" n="${2:-2000}" comp cmd len
    [[ -s "${fq}" ]] || return 1
    comp="$(detect_compression "${fq}")"
    cmd="$(_decompress_cmd "${comp}")"
    len=$(${cmd} "${fq}" 2>/dev/null | head -n $(( n * 4 )) \
          | awk 'NR%4==2 { if (length($0) > m) m = length($0) } END { print m+0 }')
    [[ "${len:-0}" -gt 0 ]] || return 1
    echo "${len}"
}

# 可选：把实测值写回 config（默认不做，仅在 --write-config 时调用）
# 用带标记的块整体替换，保留用户手写的内容与注释；重复执行不会累积。
write_measured_block_to_config() {
    local cfg="${CHROMAFLOW_CONFIG:-}" df tmp
    [[ -f "${cfg}" ]] || { log_warn "Config not found; skipping write-back"; return 1; }
    df="$(derived_values_path)"
    [[ -f "${df}" ]] || { log_warn "No measured values available; skipping write-back"; return 1; }
    tmp="${cfg}.chromaflow.tmp"
    awk '/^# >>> chromaflow measured >>>/{skip=1} !skip{print} /^# <<< chromaflow measured <<</{skip=0}' \
        "${cfg}" > "${tmp}"
    {
        echo ""
        echo "# >>> chromaflow measured >>>"
        echo "# Generated by chromaflow reference --write-config; replaced on rerun."
        echo "# Without --write-config these values remain in .chromaflow/derived.env."
        echo "# Both locations are supported; this copy is for project archiving."
        # 只写真正的参数（白名单），并且加引号：方法名等元数据含括号，
        # 不加引号会让 config 变成非法 bash。
        local k v
        for k in MACS_GSIZE READ_LENGTH GENOME_GC; do
            v="$(grep -m1 "^${k}=" "${df}" 2>/dev/null | cut -d= -f2-)"
            [[ -n "${v}" ]] && printf '%s="%s"\n' "${k}" "${v}"
        done
        echo "# <<< chromaflow measured <<<"
    } >> "${tmp}"
    cp "${cfg}" "${cfg}.bak.$(date +%Y%m%d_%H%M%S)"
    mv "${tmp}" "${cfg}"
    log_info "Updated ${cfg} (original backed up as ${cfg}.bak.*)"
}

# ─────────────── nf-core conda 模式的已知冲突修补 ───────────────
# 背景：nf-core 各 process 的 conda spec 只锁工具版本、不锁依赖版本。conda 解析
# 时会拉最新的 numpy(2.x)/python(3.13)，而 deeptools 3.5.x 仍在用 numpy 2.0
# 已删除的 np.NAN，于是 DEEPTOOLS_COMPUTEMATRIX 在 conda 模式下必崩：
#     AttributeError: module 'numpy' has no attribute 'NAN'
# 容器模式没这问题（镜像里的依赖是构建时固定好的），所以这是 conda 模式专属。
# 做法：生成一份 nextflow 配置用 -c 传入，为受影响的 process 覆盖 conda spec，
# 补上 numpy<2 的约束。不改用户 config，也不碰 nf-core 源码。

# 从已下载的 pipeline assets 里读出 nf-core 锁定的工具版本（读不到返回空）
_nfcore_module_version() {   # _nfcore_module_version <pipeline> <tool>
    local pipeline="$1" tool="$2" base yml v
    base="${NXF_ASSETS:-${HOME}/.nextflow/assets}/nf-core/${pipeline}"
    [[ -d "${base}" ]] || return 1
    while read -r yml; do
        v="$(grep -oP "${tool}\s*[=:]+\s*\K[0-9][0-9a-zA-Z.]*" "${yml}" 2>/dev/null | head -1)"
        [[ -n "${v}" ]] && { echo "${v}"; return 0; }
    done < <(find "${base}/modules" -path "*${tool}*" -name "environment.yml" 2>/dev/null)
    return 1
}

# 半成品 conda 环境清理
# nextflow 只按目录是否存在判断环境「已建好」。一旦 createTimeout 到时把
# conda 杀掉，就会留下一个不完整的 env 目录；下次运行 nextflow 直接复用它，
# 于是在 process 里表现为「命令找不到」（我们踩过的 run_spp.R 就是这样来的）。
# 判据：完整环境一定有 conda-meta/history。
clean_partial_conda_envs() {   # clean_partial_conda_envs [envs_dir]
    local dir="${1:-${OUTPUT_DIR}/.conda_envs}" env n=0
    [[ -d "${dir}" ]] || return 0
    while read -r env; do
        [[ -n "${env}" ]] || continue
        [[ -f "${env}/conda-meta/history" ]] && continue
        log_warn "Removing incomplete Conda environment before rebuilding: $(basename "${env}")"
        rm -rf "${env}"; n=$((n+1))
    done < <(find "${dir}" -mindepth 1 -maxdepth 1 -type d -name 'env-*' 2>/dev/null)
    [[ ${n} -gt 0 ]] && echo "  Removed ${n} incomplete environments"
    return 0
}

# 生成 conda 模式的补充配置，回显其路径
# with_pins=false 时只写超时等通用设置，不注入依赖版本约束
write_nfcore_conda_patch() {   # write_nfcore_conda_patch <pipeline> [with_pins]
    local pipeline="${1:-chipseq}" with_pins="${2:-true}" cfg dt_ver dt_spec
    cfg="${OUTPUT_DIR}/.chromaflow/nfcore_conda_patch.config"
    mkdir -p "$(dirname "${cfg}")" 2>/dev/null || return 1

    # 工具版本优先跟随 nf-core 的锁定值；读不到就不锁版本（只锁出问题的依赖）
    if dt_ver="$(_nfcore_module_version "${pipeline}" deeptools)"; then
        dt_spec="bioconda::deeptools=${dt_ver}"
    else
        dt_spec="bioconda::deeptools"
    fi
    local mq_ver mq_spec
    if mq_ver="$(_nfcore_module_version "${pipeline}" multiqc)"; then
        mq_spec="bioconda::multiqc=${mq_ver}"
    else
        mq_spec="bioconda::multiqc"
    fi

    # 用 mamba 建环境：mamba 的下载器有重试/并行/断点能力，conda 自带的
    # curl 调用在慢链路上极易断（curl: (18) end of response ... bytes missing）。
    # ★ 只能通过配置里的 conda.useMamba 打开；NXF_CONDA_USE_MAMBA 这个环境变量
    #   并不存在（旧脚本 export 了它，其实一直是 conda 在跑）。
    local mamba_line=""
    if command -v mamba &>/dev/null; then
        mamba_line="    useMamba = true"
    elif command -v micromamba &>/dev/null; then
        mamba_line="    useMicromamba = true"
    fi

    cat > "${cfg}" << EOF

conda {
${mamba_line}
    createTimeout = '${NFCORE_CONDA_TIMEOUT:-3 h}'
}
EOF

    if [[ "${with_pins}" == "true" ]]; then
        cat >> "${cfg}" << EOF

process {
    withName: '.*:?DEEPTOOLS_.*' {
        conda = "${dt_spec} conda-forge::numpy=1.26 conda-forge::matplotlib-base"
    }
    withName: '.*:?MULTIQC' {
        conda = "${mq_spec} conda-forge::python=3.12 conda-forge::python-kaleido=0.2.1"
    }
}
EOF
    fi
    echo "${cfg}"
}

# ─────────────── nf-core 执行模式 (profile) ───────────────
# 解析优先级: --profile 命令行 (CHROMAFLOW_PROFILE) > config 的 NFCORE_PROFILE > singularity
resolve_nfcore_profile() {
    echo "${CHROMAFLOW_PROFILE:-${NFCORE_PROFILE:-singularity}}"
}

# 判断 profile 属于哪条执行路径: container | conda | other
# 支持逗号组合(如 "singularity,test")，按其中出现的引擎关键字判断。
nfcore_profile_kind() {
    local p="${1:-$(resolve_nfcore_profile)}"
    case ",${p}," in
        *,conda,*|*,mamba,*)                       echo conda ;;
        *singularity*|*apptainer*|*docker*|*podman*|*charliecloud*|*shifter*) echo container ;;
        *) echo other ;;
    esac
}

# profile 需要的容器引擎可执行文件（container 类才有意义）
nfcore_container_engine() {
    local p="${1:-$(resolve_nfcore_profile)}"
    case "${p}" in
        *apptainer*)   echo apptainer   ;;
        *singularity*) echo singularity ;;
        *podman*)      echo podman      ;;
        *docker*)      echo docker      ;;
        *)             echo ""          ;;
    esac
}

# ─────────────── 参考文件：注册表 + 自动识别 ───────────────
# 设计原则（避免「把脚本输出写进 config」）：
#   config 只声明你手里「本来就有」的源文件（源 fasta / 源 gff3 / 注释表 / motif 库）；
#   由 chromaflow 生成的派生文件（combined.fa / .fai / chrom.sizes / .gtf / *.agat.gff3）
#   由生成它们的步骤写进 manifest，其余步骤靠 manifest + 目录扫描自动找到。
# 解析优先级: config 显式赋值 > manifest 注册 > 目录自动识别 > 兜底命名(可能尚不存在)

# 派生文件落盘目录：优先 REFERENCE_DIR；只读(如共享参考库)时退到 OUTPUT_DIR/reference
derived_ref_dir() {
    if [[ -n "${REFERENCE_DIR:-}" && -d "${REFERENCE_DIR}" && -w "${REFERENCE_DIR}" ]]; then
        echo "${REFERENCE_DIR}"
    elif [[ -n "${REFERENCE_DIR:-}" && ! -e "${REFERENCE_DIR}" ]] && mkdir -p "${REFERENCE_DIR}" 2>/dev/null; then
        echo "${REFERENCE_DIR}"
    else
        mkdir -p "${OUTPUT_DIR}/reference" 2>/dev/null || true
        echo "${OUTPUT_DIR}/reference"
    fi
}

reference_manifest_path() { echo "${OUTPUT_DIR}/.chromaflow/reference.manifest"; }

# manifest_register <KEY> <PATH>  由 00 / 00b 等「生产者」调用
manifest_register() {
    local key="$1" val="$2" mf
    [[ -n "${key}" && -n "${val}" ]] || return 0
    mf="$(reference_manifest_path)"
    mkdir -p "$(dirname "${mf}")" 2>/dev/null || return 0
    if [[ -f "${mf}" ]]; then grep -v "^${key}=" "${mf}" > "${mf}.tmp" 2>/dev/null || : > "${mf}.tmp"
    else : > "${mf}.tmp"; fi
    printf '%s=%s\n' "${key}" "${val}" >> "${mf}.tmp"
    mv "${mf}.tmp" "${mf}"
    log_info "Registered ${key} → ${val}"
}

# 读取 manifest：只填「当前为空」的变量（config 显式声明始终优先）
load_reference_manifest() {
    local mf k v cur; mf="$(reference_manifest_path)"
    [[ -f "${mf}" ]] || return 0
    while IFS='=' read -r k v; do
        [[ -n "${k}" && "${k}" != \#* ]] || continue
        [[ "${k}" =~ ^[A-Z][A-Z0-9_]*$ ]] || continue
        [[ -e "${v}" ]] || continue                       # 文件已删/搬走 → 交给automatic detection
        cur="$(eval "echo \"\${${k}:-}\"")"
        [[ -n "${cur}" ]] && continue
        declare -g "${k}=${v}"
    done < "${mf}"
}

# 在参考目录（及派生目录）里按模式找最新的一个文件
_ref_search_dirs() {
    local d
    { [[ -n "${REFERENCE_DIR:-}" ]] && echo "${REFERENCE_DIR}"
      [[ -n "${OUTPUT_DIR:-}"    ]] && echo "${OUTPUT_DIR}/reference"
      [[ -n "${INPUT_DIR:-}"     ]] && echo "${INPUT_DIR}/reference"
    } 2>/dev/null | awk 'NF && !seen[$0]++'
}

_find_ref() {   # _find_ref <glob> [glob...] → 第一个命中模式、且未压缩的最新文件
    local pat dir f
    for pat in "$@"; do
        while read -r dir; do
            [[ -d "${dir}" ]] || continue
            while read -r f; do
                [[ -n "${f}" ]] || continue
                # 跳过压缩文件：gff3/gtf/chrom.sizes 等下游要求纯文本
                [[ "$(detect_compression "${f}")" == "none" ]] || continue
                echo "${f}"; return 0
            done < <(find -L "${dir}" -maxdepth 2 -type f -name "${pat}" -printf '%T@\t%p\n' 2>/dev/null \
                     | sort -rn | cut -f2-)
        done < <(_ref_search_dirs)
    done
    return 1
}

# 自动识别派生参考文件（只填空变量）
autodetect_reference_files() {
    local f cands=()
    # 基因组 FASTA：按「combined > 已建索引 > nuclear」打分挑选
    if [[ -z "${GENOME_COMBINED:-}" ]]; then
        local best="" best_score=-1 score dir x
        while read -r dir; do
            [[ -d "${dir}" ]] || continue
            while read -r x; do [[ -n "${x}" ]] && cands+=("${x}"); done < <(
                find -L "${dir}" -maxdepth 2 -type f \( -name "*.fa" -o -name "*.fasta" \) 2>/dev/null | sort)
        done < <(_ref_search_dirs)
        for x in ${cands[@]+"${cands[@]}"}; do
            [[ "$(detect_compression "${x}")" == "none" ]] || continue   # .fa 里若是压缩内容则跳过
            score=0
            [[ "${x}" == *combined* ]] && score=$((score + 4))
            [[ -f "${x}.fai"       ]] && score=$((score + 2))
            [[ "${x}" == *_nuclear.fa ]] && score=$((score + 1))
            if [[ ${score} -gt ${best_score} ]]; then best="${x}"; best_score=${score}; fi
        done
        [[ -n "${best}" ]] && GENOME_COMBINED="${best}"
    fi
    # GFF3：AGAT/规范化产物优先，其次任意 gff3（→ 未跑 00b 也能注释）
    [[ -z "${GENOME_GFF3:-}" ]] && f=$(_find_ref "*.normalized.gff3" "*.agat.gff3" "*.gff3" "*.gff") \
        && GENOME_GFF3="${f}"
    [[ -z "${GENOME_GTF:-}"  ]] && f=$(_find_ref "*.gtf")                     && GENOME_GTF="${f}"
    [[ -z "${CHROM_SIZES:-}" ]] && f=$(_find_ref "*.chrom.sizes" "*.chromsizes" "*.genome") \
        && CHROM_SIZES="${f}"
    [[ -z "${BLACKLIST_BED:-}"   ]] && f=$(_find_ref "*blacklist*.bed")       && BLACKLIST_BED="${f}"
    [[ -z "${ANNOTATION_INFO:-}" ]] && f=$(_find_ref "*annotation_info*.txt" "*annotation_info*.tsv") \
        && ANNOTATION_INFO="${f}"
    [[ -z "${MEME_DB:-}" ]] && f=$(_find_ref "*.meme")                        && MEME_DB="${f}"
    return 0
}

# chrom.sizes 按需生成（缺就从 .fai 现算，不必让用户先跑 00 或写进 config）
ensure_chrom_sizes() {
    if [[ -s "${CHROM_SIZES:-}" ]]; then echo "${CHROM_SIZES}"; return 0; fi
    local fai="${GENOME_COMBINED:-}.fai" out
    [[ -s "${GENOME_COMBINED:-}" ]] || { echo ""; return 1; }
    if [[ ! -s "${fai}" ]] && command -v samtools &>/dev/null; then
        samtools faidx "${GENOME_COMBINED}" 2>/dev/null || true
    fi
    [[ -s "${fai}" ]] || { echo ""; return 1; }
    out="$(derived_ref_dir)/$(basename "${GENOME_COMBINED%.*}").chrom.sizes"
    cut -f1,2 "${fai}" > "${out}" || { echo ""; return 1; }
    CHROM_SIZES="${out}"; manifest_register CHROM_SIZES "${out}" >/dev/null
    echo "${CHROM_SIZES}"
}

# 细胞器名列表（ORGANELLE_SOURCES 未定义/为空时安全返回空）
organelle_names() {
    local e
    [[ ${#ORGANELLE_SOURCES[@]} -eq 0 ]] && return 0
    for e in "${ORGANELLE_SOURCES[@]}"; do [[ -n "${e}" ]] && echo "${e%%:*}"; done
}

# ─────────────── 配置加载 ───────────────
# 用法: load_config <config_path>
# 解析优先级：显式 -c 参数 > 环境变量 CHROMAFLOW_CONFIG > ./chromaflow.config
resolve_config() {
    local explicit="${1:-}"
    if [[ -n "${explicit}" ]]; then
        echo "${explicit}"; return
    fi
    if [[ -n "${CHROMAFLOW_CONFIG:-}" ]]; then
        echo "${CHROMAFLOW_CONFIG}"; return
    fi
    if [[ -f "./chromaflow.config" ]]; then
        echo "$(pwd)/chromaflow.config"; return
    fi
    die "Project config not found. Use -c <config>, or add chromaflow.config to the project directory."
}

load_config() {
    local cfg
    cfg="$(resolve_config "${1:-}")"
    [[ -f "${cfg}" ]] || die "Config file does not exist: ${cfg}"
    CHROMAFLOW_CONFIG="$(cd "$(dirname "${cfg}")" && pwd)/$(basename "${cfg}")"

    # 两遍 source 解决「配置引用派生路径」的先后问题：
    #   1) 先取锚点（PROJECT_DIR 等），允许 REFERENCE_DIR-相对行暂时展开错误
    #   2) derive_paths 定义 REFERENCE_DIR/OUTPUT_DIR/... 
    #   3) 再 source 一遍：此时 ${REFERENCE_DIR}/... 等能正确展开
    set +u
    # shellcheck disable=SC1090
    source "${CHROMAFLOW_CONFIG}"
    derive_paths
    # shellcheck disable=SC1090
    source "${CHROMAFLOW_CONFIG}"
    derive_paths
    set -u
    validate_config
    # CLI --peakset takes precedence; project default is opt-in.
    export CHROMAFLOW_PEAKSET="${CHROMAFLOW_PEAKSET:-${DEFAULT_PEAKSET:-idr}}"
    export TRACK_GENE_HEIGHT="${TRACK_GENE_HEIGHT:-auto}"
    export TRACK_MOTIF_HEIGHT="${TRACK_MOTIF_HEIGHT:-auto}"
}

# ─────────────── 路径派生：一切从 config 的 3 个锚点展开 ───────────────
# config 只需给出：PROJECT_DIR（可选）、INPUT_DIR、OUTPUT_DIR、REFERENCE_DIR
# 其余全部在此派生，脚本正文不再拼绝对路径。
derive_paths() {
    # 锚点默认值：若 config 只给了 PROJECT_DIR，则 inputs/outputs 落在其下
    : "${PROJECT_DIR:=$(pwd)}"
    : "${INPUT_DIR:=${PROJECT_DIR}/inputs}"
    : "${OUTPUT_DIR:=${PROJECT_DIR}/outputs}"
    : "${REFERENCE_DIR:=${INPUT_DIR}/reference}"
    : "${RAW_DATA_DIR:=${INPUT_DIR}/raw}"

    # ---- 输出侧（全部在 OUTPUT_DIR 下，与代码/输入彻底隔离）----
    RESULTS_DIR="${OUTPUT_DIR}/results"
    LOG_DIR="${OUTPUT_DIR}/logs"
    SAMPLESHEET_DIR="${OUTPUT_DIR}/samplesheets"
    WORK_DIR="${OUTPUT_DIR}/work"               # nextflow work
    NXF_SINGULARITY_CACHEDIR="${OUTPUT_DIR}/.singularity_cache"
    export NXF_SINGULARITY_CACHEDIR SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR}"
    # Apptainer 1.x 读的是自己的变量名，一并指向同一缓存，避免镜像被下载两份
    export NXF_APPTAINER_CACHEDIR="${NXF_SINGULARITY_CACHEDIR}"
    export APPTAINER_CACHEDIR="${NXF_SINGULARITY_CACHEDIR}"

    # ---- 02 nf-core 运行参数（config 未写时的兜底，避免 set -u 下未定义变量报错）----
    # NFCORE_PROFILE 决定 nextflow 的 -profile，也决定 `chromaflow run` 走哪条执行路径：
    #   singularity / apptainer / docker / podman → 容器模式 (02_run_nfcore.sh)
    #   conda / mamba                             → conda 模式 (02_run_nfcore_conda.sh)
    : "${NFCORE_PROFILE:=singularity}"
    : "${NFCORE_CHIPSEQ_VERSION:=2.1.0}"
    : "${NFCORE_ATACSEQ_VERSION:=2.1.2}"
    : "${MAX_CPUS:=16}"
    : "${MAX_MEMORY:=64.GB}"
    : "${MAX_TIME:=72.h}"
    # 留空是正常状态：00 会用 khmer 实测并写进 derived.env；
    # 实在测不到时 02 会改为把 --read_length 交给 nf-core，由它内部的
    # KHMER_UNIQUEKMERS 自己算。两条路都不需要用户手填。
    : "${MACS_GSIZE:=}"

    # ---- 参考文件解析 ----
    # 优先级: config 显式赋值 > manifest(由 00/00b 登记) > 目录自动识别 > 兜底命名
    # config 只需声明「源文件」；派生文件不必写进 config（写了也依然生效，向后兼容）。

    # 源文件（用户手里本来就有的东西）。兼容旧变量名 PHYTOZOME_*/GENOME_GFF3_RAW。
    : "${GENOME_FASTA_SOURCE:=${PHYTOZOME_GENOME_GZ:-}}"                      # 源核基因组 .fa/.fa.gz
    : "${GENOME_GFF3_SOURCE:=${PHYTOZOME_GFF3_GZ:-${GENOME_GFF3_RAW:-}}}"     # 源注释 .gff3/.gff3.gz

    # 细胞器：未定义时给一个空数组，保证 set -u 下所有 ${ORGANELLE_SOURCES[@]} 安全
    declare -p ORGANELLE_SOURCES >/dev/null 2>&1 || ORGANELLE_SOURCES=()

    load_reference_manifest        # 步骤Registered过的派生文件
    load_derived_values            # 步骤已measured过的标量 (MACS_GSIZE / READ_LENGTH / GENOME_GC)
    autodetect_reference_files     # 目录扫描（*combined*.fa / *.agat.gff3 / *.gtf / *.chrom.sizes ...）

    # 兜底命名（文件可能尚未生成；00 会按这些名字产出并登记）
    : "${GENOME_NUCLEAR:=$(derived_ref_dir)/${SPECIES_SHORT}_nuclear.fa}"
    : "${GENOME_COMBINED:=$(derived_ref_dir)/${SPECIES_SHORT}_combined.fa}"
    : "${GENOME_GFF3:=$(derived_ref_dir)/${SPECIES_SHORT}_genes.gff3}"
    : "${GENOME_GFF3_RAW:=${GENOME_GFF3_SOURCE}}"                             # 00b 的默认输入
    : "${GENOME_GTF:=$(derived_ref_dir)/${SPECIES_SHORT}_genes.gtf}"
    : "${CHROM_SIZES:=$(derived_ref_dir)/${SPECIES_SHORT}_combined.chrom.sizes}"
    : "${BLACKLIST_BED:=$(derived_ref_dir)/${SPECIES_SHORT}_blacklist.bed}"
    : "${ANNOTATION_INFO:=}"                    # Phytozome annotation_info(含GO)，可留空
    : "${MEME_DB:=}"                            # motif 库（04 tomtom 用），可留空

    # 细胞器 fasta（由 ORGANELLE_SOURCES 派生；如已在 config 显式给出则不覆盖）
    if [[ ${#ORGANELLE_SOURCES[@]} -gt 0 ]]; then
        for _entry in "${ORGANELLE_SOURCES[@]}"; do
            [[ -n "${_entry}" ]] || continue
            local _name="${_entry%%:*}"
            local _upper; _upper=$(echo "${_name}" | tr '[:lower:]' '[:upper:]')
            [[ -z "$(eval echo "\${GENOME_${_upper}:-}")" ]] && \
                declare -g "GENOME_${_upper}=$(derived_ref_dir)/${_name}.fa"
        done
    fi

    # 测序读长（02 传给 nf-core；00 用于 unique-kmers 提示）
    : "${READ_LENGTH:=}"            # 留空=由 01 从 fastq measured（见 derived.env）

    # 下游默认参数（可被 config 覆盖）
    : "${PROMOTER_UP:=1000}"        # TSS 上游
    : "${PROMOTER_DOWN:=200}"       # TSS 下游
    : "${TARGET_DIST_KB:=3}"        # 距离过滤Target genes阈值 (kb)
    : "${MAKE_PLOTS:=true}"         # 是否用 matplotlib 出图

    # deepTools 信号可视化（step 06）
    : "${DEEPT_BIN_SIZE:=10}"       # computeMatrix binSize
    : "${SIGNAL_PEAK_WIN:=2000}"    # peak summit 两侧窗口 (bp)
    : "${SIGNAL_TSS_WIN:=3000}"     # TSS 两侧窗口 (bp)
    : "${THREADS:=${MAX_CPUS:-8}}"  # deepTools 并行核数

    # 单基因浏览器轨道图（step 07, pyGenomeTracks）
    : "${TRACK_FLANK:=3000}"        # 基因两侧延伸 (bp)
    : "${TRACK_MAX:=}"              # bigWig 轨道Shared y-axis maximum（空=自动按区域计算）

    # ChIPseeker 注释(05) 与 motif 整合(06)
    : "${KEGG_ORGANISM:=}"          # KEGG 物种码（衣藻=cre；空=跳过 KEGG）
    : "${FIMO_PVAL:=1e-4}"          # FIMO motif 命中 p 值阈值（匹配质量，非富集）
    : "${MEME_SEARCHSIZE:=0}"       # MEME 参与搜索的序列字符数上限(0=全部；序列多时可设值限速)
    # 候选"合适位置"特征集（逗号分隔）。可选值见 06 输出的 feature 列：
    #   Promoter / 5' UTR / 3' UTR / 1st Exon / Other Exon / 1st Intron / Other Intron / Downstream / Distal Intergenic
    : "${TARGET_FEATURES:=Promoter}"

    # 样本表与 FASTQ 配对（01）
    : "${SAMPLES_FILE:=}"           # 留空=自动在项目目录下找 samples[_type].txt/tsv/csv
    : "${FASTQ_PAIR_PATTERN:=}"     # 留空=automatic detection R1/R2 记号；也可写 "_good_{R}.fq.gz" 这类模板
    : "${FASTQ_SINGLE_END:=false}"  # true=forced single-end，不做Pair

    TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
}

# 自动定位样本表：显式 > config SAMPLES_FILE > 项目/输入目录下的常见命名
# 用法: resolve_samples_file <pipeline> [显式路径]
resolve_samples_file() {
    local pipeline="${1:-chipseq}" explicit="${2:-}" f d
    if [[ -n "${explicit}" ]]; then echo "${explicit}"; return 0; fi
    if [[ -n "${SAMPLES_FILE:-}" ]]; then echo "${SAMPLES_FILE}"; return 0; fi
    for d in "${PROJECT_DIR}" "${INPUT_DIR}" "${PROJECT_DIR}/samples" "$(pwd)"; do
        [[ -d "${d}" ]] || continue
        for f in "samples_${pipeline}.txt" "samples_${pipeline}.tsv" "samples_${pipeline}.csv" \
                 "samples.txt" "samples.tsv" "samples.csv"; do
            [[ -f "${d}/${f}" ]] && { echo "${d}/${f}"; return 0; }
        done
    done
    echo "${PROJECT_DIR}/samples_${pipeline}.txt"   # 兜底（用于报错信息）
    return 1
}

validate_config() {
    [[ -n "${SPECIES_SHORT:-}" ]] || die "Config is missing SPECIES_SHORT"
    [[ -n "${OUTPUT_DIR:-}"    ]] || die "Config is missing OUTPUT_DIR (or PROJECT_DIR)"
}

# 软加载：找不到/不合法时返回非 0 而不退出（供 doctor 等使用）
try_load_config() {
    local cfg=""
    if [[ -n "${1:-}" ]]; then cfg="$1"
    elif [[ -n "${CHROMAFLOW_CONFIG:-}" ]]; then cfg="${CHROMAFLOW_CONFIG}"
    elif [[ -f "./chromaflow.config" ]]; then cfg="$(pwd)/chromaflow.config"; fi
    [[ -n "${cfg}" && -f "${cfg}" ]] || return 1
    CHROMAFLOW_CONFIG="$(cd "$(dirname "${cfg}")" && pwd)/$(basename "${cfg}")"
    set +u
    # shellcheck disable=SC1090
    source "${CHROMAFLOW_CONFIG}"; derive_paths
    # shellcheck disable=SC1090
    source "${CHROMAFLOW_CONFIG}"; derive_paths
    set -u
    return 0
}

# 基础目录初始化（始终需要）
init_dirs() {
    mkdir -p "${REFERENCE_DIR}" "${LOG_DIR}" "${SAMPLESHEET_DIR}" \
             "${RAW_DATA_DIR}" "${NXF_SINGULARITY_CACHEDIR}" "${WORK_DIR}"
}

# 按需初始化某 pipeline 的输出目录
init_pipeline_dirs() {
    local pipeline="${1:?Usage: init_pipeline_dirs <chipseq|atacseq|dapseq>}"
    mkdir -p "${RESULTS_DIR}/downstream/${pipeline}"
}

# ─────────────── nf-core 输出路径帮助函数 ───────────────
get_nfcore_peak_dir() {
    local pipeline="${1}" peak_type="${2:-narrow_peak}"
    echo "${RESULTS_DIR}/${pipeline}/bwa/merged_library/macs3/${peak_type}"
}

get_consensus_peaks() {
    local pipeline="${1}" peak_dir found
    peak_dir=$(get_nfcore_peak_dir "${pipeline}")
    found=$(find "${peak_dir}/consensus" -name "*.consensus_peaks.bed" 2>/dev/null | head -1)
    if [[ -n "${found}" ]]; then echo "${found}"
    else echo "${peak_dir}/consensus/consensus_peaks.bed"; fi
}

# 给 consensus peaks 回填 signalValue（与各重复 narrowPeak 求交取 max），产出带富集值的 bed。
#
# 为什么必须回填：nf-core 的 consensus_peaks.bed 是这样生成的
#     awk 'FNR>1 {print $1,$2,$3,$4,"0","+"}' consensus_peaks.boolean.txt
# 第 5 列是**硬编码的字符串 "0"**，本来就不含任何富集信息。不回填的话，
# 05 的 peak_score、06 的排序、07 的 signalValue 标签全是 0。
#
# 结果缓存到 idr 目录，避免每次重算。输出: chr start end name signalValue strand
# 注意：本函数在 $(...) 中调用，所有提示必须走 stderr。
get_consensus_with_signal() {
    local pipeline="${1}"
    local cons; cons=$(get_consensus_peaks "${pipeline}")
    local idr_dir="${RESULTS_DIR}/downstream/${pipeline}/idr"
    local out="${idr_dir}/consensus_with_signal.bed"
    mkdir -p "${idr_dir}"
    [[ -s "${cons}" ]] || { echo "${cons}"; return 0; }
    if [[ -s "${out}" && "${out}" -nt "${cons}" ]]; then echo "${out}"; return 0; fi

    # 广搜各重复 narrowPeak（排除 consensus 子目录），兼容不同目录层级
    mapfile -t NP < <(find "${RESULTS_DIR}/${pipeline}" -name "*_peaks.narrowPeak" 2>/dev/null | grep -v "/consensus/" | sort)
    if [[ ${#NP[@]} -eq 0 ]]; then
        log_warn "Replicate narrowPeak files not found; cannot populate consensus signalValue" >&2
        echo "         → Downstream peak_score will be zero (nf-core consensus BED column 5 is zero)" >&2
        echo "${cons}"; return 0
    fi

    local allnp="${idr_dir}/.tmp_allnp_$$.bed"; cat "${NP[@]}" > "${allnp}"
    # 用一条数据行(第2、3列为整数)确定 consensus 列宽，避开表头
    local aw; aw=$(awk '$2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {print NF; exit}' "${cons}")
    if [[ -z "${aw}" ]]; then
        rm -f "${allnp}"
        log_warn "Unrecognized consensus file format: ${cons}" >&2
        echo "${cons}"; return 0
    fi

    local ok=1
    if command -v bedtools &>/dev/null; then
        if bedtools intersect -a "${cons}" -b "${allnp}" -wa -wb 2>/dev/null \
            | awk -v aw="${aw}" 'BEGIN{OFS="\t"} $2 ~ /^[0-9]+$/ {
                  sig=$(aw+7)+0; key=$1"\t"$2"\t"$3"\t"$4;
                  if(sig>m[key]) m[key]=sig }
                  END{ for(k in m){ split(k,a,"\t"); nm=(a[4]==""?"peak":a[4]); print a[1],a[2],a[3],nm,m[k],"+" } }' \
            | sort -k1,1 -k2,2n > "${out}" 2>/dev/null && [[ -s "${out}" ]]; then
            ok=0
        fi
    else
        # 无 bedtools 时用 python3 做区间重叠（05 的 ChromaFlow_R 环境通常没有 bedtools）
        log_info "bedtools not found; using Python to populate consensus signalValue" >&2
        if python3 "${CHROMAFLOW_HOME}/lib/steps/support/consensus_signal.py" \
               "${cons}" "${allnp}" "${out}" >&2 && [[ -s "${out}" ]]; then
            ok=0
        fi
    fi
    rm -f "${allnp}"

    if [[ ${ok} -eq 0 ]]; then
        echo "${out}"
    else
        rm -f "${out}"          # 别把半成品留在缓存位置，否则下次会被当成有效结果
        log_warn "Could not populate consensus signalValue; using original consensus (peak_score is zero)" >&2
        echo "         Requires bedtools or python3; alternatively use --peakset idr" >&2
        echo "${cons}"
    fi
}

# 解析 --peakset：决定下游用哪套 peak + 输出目录后缀，实现"跳过/使用 IDR"的开关。
# 设置全局: PEAKSET_FILE, PEAKSET_SUFFIX
#   idr | idr:0.01   → IDR<0.01 (recommended)，后缀 ""
#   idr:0.05|0.1     → 对应阈值文件，后缀 "_idr<阈值>"
#   consensus        → nf-core consensus(回填 signalValue)，后缀 "_consensus"
#   <文件路径>       → 直接用该文件，后缀 "_custom"
resolve_peakset() {
    local pipeline="${1}" spec="${2:-idr}"
    local idr_dir="${RESULTS_DIR}/downstream/${pipeline}/idr"
    case "${spec}" in
        idr|idr:*)
            if [[ -f "${idr_dir}/idr.status" && "$(cat "${idr_dir}/idr.status")" != complete ]]; then
                die 'The latest IDR run did not complete; rerun IDR before using its peak sets.'
            fi
            local thr="${spec#idr}"; thr="${thr#:}"; thr="${thr:-0.01}"
            PEAKSET_FILE="${idr_dir}/idr_peaks_threshold_${thr}.bed"
            [[ -s "${PEAKSET_FILE}" ]] || die "Requested IDR peak set is missing or empty: ${PEAKSET_FILE}. Run idr or explicitly choose --peakset consensus."
            if [[ "${thr}" == "0.01" ]]; then PEAKSET_SUFFIX=""; else PEAKSET_SUFFIX="_idr${thr}"; fi
            ;;
        consensus)
            PEAKSET_FILE="$(get_consensus_with_signal "${pipeline}")"
            PEAKSET_SUFFIX="_consensus"
            ;;
        *)
            PEAKSET_FILE="${spec}"; PEAKSET_SUFFIX="_custom"
            [[ -s "${PEAKSET_FILE}" ]] || die "Peak-set file does not exist: ${PEAKSET_FILE}"
            ;;
    esac
    export PEAKSET_FILE PEAKSET_SUFFIX
}

# 打开日志（tee 到 OUTPUT_DIR/logs）
open_log() {
    [[ "${CHROMAFLOW_LOG_MANAGED:-0}" != 1 ]] || return 0
    local step="$1" pipeline="${2:-}"
    local name="${step}${pipeline:+_${pipeline}}_${TIMESTAMP}_${BASHPID}.log"
    mkdir -p "${LOG_DIR}"
    exec > >(tee -a "${LOG_DIR}/${name}") 2>&1
    echo "[LOG] ${LOG_DIR}/${name}"
}
