#!/usr/bin/env bash
# =============================================================================
# md_watch.sh - progress of the docked-pose MD replicates, with an ETA.
#
# Reports for every compound x replicate under 05_MD_v2:
#   stage (not started / NVT / NPT / production / COMPLETE)
#   percent done and simulated time for the current stage
#   whether a gmx mdrun process is actually alive in that directory
#   measured ns/day, and an estimated finish time for the whole campaign
#
# Usage
#   ./md_watch.sh                 one-shot report
#   ./md_watch.sh -w              refresh every 60 s (Ctrl-C to stop)
#   ./md_watch.sh -w 30           refresh every 30 s
#   ./md_watch.sh -s              one short line, for a prompt or a cron mail
#   ./md_watch.sh -r /other/path  point at a different 05_MD_v2
#
# Safe to run while the simulations are going. It only reads files.
# =============================================================================
set -u

MD="${MD_ROOT:-$HOME/Research/05_MD_v2}"
WATCH=0; INTERVAL=60; SHORT=0
while [ $# -gt 0 ]; do
    case "$1" in
        -w|--watch) WATCH=1; [[ "${2:-}" =~ ^[0-9]+$ ]] && { INTERVAL=$2; shift; } ;;
        -s|--short) SHORT=1 ;;
        -r|--root)  MD="$2"; shift ;;
        -h|--help)  sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
    shift
done

COMPOUNDS=(LOA22-B1 LOA22-B2 LOA22-B3)
REPS=(1 2 3)
NVT_STEPS=250000        # 500 ps at dt 0.002
NPT_STEPS=1000000       # 2 ns
PROD_STEPS=50000000     # 100 ns
PROD_NS=100
FALLBACK_NSDAY=120.307  # what rep1 actually achieved on 117,185 atoms

# ---------------------------------------------------------------- helpers
newest_log() {   # $1 = dir, $2 = deffnm. Handles md_prod.log and md_prod.partNNNN.log
    ls -t "$1/$2".log "$1/$2".part[0-9]*.log 2>/dev/null | head -1
}

have_gro() {     # $1 = dir, $2 = deffnm. True for md_prod.gro OR md_prod.partNNNN.gro
    compgen -G "$1/$2.gro" >/dev/null 2>&1 && return 0
    compgen -G "$1/$2.part[0-9]*.gro" >/dev/null 2>&1
}

live_pid() {     # $1 = dir, $2 = deffnm -> pid of a gmx mdrun running there
    local dir want pid cwd cmd
    dir=$(readlink -f "$1" 2>/dev/null) || return
    want="$2"
    while read -r pid; do
        [ -z "$pid" ] && continue
        cwd=$(readlink -f "/proc/$pid/cwd" 2>/dev/null) || continue
        [ "$cwd" = "$dir" ] || continue
        cmd=$(tr '\0' ' ' < "/proc/$pid/cmdline" 2>/dev/null)
        case "$cmd" in *"-deffnm $want"*) echo "$pid"; return ;; esac
    done < <(pgrep -f "gmx[^ ]* mdrun" 2>/dev/null)
}

stage_progress() {   # $1 = logfile, $2 = target steps -> "pct ns"
    local f=$1 target=$2 step
    [ -n "$f" ] && [ -f "$f" ] || { echo "0.0 0.00"; return; }
    step=$(grep -A1 -E "^ *Step +Time" "$f" 2>/dev/null | tail -1 | awk '{print $1}')
    [[ "$step" =~ ^[0-9]+$ ]] || { echo "0.0 0.00"; return; }
    awk -v s="$step" -v t="$target" 'BEGIN{printf "%.1f %.2f\n", 100*s/t, s*0.002/1000}'
}

perf_of() {      # measured ns/day from a finished stage
    [ -n "${1:-}" ] && [ -f "$1" ] || return
    grep -E "^Performance:" "$1" 2>/dev/null | tail -1 | awk '{print $2}'
}

# ---------------------------------------------------------------- report
report() {
    local done_ns=0 target_ns=$(( ${#COMPOUNDS[@]} * ${#REPS[@]} * PROD_NS ))
    local live=0 nsday="" out=""

    for C in "${COMPOUNDS[@]}"; do
        out+=$'\n'"── $C ───────────────────────────────────────────"$'\n'
        for R in "${REPS[@]}"; do
            local D="$MD/$C/rep$R" label stage detail lg pid pct ns p
            label=$(printf "rep%-2s" "$R")

            if [ ! -d "$D" ]; then
                stage="not started"; detail="-"
            elif have_gro "$D" md_prod; then
                lg=$(newest_log "$D" md_prod); p=$(perf_of "$lg")
                [ -n "$p" ] && nsday="$p"
                stage="production COMPLETE"; detail="100.0%   ${p:-?} ns/day"
                done_ns=$(awk -v a="$done_ns" -v n="$PROD_NS" 'BEGIN{printf "%.1f", a+n}')
            elif [ -f "$D/md_prod.tpr" ]; then
                lg=$(newest_log "$D" md_prod); pid=$(live_pid "$D" md_prod)
                read -r pct ns <<< "$(stage_progress "$lg" "$PROD_STEPS")"
                if [ -n "$pid" ]; then stage="production RUNNING"; live=1
                else stage="production STALLED"; fi
                detail="${ns} / ${PROD_NS} ns   (${pct}%)"
                done_ns=$(awk -v a="$done_ns" -v n="$ns" 'BEGIN{printf "%.1f", a+n}')
            elif have_gro "$D" npt; then
                stage="NPT done, prod not built"; detail="waiting on grompp"
            elif [ -f "$D/npt.tpr" ]; then
                lg=$(newest_log "$D" npt); pid=$(live_pid "$D" npt)
                read -r pct ns <<< "$(stage_progress "$lg" "$NPT_STEPS")"
                [ -n "$pid" ] && { stage="NPT RUNNING"; live=1; } || stage="NPT STALLED"
                detail="${ns} / 2.0 ns   (${pct}%)"
            elif have_gro "$D" nvt; then
                stage="NVT done, NPT not built"; detail="waiting on grompp"
            elif [ -f "$D/nvt.tpr" ]; then
                lg=$(newest_log "$D" nvt); pid=$(live_pid "$D" nvt)
                read -r pct ns <<< "$(stage_progress "$lg" "$NVT_STEPS")"
                [ -n "$pid" ] && { stage="NVT RUNNING"; live=1; } || stage="NVT STALLED"
                detail="${ns} / 0.5 ns   (${pct}%)"
            else
                stage="build present, not begun"; detail="-"
            fi
            out+=$(printf "    %-6s %-28s %s" "$label" "$stage" "$detail")$'\n'
        done
    done

    [ -z "$nsday" ] && nsday=$FALLBACK_NSDAY
    local pct_all remain eta
    pct_all=$(awk -v d="$done_ns" -v t="$target_ns" 'BEGIN{printf "%.1f", 100*d/t}')
    remain=$(awk -v d="$done_ns" -v t="$target_ns" 'BEGIN{printf "%.1f", t-d}')
    local hours; hours=$(awk -v r="$remain" -v p="$nsday" 'BEGIN{printf "%.1f", 24*r/p}')
    local mins; mins=$(awk -v h="$hours" 'BEGIN{printf "%d", h*60}')
    eta=$(date -d "+${mins} minutes" '+%a %d %b %H:%M' 2>/dev/null \
          || date -v+"${mins}"M '+%a %d %b %H:%M' 2>/dev/null || echo "?")

    if [ "$SHORT" -eq 1 ]; then
        printf "MD %s/%s ns (%s%%)  %s ns/day  ETA %s\n" \
               "$done_ns" "$target_ns" "$pct_all" "$nsday" "$eta"
        return
    fi

    echo "════════════════════════════════════════════════════════"
    echo "  DOCKED-POSE MD PROGRESS      $(date '+%Y-%m-%d %H:%M:%S')"
    echo "  root: $MD"
    echo "════════════════════════════════════════════════════════"
    printf '%s' "$out"
    echo
    echo "────────────────────────────────────────────────────────"
    printf "  TOTAL    %s / %s ns  (%s%%)\n" "$done_ns" "$target_ns" "$pct_all"
    printf "  RATE     %s ns/day\n" "$nsday"
    printf "  REMAIN   %s ns, about %s hours\n" "$remain" "$hours"
    printf "  ETA      %s\n" "$eta"
    echo "────────────────────────────────────────────────────────"

    echo "  load  : $(cut -d' ' -f1-3 /proc/loadavg)   cores: $(nproc)"
    echo "  memory: $(free -h 2>/dev/null | awk '/^Mem:/{print $3" used of "$2}')"
    if command -v nvidia-smi >/dev/null 2>&1; then
        nvidia-smi --query-gpu=utilization.gpu,memory.used,memory.total,temperature.gpu \
            --format=csv,noheader 2>/dev/null \
            | awk -F', ' '{printf "  GPU   : %s util, %s of %s, %s\n",$1,$2,$3,$4}'
    fi

    if [ "$live" -eq 0 ] && [ "$pct_all" != "100.0" ]; then
        echo
        echo "  !! No gmx mdrun process is running, and the campaign is not finished."
        echo "     Check:  ps aux | grep mdrun      tmux attach -t loa22"
        echo "     Resume: cd <the stalled replicate> && \\"
        echo "             gmx mdrun -deffnm md_prod -cpi md_prod.cpt -append \$GMXRUN"
    fi
    echo "════════════════════════════════════════════════════════"
}

if [ "$WATCH" -eq 1 ]; then
    trap 'echo; echo "stopped watching"; exit 0' INT
    while true; do clear; report; sleep "$INTERVAL"; done
else
    report
fi
