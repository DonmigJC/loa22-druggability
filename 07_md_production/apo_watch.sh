#!/usr/bin/env bash
# =============================================================================
# apo_watch.sh - progress of the three ligand-free (apo) v2 replicates.
#
# The apo set is the last thing standing between the manuscript and the
# journal's three-replica condition, and nobody is driving it by hand: cron
# fires run_queue_resumable.sh every 15 minutes, that script takes an flock and
# refuses to start while any gmx process is live, and when the running
# replicate finishes it calls build_apo_v2.sh, which loops rep1..rep3 at seeds
# 42, 1337, 24601 and skips any replicate that already has md_prod.gro.
#
# So the two things worth watching are the replicate in front of you and the
# health of that chain. This reports both.
#
# Usage
#   ./apo_watch.sh                one-shot report
#   ./apo_watch.sh -w             refresh every 60 s (Ctrl-C to stop)
#   ./apo_watch.sh -w 300         refresh every 5 min
#   ./apo_watch.sh -s             one short line, for a prompt or cron mail
#   ./apo_watch.sh -r /other/path point at a different 06_MD_apo_v2
#
# Read-only. Safe to run as often as you like while the simulations go.
#
# The rate is self-calibrating: each run drops a sample in .apo_watch_state
# and the next run measures ns/day from the gap. The first run after a long
# pause falls back to 113.3 ns/day, which is what apo rep1 actually achieved.
# =============================================================================
set -u

APO="${APO_ROOT:-$HOME/Research/06_MD_apo_v2}"
MD="${MD_ROOT:-$HOME/Research/05_MD_v2}"
WATCH=0; INTERVAL=60; SHORT=0

while [ $# -gt 0 ]; do
    case "$1" in
        -w|--watch) WATCH=1; [[ "${2:-}" =~ ^[0-9]+$ ]] && { INTERVAL=$2; shift; } ;;
        -s|--short) SHORT=1 ;;
        -r|--root)  APO="$2"; shift ;;
        -h|--help)  sed -n '2,29p' "$0"; exit 0 ;;
        *) echo "unknown option: $1" >&2; exit 1 ;;
    esac
    shift
done

REPS=(1 2 3)
SEEDS=(42 1337 24601)
NVT_STEPS=250000        # 500 ps at dt 0.002
NPT_STEPS=1000000       # 2 ns
PROD_STEPS=50000000     # 100 ns
PROD_NS=100
EQUIL_NS=2.5            # 0.5 ns NVT + 2 ns NPT, per replicate
FALLBACK_NSDAY=113.3    # measured on apo rep1, 10-11 Nov 2026
STATE="$APO/.apo_watch_state"

# ---------------------------------------------------------------- helpers
newest_log() {   # dir, deffnm -> most recent log, including .partNNNN.log
    ls -t "$1/$2".log "$1/$2".part[0-9]*.log 2>/dev/null | head -1
}

have_gro() {     # dir, deffnm -> true if the stage finished
    compgen -G "$1/$2.gro" >/dev/null 2>&1 && return 0
    compgen -G "$1/$2.part[0-9]*.gro" >/dev/null 2>&1
}

live_pid() {     # dir, deffnm -> pid of a gmx mdrun running there
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

stage_progress() {   # logfile, target steps -> "pct ns"
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

hms() {          # hours -> "2 d 04 h 11 m"
    awk -v h="$1" 'BEGIN{
        if (h < 0) h = 0
        d = int(h/24); r = h - 24*d; hh = int(r); mm = int((r-hh)*60+0.5)
        if (mm == 60) { mm = 0; hh++ }
        if (d > 0) printf "%d d %02d h %02d m", d, hh, mm
        else       printf "%d h %02d m", hh, mm
    }'
}

# ---------------------------------------------------------------- the chain
queue_health() {
    local qpid lockholder cronline gmxpid gmxcwd
    qpid=$(pgrep -f "run_queue_resumable.sh" 2>/dev/null | tr '\n' ' ')
    cronline=$(crontab -l 2>/dev/null | grep -c "run_queue_resumable.sh")
    gmxpid=$(pgrep -f "gmx[^ ]* mdrun" 2>/dev/null | head -1)
    gmxcwd=""
    [ -n "$gmxpid" ] && gmxcwd=$(readlink -f "/proc/$gmxpid/cwd" 2>/dev/null)

    if [ -n "$qpid" ]; then
        echo "  queue    : ALIVE (pid${qpid% })"
    else
        echo "  queue    : not running"
    fi
    if [ "$cronline" -ge 1 ]; then
        echo "  cron     : $cronline entry/entries present, fires every 15 min"
    else
        echo "  cron     : !! NO CRON ENTRY - nothing will restart the queue"
    fi
    if [ -n "$gmxpid" ]; then
        echo "  mdrun    : pid $gmxpid in ${gmxcwd/#$HOME/~}"
    else
        echo "  mdrun    : none live"
    fi

    # will the queue actually reach the apo stage?
    local blocked=0 c r
    for c in LOA22-B1 LOA22-B2 LOA22-B3; do
        for r in rep1 rep2 rep3; do
            have_gro "$MD/$c/$r" md_prod || blocked=$((blocked+1))
        done
    done
    if [ "$blocked" -eq 0 ]; then
        echo "  docked   : all 9 complete, so the queue falls through to the apo stage"
    else
        echo "  docked   : $blocked of 9 incomplete - the queue will do those first"
    fi
    if [ -z "$qpid" ] && [ -z "$gmxpid" ]; then
        echo
        echo "  Nothing is running. To restart the chain by hand:"
        echo "    setsid nohup $MD/run_queue_resumable.sh >> $MD/queue.log 2>&1 < /dev/null &"
    fi
}

# ---------------------------------------------------------------- report
report() {
    local done_ns=0 target_ns=$(( ${#REPS[@]} * PROD_NS ))
    local live=0 nsday="" out="" running_rep="" remaining_equil=0 i=0

    for R in "${REPS[@]}"; do
        local D="$APO/rep$R" stage detail lg pid pct ns p seed
        seed=${SEEDS[$i]}; i=$((i+1))

        if [ ! -d "$D" ]; then
            stage="not started"; detail="seed $seed, builds when the queue reaches it"
            remaining_equil=$(awk -v a="$remaining_equil" -v e="$EQUIL_NS" 'BEGIN{print a+e}')
        elif have_gro "$D" md_prod; then
            lg=$(newest_log "$D" md_prod); p=$(perf_of "$lg")
            [ -n "$p" ] && nsday="$p"
            stage="production COMPLETE"; detail="100.0%   ${p:-?} ns/day   seed $seed"
            done_ns=$(awk -v a="$done_ns" -v n="$PROD_NS" 'BEGIN{printf "%.2f", a+n}')
        elif [ -f "$D/md_prod.tpr" ]; then
            lg=$(newest_log "$D" md_prod); pid=$(live_pid "$D" md_prod)
            read -r pct ns <<< "$(stage_progress "$lg" "$PROD_STEPS")"
            if [ -n "$pid" ]; then
                stage="production RUNNING"; live=1; running_rep="rep$R"
            else
                stage="production STALLED"
            fi
            detail="${ns} / ${PROD_NS} ns   (${pct}%)   seed $seed"
            done_ns=$(awk -v a="$done_ns" -v n="$ns" 'BEGIN{printf "%.2f", a+n}')
        elif have_gro "$D" npt; then
            stage="NPT done, prod not built"; detail="waiting on grompp"
        elif [ -f "$D/npt.tpr" ]; then
            lg=$(newest_log "$D" npt); pid=$(live_pid "$D" npt)
            read -r pct ns <<< "$(stage_progress "$lg" "$NPT_STEPS")"
            [ -n "$pid" ] && { stage="NPT RUNNING"; live=1; running_rep="rep$R"; } \
                          || stage="NPT STALLED"
            detail="${ns} / 2.0 ns   (${pct}%)   equilibration"
        elif have_gro "$D" nvt; then
            stage="NVT done, NPT not built"; detail="waiting on grompp"
        elif [ -f "$D/nvt.tpr" ]; then
            lg=$(newest_log "$D" nvt); pid=$(live_pid "$D" nvt)
            read -r pct ns <<< "$(stage_progress "$lg" "$NVT_STEPS")"
            [ -n "$pid" ] && { stage="NVT RUNNING"; live=1; running_rep="rep$R"; } \
                          || stage="NVT STALLED"
            detail="${ns} / 0.5 ns   (${pct}%)   equilibration"
        else
            stage="build present, not begun"; detail="seed $seed"
        fi
        out+=$(printf "    %-6s %-28s %s\n" "rep$R" "$stage" "$detail")$'\n'
    done

    # ---- self-calibrating rate from the previous invocation ---------------
    local now prev_t prev_ns measured="" dt dn
    now=$(date +%s)
    if [ -f "$STATE" ]; then
        read -r prev_t prev_ns < "$STATE" 2>/dev/null || true
        if [[ "${prev_t:-}" =~ ^[0-9]+$ ]]; then
            dt=$((now - prev_t))
            dn=$(awk -v a="$done_ns" -v b="${prev_ns:-0}" 'BEGIN{printf "%.4f", a-b}')
            # A sample spanning a replicate changeover reads far too low: the
            # production counter stalls for the ~30 min of em, NVT and NPT that
            # the next replicate spends equilibrating. Require a 10 min window
            # and discard any rate outside a plausible band for this hardware,
            # falling back to the ns/day a finished replicate actually achieved.
            if [ "$dt" -ge 600 ] && awk -v d="$dn" 'BEGIN{exit !(d>0.01)}'; then
                measured=$(awk -v d="$dn" -v t="$dt" 'BEGIN{printf "%.1f", d*86400/t}')
                awk -v m="$measured" 'BEGIN{exit !(m>=70 && m<=250)}' || measured=""
            fi
        fi
    fi
    printf "%s %s\n" "$now" "$done_ns" > "$STATE" 2>/dev/null || true

    local rate="${measured:-${nsday:-$FALLBACK_NSDAY}}"
    local rate_src="measured since last check"
    [ -z "$measured" ] && rate_src="from a completed replicate"
    [ -z "$measured" ] && [ -z "$nsday" ] && rate_src="fallback, apo rep1"

    local remain_ns hours mins eta pct_total
    remain_ns=$(awk -v t="$target_ns" -v d="$done_ns" 'BEGIN{printf "%.1f", t-d}')
    hours=$(awk -v r="$remain_ns" -v e="$remaining_equil" -v n="$rate" \
                'BEGIN{ if (n<=0) n=1; printf "%.3f", (r+e)*24/n }')
    # date(1) will not take fractional hours in a relative item, so pass minutes
    mins=$(awk -v h="$hours" 'BEGIN{printf "%d", h*60+0.5}')
    eta=$(date -d "+${mins} minutes" "+%a %d %b %H:%M" 2>/dev/null \
          || date -v"+${mins}M" "+%a %d %b %H:%M" 2>/dev/null || echo "?")
    pct_total=$(awk -v d="$done_ns" -v t="$target_ns" 'BEGIN{printf "%.1f", 100*d/t}')

    if [ "$SHORT" -eq 1 ]; then
        printf "apo %s/%s ns (%s%%) %s %s ns/day ETA %s\n" \
            "$done_ns" "$target_ns" "$pct_total" \
            "${running_rep:-idle}" "$rate" "$eta"
        return
    fi

    echo "════════════════════════════════════════════════════════════════"
    printf "  APO v2 REPLICATES        %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
    printf "  root: %s\n" "${APO/#$HOME/~}"
    echo "════════════════════════════════════════════════════════════════"
    echo
    printf "%s" "$out"
    echo
    echo "────────────────────────────────────────────────────────────────"
    printf "  TOTAL    %s / %s ns  (%s%%)\n" "$done_ns" "$target_ns" "$pct_total"
    printf "  RATE     %s ns/day   (%s)\n" "$rate" "$rate_src"
    printf "  REMAIN   %s ns production" "$remain_ns"
    awk -v e="$remaining_equil" 'BEGIN{ if (e>0) printf " + %.1f ns equilibration", e }'
    echo
    printf "           about %s\n" "$(hms "$hours")"
    printf "  ETA      %s\n" "$eta"
    echo "────────────────────────────────────────────────────────────────"
    echo "  THE CHAIN"
    queue_health
    echo "────────────────────────────────────────────────────────────────"
    if command -v nvidia-smi >/dev/null 2>&1; then
        printf "  GPU   : %s\n" "$(nvidia-smi --query-gpu=utilization.gpu,memory.used,temperature.gpu \
                                   --format=csv,noheader 2>/dev/null | head -1)"
    fi
    printf "  load  : %s   cores: %s\n" \
        "$(cut -d' ' -f1-3 /proc/loadavg)" "$(nproc)"
    printf "  disk  : %s free on %s\n" \
        "$(df -h "$APO" 2>/dev/null | awk 'NR==2{print $4}')" \
        "$(df -h "$APO" 2>/dev/null | awk 'NR==2{print $6}')"

    # ---- the only thing that matters for the manuscript ------------------
    echo "────────────────────────────────────────────────────────────────"
    local ncomplete=0
    for R in "${REPS[@]}"; do have_gro "$APO/rep$R" md_prod && ncomplete=$((ncomplete+1)); done
    if [ "$ncomplete" -eq 3 ]; then
        echo "  ALL THREE APO REPLICATES COMPLETE."
        echo "  Next, and it is the last data step before writing:"
        echo "    bash ~/Desktop/Loa22_Research/02_scripts/analyse_conformational_stability_v2.sh apo"
        echo "    python3 .../summarize_conformational_stability.py ~/Research"
        echo "  That fills the apo row of Table S12, the section 2.6 baseline,"
        echo "  the section 3.4 caveat and the cover-letter Overview sentence."
    else
        printf "  %d of 3 apo replicates complete. Nothing to do but wait.\n" "$ncomplete"
        echo "  Do not start mdrun by hand: the queue holds an flock and will"
        echo "  refuse to run while a gmx process is live, which would stall it."
    fi
    echo "════════════════════════════════════════════════════════════════"
}

if [ "$WATCH" -eq 1 ]; then
    trap 'echo; exit 0' INT
    while :; do clear; report; sleep "$INTERVAL"; done
else
    report
fi
