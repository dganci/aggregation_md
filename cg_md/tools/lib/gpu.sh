# GPU slot allocation for the concurrent entry pool.
#
# Why locks and not "entry index modulo number of GPUs": the pool in
# run_batch.sh starts a new entry as soon as ANY running one finishes, and
# entries do not finish in the order they started - 5x1-108 and 5x1-263 differ
# by an order of magnitude in size, and the adaptive stop condition fires at
# different times. With index arithmetic, entry 4 would take device 0 while
# entry 0 is still on it, and two mdrun processes would share one A100 while
# another sat idle. A lock is the only thing that tracks what is ACTUALLY free.
#
# Usage:  GPUS=4  (or GPUS=0,1,2,3) in the environment.

GPU_LOCK_DIR="${GPU_LOCK_DIR:-${TMPDIR:-/tmp}/cg_md_gpu_locks.$$}"

# Expands GPUS into the array of device ids to hand out.
gpu_devices() {
    local spec="${GPUS:-}"
    [[ -z "$spec" ]] && return 1
    if [[ "$spec" =~ ^[0-9]+$ ]]; then
        local i
        for ((i = 0; i < spec; i++)); do printf '%s\n' "$i"; done
    else
        printf '%s\n' "${spec//,/$'\n'}"
    fi
}

# gpu_run <log-fd> <entry-name> <command...>
# Holds one device for the whole command and releases it on exit.
gpu_run() {
    local fd="$1" name="$2"; shift 2
    local -a devices
    mapfile -t devices < <(gpu_devices)
    if [[ "${#devices[@]}" -eq 0 ]]; then
        echo "GPUS='$GPUS' names no device" >&"$fd"
        return 2
    fi
    mkdir -p "$GPU_LOCK_DIR"

    local dev lockfd
    while :; do
        for dev in "${devices[@]}"; do
            exec {lockfd}>"$GPU_LOCK_DIR/gpu.$dev"
            if flock -n "$lockfd"; then
                echo "    $name: GPU $dev" >&"$fd"
                # The subshell inherits the locked descriptor and nothing else
                # holds it, so the lock lifetime is exactly the command's.
                ( export CUDA_VISIBLE_DEVICES="$dev"
                  # GROMACS numbers devices within what it can SEE, so after the
                  # filter above the only visible device is index 0. Passing the
                  # original id here would ask for a device that is not there.
                  export GMX_GPU_ID=0
                  "$@" ) >&"$fd" 2>&"$fd"
                local code=$?
                exec {lockfd}>&-
                return $code
            fi
            exec {lockfd}>&-
        done
        # Every device busy: the pool should prevent this, so it means JOBS
        # exceeds the number of GPUs. Wait rather than oversubscribe silently.
        sleep 5
    done
}
