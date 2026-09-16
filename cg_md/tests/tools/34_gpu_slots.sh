# shellcheck shell=bash
# GPU assignment across the concurrent entry pool.
#
# The invariant that matters: two live entries must NEVER see the same device.
# What breaks "entry index modulo number of GPUs" is that entries finish out of
# order - 5x1-108 and 5x116-347 differ by a factor of ten in size, and the
# adaptive stop condition fires at different times - so the fifth entry can
# start while the first is still running, and would take its device.
source "$ROOT/tools/lib/gpu.sh"

check "gpu_devices expands a count" "$(GPUS=4 gpu_devices | tr '\n' ' ')" "0 1 2 3 "
check "gpu_devices expands a list"  "$(GPUS=0,2,5 gpu_devices | tr '\n' ' ')" "0 2 5 "
check "without GPUS it expands nothing" "$(GPUS= gpu_devices | tr '\n' ' ')" ""

export GPU_LOCK_DIR="$TMP/gpulocks"
rm -rf "$GPU_LOCK_DIR"
: > "$TMP/gpu_used.txt"
(
  export GPUS=2
  for n in 1 2 3 4; do
    ( gpu_run 1 "j$n" bash -c \
        "echo \"\$CUDA_VISIBLE_DEVICES \$\$ start\" >> '$TMP/gpu_used.txt'
         sleep 1
         echo \"\$CUDA_VISIBLE_DEVICES \$\$ end\" >> '$TMP/gpu_used.txt'" \
    ) >/dev/null 2>&1 &
  done
  wait
)

check "every job was handed a device" \
      "$(grep -c start "$TMP/gpu_used.txt" | tr -d ' ')" "4"
check "only the declared devices were used" \
      "$(awk '{print $1}' "$TMP/gpu_used.txt" | sort -u | tr '\n' ' ')" "0 1 "
check "no two jobs overlapped on one device" \
      "$(awk '{if($3=="start"){c[$1]++; if(c[$1]>1) bad=1} else c[$1]--}
              END{print bad ? "COLLISION" : "ok"}' "$TMP/gpu_used.txt")" "ok"
