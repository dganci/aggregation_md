# shellcheck shell=bash
# Assegnazione delle GPU alle entry concorrenti.
#
# L'invariante che conta: due entry vive non devono MAI vedere lo stesso
# dispositivo. Il caso che rompe l'aritmetica "indice modulo numero di GPU" e'
# che le entry finiscono in ordine sparso - 5x1-108 e 5x1-263 differiscono di un
# ordine di grandezza in dimensione - quindi la quinta puo' partire mentre la
# prima e' ancora viva, e prenderebbe il suo stesso dispositivo.
source "$ROOT/tools/lib/gpu.sh"

check "gpu_devices espande un conteggio" "$(GPUS=4 gpu_devices | tr '\n' ' ')" "0 1 2 3 "
check "gpu_devices espande una lista"    "$(GPUS=0,2,5 gpu_devices | tr '\n' ' ')" "0 2 5 "
check "senza GPUS non espande nulla"     "$(GPUS= gpu_devices | tr '\n' ' ')" ""

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

check "ogni lavoro ha ricevuto un dispositivo" \
      "$(grep -c start "$TMP/gpu_used.txt" | tr -d ' ')" "4"
check "sono stati usati solo i dispositivi dichiarati" \
      "$(awk '{print $1}' "$TMP/gpu_used.txt" | sort -u | tr '\n' ' ')" "0 1 "
check "nessuna sovrapposizione sullo stesso dispositivo" \
      "$(awk '{if($3=="start"){c[$1]++; if(c[$1]>1) bad=1} else c[$1]--}
              END{print bad ? "COLLISIONE" : "ok"}' "$TMP/gpu_used.txt")" "ok"
