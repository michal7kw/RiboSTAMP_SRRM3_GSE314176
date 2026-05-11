#!/bin/bash
# Detach Surface A from the parent shell so it survives wsl session exit.
LOG=/tmp/surface_a.log
echo "" > "${LOG}"

setsid bash /mnt/d/Github/SRF/Elisa/RiboSTAMP_SRRM3_GSE314176/pc_pipeline/sr_junction_psi/_run_surface_a.sh \
    >> "${LOG}" 2>&1 < /dev/null &

PID=$!
echo "Detached PID: ${PID}"
echo "Logging to: ${LOG}"

# Brief sleep to let it start
sleep 3
if ps -p "${PID}" > /dev/null 2>&1; then
    echo "Process is running."
else
    echo "WARN: Process ${PID} died immediately. Tail of log:"
    tail -20 "${LOG}"
fi
