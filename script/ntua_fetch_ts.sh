#! /bin/bash

HYPATIA_IP="147.102.106.110"
PARASKEVAS_IP="147.102.106.113"
TS_PATH=/home/bpe/data/time-series

if test "$#" -lt 1 ; then
  echo "[ERROR] Need to provide one or more station names"
  exit 1
fi

TM_STAMP=$(date +"%Y%m%d%H%M")
for station in "$@" ; do
  ts_hyp=${station}.cts.hypatia.${TM_STAMP}
  ts_par=${station}.cts.paraskevas.${TM_STAMP}
  scp bpe@${HYPATIA_IP}:${TS_PATH}/${station}/${station}.cts $ts_hyp
  scp bpe@${PARASKEVAS_IP}:${TS_PATH}/${station}/${station}.cts $ts_par
done
