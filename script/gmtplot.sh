#! /usr/bin/bash

echoerr() { echo "$@" 1>&2; }

##  Use YYYY-MM-DD format (for plotting)
USE_YMD_FORMAT=0
##  Do not plot outliers
REMOVE_OUTLIERS=0
##  current pid
pid=$(echo $$)
##  number of ticks on Y-axis
Y_TICKS_NR=6
## output file
ps=koko.ps
rm $ps 2>/dev/null

while test $# -gt 0
do
    key="$1"
    case $key in
        -i|--input)
        cts_file="$2"
        shift 2
        ;;
        -y|--ymd)
        USE_YMD_FORMAT=1
        shift
        ;;
        -n|--no-outliers)
        REMOVE_OUTLIERS=1
        shift
        ;;
        *) # unknown option
        echoerr "Fuck is this: \"$key\"? skipped ..."
        shift
        ;;
    esac
done

##  TODO check that cts file exists

## tmps file name
cts_np=$(basename $cts_file)
tmpf=".${cts_np}.${pid}"
cat ${cts_file} > ${tmpf}
cts_file=${tmpf}

## clear outlier if needed
if test "${REMOVE_OUTLIERS}" -eq 1
then
    sed '/o/d' ${cts_file} > ${cts_file}.tmp
    mv ${cts_file}.tmp ${cts_file}
fi

## transform to yyyy-mm-dd if needed, and get min/max values
if test "${USE_YMD_FORMAT}" -eq 1
then
    gmt set FORMAT_DATE_IN yyyy-mm-dd
    gmt set FORMAT_DATE_MAP o
    gmt set FORMAT_TIME_PRIMARY_MAP abbreviated
    gmt set FORMAT_DATE_OUT yyyy-o-dd
    if ! ./mjd2ymdhms.py -f ${cts_file} -s 1 -d ',' > ${cts_file}.tmp
    then
        echoerr "[ERROR] Failed to transform dates in input file."
        exit 1
    fi
    mv ${cts_file}.tmp ${cts_file}
    wesn=( $(gmt info -fT -I0.01 -C @"${cts_file}" --FORMAT_DATE_IN=yyyy-mm-dd -h --IO_N_HEADER_RECS=1) )
else
    wesn=( $(gmt info -I0.01 -C @"${cts_file}" -h --IO_N_HEADER_RECS=1) )
fi

##  Limits for x,y,z / n,e,u plots
xmin="${wesn[2]}"
xmax="${wesn[3]}"
Rx="-R${wesn[0]}/${wesn[1]}/${xmin}/${xmax}"
ymin="${wesn[6]}"
ymax="${wesn[7]}"
Ry="-R${wesn[0]}/${wesn[1]}/${ymin}/${ymax}"
zmin="${wesn[10]}"
zmax="${wesn[11]}"
Rz="-R${wesn[0]}/${wesn[1]}/${zmin}/${zmax}"

# x/y limits array
R=("${Rx}" "${Ry}" "${Rz}")
# file column array
Y=("2" "5" "8")
# y annotation every
tcks1=$(./annotevery.py ${xmin} ${xmax} ${Y_TICKS_NR})
tcks2=$(./annotevery.py ${ymin} ${ymax} ${Y_TICKS_NR})
tcks3=$(./annotevery.py ${zmin} ${zmax} ${Y_TICKS_NR})
YAI=("${tcks1}" "${tcks2}" "${tcks3}")
# x/y annotation ticks on axis
XAT=("WseN" "Wsen" "WSen")

if test "${USE_YMD_FORMAT}" -eq 1
then
    # x annotation interval (annotate every _ month)
    osw=$(./annotevery.py ${wesn[0]} ${wesn[1]} --time)
fi

gmt set FONT_ANNOT_PRIMARY +8p
gmt set PS_CHAR_ENCODING ISOLatin1+

it=2
for comp in up east north
do
    yshift="2.2"
    ycol="${Y[$it]}"
    if test "${it}" -eq 2
    then
        Osw=
        Ysh=
    else
        Osw="-O"
        Ysh="-Y${yshift}i"
    fi
    if test "${USE_YMD_FORMAT}" -eq 1
    then
        gmt psbasemap "${R[$it]}" -JX10i/2i -K ${Osw} \
            -Bsx1Y -Bpxa${osw}Of1o \
            -Bpy"${YAI[$it]}"+l${comp} -B"${XAT[$it]}"+t"Time Series"+glightgreen \
            ${Ysh} >> $ps
    else
        gmt psbasemap "${R[$it]}" -JX10i/2i -K ${Osw} \
            -Bsx100 \
            -Bpy"${YAI[$it]}"+l${comp} -B"${XAT[$it]}"+t"Time Series"+glightgreen \
            ${Ysh} >> $ps
    fi

    cat ${cts_file} | \
        awk -F "," -v y=$ycol '{print $1,$y}' | \
        gmt psxy "${R[$it]}" -J \
        -Wthin,red  -Gblack -Sc.05 \
        -h --IO_N_HEADER_RECS=1 -O -K >> $ps

    let it=it-1
done
