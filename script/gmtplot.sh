#! /usr/bin/bash

echoerr() { echo "$@" 1>&2; }

##  Use YYYY-MM-DD format (for plotting)
USE_YMD_FORMAT=0
##  Do not plot outliers
REMOVE_OUTLIERS=0
##  current pid
pid=$(echo $$)
## output file
ps=koko.ps
rm $ps > /dev/null

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
tmpf=".${cts_file}.${pid}"
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
    ./mjd2ymdhms.py -f ${cts_file} -s 1 -d ',' > ${cts_file}.tmp
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
R=("${Rx}" "${Ry}" "${Rz}")

gmt set FONT_ANNOT_PRIMARY +10p
gmt set PS_CHAR_ENCODING ISOLatin1+
echo $Rx

it=0
for comp in north east up
do
    Yshift=$(( ${it} * 3 ))
    if test "${USE_YMD_FORMAT}" -eq 1
    then
        gmt psbasemap "${R[$it]}" -JX9i/2i -O -K -Bsx1Y -Bpxa6Of1o -Bpy0.05 \
            -BWSen+t"Time Series"+glightgreen -Y${Yshift}i >> $ps
    else
        gmt psbasemap "${R[$it]}" -JX9i/2i -O -K -Bsx100 -Bpy0.05 \
            -BWSen+t"Time Series"+glightgreen -Y${Yshift}i >> $ps
    fi

    cat ${cts_file} | \
        awk -F "," '{print $1,$2}' | \
        gmt psxy "${R[$it]}" -J \
        -Wthin,red  -Gblack -Sc.05 \
        -h --IO_N_HEADER_RECS=1 -O -K >> $ps

    let it=it+1
done
cat ${cts_file} | head -2 | awk -F "," '{print $1,$2}' | gmt psxy "${R[2]}" -J -Wthin,red  -Gblack -Sc.05  -h --IO_N_HEADER_RECS=1 -O >> $ps

echo "cts file was: ${cts_file}"
