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
## Append model (lines) to the plots
PLOT_MODEL=0
## PLOT EVENTS
PLOT_EVENTS=0
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
        -m|--model-ascii)
        PLOT_MODEL=1
        MODEL_INFILE="$2"
        shift 2
        ;;
        -e|--plot-events)
        PLOT_EVENTS=1
        shift
        ;;
        -q|--noa-catalogue)
        NOA_CAT_FILE="$2"
        shift 2
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
    gmt set FORMAT_DATE_OUT yyyy-mm-dd
    if ! ./mjd2ymdhms.py -f ${cts_file} -s 1 -d ',' > ${cts_file}.tmp
    then
        echoerr "[ERROR] Failed to transform dates in input file."
        rm ${cts_file} 2>/dev/null
        exit 1
    fi
    mv ${cts_file}.tmp ${cts_file}
    wesn=( $(gmt info -fT -I0.01 -C @"${cts_file}" --FORMAT_DATE_IN=yyyy-mm-dd -h --IO_N_HEADER_RECS=1) )
else
    wesn=( $(gmt info -I0.01 -C @"${cts_file}" -h --IO_N_HEADER_RECS=1) )
fi
echo "[DEBUG] Here are the limits: ${wesn[0]} and ${wesn[1]}"
this_year=$(echo ${wesn[1]} | awk -F"-" '{print $1}')
next_year=$(echo ${wesn[1]} | awk -F"-" '{print $1+1}')
#wesn[1]="${next_year}-Jan-01T00:00:00"
#wesn[0]="2020-Feb-01T00:00:00"
#wesn[1]="2020-Apr-31T00:00:00"
echo "[DEBUG] Here are the limits: ${wesn[0]} and ${wesn[1]}"

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

## If we want to append the model lines, prepare the input data
if test "${PLOT_EVENTS}" -eq 1
then
    APPND_CMD="-v .evn.dat"
    if [ ! -z ${NOA_CAT_FILE+x} ]
    then
        APPND_CMD="${APPND_CMD} -q ${NOA_CAT_FILE}"
        echo "[DEBUG] Identifing earthquakes with a catalogues files may take a little longer"
        echo "        Hold on ..."
    fi
fi
if test "${PLOT_MODEL}" -eq 1
then
    if test "${USE_YMD_FORMAT}" -eq 1
    then
        if ! ./tsmodel.py -f "${MODEL_INFILE}" -s "${wesn[0]}" -e "${wesn[1]}" ${APPND_CMD} > .mld.dat
        then
            echo "[ERROR] Failed to parse model file"
            exit 1
        fi
    else
        if ! ./tsmodel.py -f "${MODEL_INFILE}" -s "${wesn[0]}" -e "${wesn[1]}" -m  ${APPND_CMD} > .mld.dat
        then
            echo "[ERROR] Failed to parse model file"
            exit 1
        fi
    fi
fi

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
        echo "R option=${R[$it]}"
        gmt psbasemap "${R[$it]}" -JX10i/2i -K ${Osw} \
            -Bsxpa7Rf1d -Bpxa1OS \
            -Bpy"${YAI[$it]}"+l${comp} -B"${XAT[$it]}"+t"Time Series"+glightgreen \
            ${Ysh} >> $ps
    else
        gmt psbasemap "${R[$it]}" -JX10i/2i -K ${Osw} \
            -Bsx100 \
            -Bpy"${YAI[$it]}"+l${comp} -B"${XAT[$it]}"+t"Time Series"+glightgreen \
            ${Ysh} >> $ps
    fi
    
    if test "${PLOT_EVENTS}" -eq 1
    then
        echo "[DEBUG] Ploting events ..."
        c_min=$(echo "${R[$it]}" | awk -F "/" '{print $3}')
        c_max=$(echo "${R[$it]}" | awk -F "/" '{print $4}')
        evncols=("brown" "orange" "yellow")
        evncit=0
        eq_text_dy=(0.0 0.3)
        for evn in j e v
        do
            if grep " $evn" .evn.dat > .evn.dat${evn} 2>/dev/null
            then
                eq_it=0
                while read line
                do
                    #t=$(echo $line | sed "s/${evn}//g" | sed "s/ $//g" | tr ' ' 'T')
                    t=$(echo $line | awk '{print $1"T"$2}' | sed 's/\..*$//g')
                    magn=$(echo $line | awk '{print $NF}')
                    echo "${t} ${c_min}" >  .evn.dat${evn}${evn}
                    echo "${t} ${c_max}" >> .evn.dat${evn}${evn}
                    cat .evn.dat${evn}${evn} | gmt psxy "${R[$it]}" -J -W0.3p,"${evncols[$evncit]}" -O -K >> $ps
                    if test "${evn}" = "e" && test "${comp}" = "up"
                    then
                        magn_upper=$(echo "$magn" | tr 'a-z' 'A-Z')
                        if test "${magn_upper}" != "NONE"
                        then
                            crd=$(echo "$c_min + ($c_max - $c_min)/8" | bc -l)
                            dy_idx=$(echo $((eq_it%2)))
                            echo "${t} ${crd} ${magn}" | \
                            gmt pstext "${R[$it]}" -J -Dj0c/"${eq_text_dy[$dy_idx]}"c \
                            -F+a0+f7,Helvetica+jCT -O -K >> $ps
                            let eq_it=eq_it+1
                        else
                            echo "[DEBUG] No magnitude for earthquake at $t; omitting ..."
                        fi
                    fi
                done < .evn.dat${evn}
            fi
            let evncit=evncit+1
        done
    fi

    cat ${cts_file} | \
        awk -F "," -v y=$ycol '{print $1,$y}' | \
        gmt psxy "${R[$it]}" -J \
        -Wthin,red  -Gblack -Sc.13 \
        -h --IO_N_HEADER_RECS=1 -O -K >> $ps

    if test "${PLOT_MODEL}" -eq 1
    then
        mcol=$((2+${it}))
        cat .mld.dat | \
        awk -v y=$mcol '{print $1,$y}' | \
        gmt psxy "${R[$it]}" -J \
        -Wthin,blue -O -K >> $ps
    fi

    let it=it-1
done

rm ${cts_file} 2>/dev/null
exit 0
