if [ $# -lt 1 ]; then
    echo "  "
    echo "  ./run.sh NJETS"
    echo "  "
    exit -1
fi


LUMINOSITY=19.468

NJETS=$1

CHANNELS=" EE MuMu EMu MuE "

SAMPLES="          \
ggWWto2L           \
WWTo2L2Nu          \
WZ                 \
ZZ                 \
TTbar              \
TW                 \
DY                 \
WgammaStar         \
WgammaNoStar       \
WJetsFakes\_Total  \
DataRun2012\_Total \
DYtautau           \
HWW125             \
Zgamma             \
"


rm -rf rootfiles/${NJETS}jet


# Loop
for CHANNEL in $CHANNELS; do

    for SAMPLE in $SAMPLES; do 
	
	root -l -b -q "runLatinosTreeScript.C($LUMINOSITY,$NJETS,\"$CHANNEL\",\"$SAMPLE\")"
  
    done

    OUTPATH=rootfiles/${NJETS}jet/${CHANNEL}

    # Merge some MC files
    hadd ${OUTPATH}/WW.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu.root 
    hadd ${OUTPATH}/Wgamma.root ${OUTPATH}/WgammaNoStar.root ${OUTPATH}/WgammaStar.root

done


# Make the SF, OF and All directories
PATHSF=rootfiles/${NJETS}jet/SF
PATHOF=rootfiles/${NJETS}jet/OF
PATHAll=rootfiles/${NJETS}jet/All

mkdir ${PATHSF}
mkdir ${PATHOF}
mkdir ${PATHAll}
    
PATHEMu=rootfiles/${NJETS}jet/EMu
PATHMuE=rootfiles/${NJETS}jet/MuE
PATHMuMu=rootfiles/${NJETS}jet/MuMu
PATHEE=rootfiles/${NJETS}jet/EE

ALLSAMPLES="       \
ggWWto2L           \
WWTo2L2Nu          \
WZ                 \
ZZ                 \
TTbar              \
TW                 \
DY                 \
WgammaStar         \
WgammaNoStar       \
WJetsFakes_Total   \
DataRun2012_Total  \
DYtautau           \
HWW125             \
WW                 \
Wgamma             \
Zgamma             \
"

for SAMPLE in $ALLSAMPLES; do 
 
    hadd ${PATHSF}/${SAMPLE}.root  ${PATHMuMu}/${SAMPLE}.root ${PATHEE}/${SAMPLE}.root
    hadd ${PATHOF}/${SAMPLE}.root  ${PATHMuE}/${SAMPLE}.root  ${PATHEMu}/${SAMPLE}.root
    hadd ${PATHAll}/${SAMPLE}.root ${PATHOF}/${SAMPLE}.root   ${PATHSF}/${SAMPLE}.root
 
done
