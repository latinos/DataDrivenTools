CHANNELS=" All SF OF EE MuMu EMu MuE "

SAMPLES="         \
ggWWto2L          \
WWTo2L2Nu         \
WZ                \
ZZ                \
TTbar             \
TW                \
DY                \
WgammaStar        \
WgammaNoStar      \
WJetsFakes_Total  \
DataRun2012_Total \
DYtautau          \
HWW125            \
"


rm -rf rootfiles/3jet
mkdir rootfiles/3jet


# Loop
for CHANNEL in $CHANNELS; do

    OUTPATH=rootfiles/3jet/${CHANNEL}

    mkdir ${OUTPATH}

    for SAMPLE in $SAMPLES; do 
	
	hadd ${OUTPATH}/${SAMPLE}.root \
	    rootfiles/0jet/${CHANNEL}/${SAMPLE}.root \
	    rootfiles/1jet/${CHANNEL}/${SAMPLE}.root \
	    rootfiles/2jet/${CHANNEL}/${SAMPLE}.root 
  
    done


    # Merge some MC files
    hadd ${OUTPATH}/WW.root     ${OUTPATH}/ggWWto2L.root     ${OUTPATH}/WWTo2L2Nu.root 
    hadd ${OUTPATH}/Wgamma.root ${OUTPATH}/WgammaNoStar.root ${OUTPATH}/WgammaStar.root

done
