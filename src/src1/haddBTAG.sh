#! /bin/bash
LOCATION=output2013MCStopUser19Nuevo_FastSim/res
LIST=$(ls $LOCATION/histo*.root)
echo "lista de archivos root"
FLAG="JP"
for FILE in $LIST; do
    if [ $FLAG == "JP" ]; then
	echo $FILE
	echo "copiando $FILE a temp"
	rm $LOCATION/temp.root
	cp $FILE $LOCATION/temp.root
	FLAG="PV"
	echo "lista de archivos root, ahora debe estar temp"
	ls $LOCATION/*.root
    else
	echo "Haciendo hadd entre temp y archivo"
	rm  $LOCATION/salida.root
	echo $FILE
	hadd -f  $LOCATION/salida.root $FILE $LOCATION/temp.root
	echo "copiando resultado adicion a temp"
	cp $LOCATION/salida.root $LOCATION/temp.root

	
    fi
    
    
done
exit 0

