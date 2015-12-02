#!/bin/bash

runnum=('F_10P_996_9999')
sourcedir=~/data2/ProductionRuns/$runnum/
destdir=~/data2/rundata/$runnum/

mkdir $destdir
echo Created directory $destdir

echo Moving files from $sourcedir ...

mv "$sourcedir"mfixconst "$destdir"mfixconst_${runnum[x]}
cp "$sourcedir"mfix.dat "$destdir"mfix_${runnum[x]}.dat

mv "$sourcedir"EP_G "$destdir"EP_G
mv "$sourcedir"T_G "$destdir"T_G
mv "$sourcedir"U_G "$destdir"U_G
mv "$sourcedir"V_G "$destdir"V_G
mv "$sourcedir"W_G "$destdir"W_G
mv "$sourcedir"RO_G "$destdir"RO_G
mv "$sourcedir"ROP_S1 "$destdir"ROP_S1
mv "$sourcedir"ROP_S2 "$destdir"ROP_S2
mv "$sourcedir"ROP_S3 "$destdir"ROP_S3

echo Finished moving files from $sourcedir to $destdir

exit
