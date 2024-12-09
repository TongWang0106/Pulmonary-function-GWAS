#! /bin/bash

# using common version!
PLINK=/working/genepi/software/bin/plink_1.90b4
GEMMA=~/bin/gemma0941

HOMEDIR=/mnt/lustre/home/working/qingdao
WORKDIR=/working/qdu
GVER=bm25e

RUN=qd2

DATADIR=/working/qdu/BMSW201511
FILE=bm_25E_fianl
FAMPED=set1.dat

GRMDIR=/working/qdu/bmsw/grm
PLKDIR=/working/qdu/bmsw/plnk

# fev.unx: N=278 & VAR N=3 & COV N=4|5 + 5PCs
# ID sex01 sme smn age ht sm3
# fev fvc evr

USE=fev
PREF=imp2
GVER=1KP3
PLKDIR=/working/qdu/impute2/imp2/$GVER
KLST=$GVER.60k.lst

DIRUSE=fev_${PREF}
cd $DIRUSE
INFDIR=/working/qdu/$RUN/$DIRUSE
DIRRESULT=out_${PREF}
OUTDIR=${PREF}_out
PDTDIR=${PREF}_pdt
VGSDIR=${PREF}_vgs2


# job8: use info and pack

USR=wnt

mkdir -p use
cd use

mkdir -p $PREF
cd $PREF

JOB=vgs2
mkdir -p $JOB
cd $JOB
USEDIR1=/working/qdu/$RUN/$DIRUSE/${PREF}_vgs2/out_plot
awk '{ print $2 }' $INFDIR/$PREF.inf >var.$USR.lst
for VAR in `cat var.$USR.lst` 
do
echo "get info for $USR pref $PREF job $JOB var $VAR "
cp $USEDIR1/${JOB}_${VAR}.out .
cp $USEDIR1/${JOB}.${VAR}.top .
cp $USEDIR1/${JOB}.${VAR}.pmh.ps .
cp $USEDIR1/${JOB}.${VAR}.qq.png .
#var
done
rm -f var.$USR.lst
cd ..

JOB=plot
mkdir -p $JOB
cd $JOB
USEDIR1=/working/qdu/$RUN/$DIRUSE/plot_use

PMV=p3a
echo "get info for $USR pref $PREF job $JOB plot $PMV "
cp $USEDIR1/gem_${PMV}.mh.pdf .
cp $USEDIR1/gem_${PMV}.mh.ps .
cp $USEDIR1/gem_${PMV}.qq.png .

awk '{ print $2 }' $INFDIR/$PREF.inf >var.$USR.lst
for VAR in `cat var.$USR.lst` 
do
echo "get info for $USR pref $PREF job $JOB var $VAR "
cp $USEDIR1/miami_vgs2/miami_${VAR}.png .
cp $USEDIR1/miami_vgs2/miami_${VAR}.ps .
#var
done
rm -f var.$USR.lst
cd ..

JOB=lcz
mkdir -p $JOB
cd $JOB
USEDIR1=/working/qdu/$RUN/$DIRUSE/${PREF}_lcz

cp $USEDIR1/topsnp.$PREF.txt .
cp $USEDIR1/*.lcz .

cd ..

JOB=psc
mkdir -p $JOB
cd $JOB
USEDIR1=/working/qdu/$RUN/$DIRUSE/${PREF}_psc/out_plot
awk '{ print $2 }' $INFDIR/$PREF.inf >var.$USR.lst
for VAR in `cat var.$USR.lst` 
do
echo "get info for $USR pref $PREF job $JOB var $VAR "
echo "CHR PathName chi2Pvalue empPvalue -log10(x2P) -log10(emP) VAR " >hd
cat hd $USEDIR1/${JOB}_${VAR}.path >${JOB}.${VAR}.path.txt
#var
done
rm -f var.$USR.lst hd
cd ..

JOB1=gemma
mkdir -p $JOB1
cd $JOB1
USEDIR1=/working/qdu/$RUN/$DIRUSE/${PREF}_out
USEDIR2=/working/qdu/$RUN/$DIRUSE/${PREF}_pdt
awk '{ print $2 }' $INFDIR/$PREF.inf >var.$USR.lst
for VAR in `cat var.$USR.lst` 
do
echo "get info for $USR pref $PREF job $JOB1 var $VAR "
cp $USEDIR1/$PREF.$VAR.hd .
cp $USEDIR1/$PREF.$VAR.gz .
cp $USEDIR1/$PREF.$VAR.1.1.assoc.txt.gz .
cp $USEDIR2/$PREF.$VAR.pmh.ps .
cp $USEDIR2/$PREF.$VAR.qq.png .
#var
done
rm -f var.$USR.lst
cd ..

cd ..

tar -czvf $USR.$PREF.tgz $PREF/*
tar -tzvf $USR.$PREF.tgz >$USR.$PREF.tgz.lst



exit



# job7a: plot_use GEMMA p3a [fev fvc evr]

cd plot_use

JOB=gem_p3a
mkdir -p $JOB
cd $JOB

REFDIR=/working/qdu/u1_imp2/plot_use
P3PSFILE=p3.gwas.tmp
P3PDFFILE=p3.gwas.pdf.tmp
P3QQFILE=p3.qq.R.tmp

RUNDIR=/working/qdu/$RUN/$DIRUSE/plot_use/$JOB

PLT=${JOB}
PLTFILE=$PLT
PDFFILE=$PLT.pdf

VAR=$JOB
VAR1=fev
VAR2=fvc
VAR3=evr

NAME1=`awk -v var=$VAR1 '$2==var { print $4 }' $INFDIR/$PREF.inf`
NAME2=`awk -v var=$VAR2 '$2==var { print $4 }' $INFDIR/$PREF.inf`
NAME3=`awk -v var=$VAR3 '$2==var { print $4 }' $INFDIR/$PREF.inf`
NUM1=`awk -v var=$VAR1 '$2==var { print $3 }' $INFDIR/$PREF.inf`
NUM2=`awk -v var=$VAR2 '$2==var { print $3 }' $INFDIR/$PREF.inf`
NUM3=`awk -v var=$VAR3 '$2==var { print $3 }' $INFDIR/$PREF.inf`

USEDIR1=/working/qdu/$RUN/$DIRUSE/${PREF}_pdt

echo "plot for var $VAR with $NAME1, $NAME2 & $NAME3  "

PROG1=GEMMA-IMPUTED
PROG2=GEMMA-IMPUTED
PROG3=GEMMA-IMPUTED

sort -nr -k3,3 $USEDIR1/$PREF.$VAR1.Rqq.tail  >u1.tail
sort -nr -k3,3 $USEDIR1/$PREF.$VAR2.Rqq.tail  >u2.tail
sort -nr -k3,3 $USEDIR1/$PREF.$VAR3.Rqq.tail  >u3.tail

cp $USEDIR1/$PREF.$VAR1.qq.lmd  u1.lmd
cp $USEDIR1/$PREF.$VAR2.qq.lmd  u2.lmd
cp $USEDIR1/$PREF.$VAR3.qq.lmd  u3.lmd

TOPSNP1=`awk 'NR==1 {print $5 }' u1.tail `
SNPLNG1=`awk 'NR==1 {print $1 }' u1.tail `
SNPLGP1=`awk 'NR==1 {print $3 }' u1.tail `
SNPNUM1=`tail -1 u1.lmd | awk '{print $2 }'`

TOPSNP2=`awk 'NR==1 {print $5 }' u2.tail `
SNPLNG2=`awk 'NR==1 {print $1 }' u2.tail `
SNPLGP2=`awk 'NR==1 {print $3 }' u2.tail `
SNPNUM2=`tail -1 u2.lmd | awk '{print $2 }'`

TOPSNP3=`awk 'NR==1 {print $5 }' u3.tail `
SNPLNG3=`awk 'NR==1 {print $1 }' u3.tail `
SNPLGP3=`awk 'NR==1 {print $3 }' u3.tail `
SNPNUM3=`tail -1 u3.lmd | awk '{print $2 }'`

sed -e "s/REPFILE/$PLTFILE/g" \
    -e "s/REPMAXY/13/g" \
    -e "s/REPPROG1/$PROG1/g" \
    -e "s/REPPROG2/$PROG2/g" \
    -e "s/REPPROG3/$PROG3/g" \
    -e "s/REPNAME1/$NAME1/g" \
    -e "s/REPNAME2/$NAME2/g" \
    -e "s/REPNAME3/$NAME3/g" \
    -e "s/REPNUM1/$NUM1/g" \
    -e "s/REPNUM2/$NUM2/g" \
    -e "s/REPNUM3/$NUM3/g" \
    -e "s/REPSNPN1/$SNPNUM1/g"  \
    -e "s/REPSNPN2/$SNPNUM2/g"  \
    -e "s/REPSNPN3/$SNPNUM3/g"  \
    -e "s/REPSNP1/$TOPSNP1/g"  \
    -e "s/REPSNP2/$TOPSNP2/g"  \
    -e "s/REPSNP3/$TOPSNP3/g"  \
    -e "s/REPPOS1/$SNPLNG1/g" \
    -e "s/REPPOS2/$SNPLNG2/g" \
    -e "s/REPPOS3/$SNPLNG3/g" \
    -e "s/REPTOP3/$SNPLGP3/g"  \
    -e "s/REPTOP1/$SNPLGP1/g"  \
    -e "s/REPTOP2/$SNPLGP2/g"  $REFDIR/$P3PSFILE >$PLTFILE.plt

sed -e "s/REPFILE/$PDFFILE/g" \
    -e "s/REPMAXY/13/g" \
    -e "s/REPPROG1/$PROG1/g" \
    -e "s/REPPROG2/$PROG2/g" \
    -e "s/REPPROG3/$PROG3/g" \
    -e "s/REPNAME1/$NAME1/g" \
    -e "s/REPNAME2/$NAME2/g" \
    -e "s/REPNAME3/$NAME3/g" \
    -e "s/REPNUM1/$NUM1/g" \
    -e "s/REPNUM2/$NUM2/g" \
    -e "s/REPNUM3/$NUM3/g" \
    -e "s/REPSNPN1/$SNPNUM1/g"  \
    -e "s/REPSNPN2/$SNPNUM2/g"  \
    -e "s/REPSNPN3/$SNPNUM3/g"  \
    -e "s/REPSNP1/$TOPSNP1/g"  \
    -e "s/REPSNP2/$TOPSNP2/g"  \
    -e "s/REPSNP3/$TOPSNP3/g"  \
    -e "s/REPPOS1/$SNPLNG1/g" \
    -e "s/REPPOS2/$SNPLNG2/g" \
    -e "s/REPPOS3/$SNPLNG3/g" \
    -e "s/REPTOP3/$SNPLGP3/g"  \
    -e "s/REPTOP1/$SNPLGP1/g"  \
    -e "s/REPTOP2/$SNPLGP2/g"  $REFDIR/$P3PDFFILE >$PDTFILE.plt

sed -e "s/REPFILE/$PLTFILE/g" \
    -e "s/REPPROG1/$PROG1/g" \
    -e "s/REPPROG2/$PROG2/g" \
    -e "s/REPPROG3/$PROG3/g" \
    -e "s/REPNAME1/$NAME1/g" \
    -e "s/REPNAME2/$NAME2/g" \
    -e "s/REPNAME3/$NAME3/g" \
    -e "s/RPRQQ1/u1/g" \
    -e "s/RPRQQ2/u2/g" \
    -e "s/RPRQQ3/u3/g"  $REFDIR/$P3QQFILE >$PLTFILE.R

JOBRUN=$VAR.sub

cat << __EOT__  >$JOBRUN
#PBS -l ncpus=1
#PBS -l mem=1gb
#PBS -l walltime=40:00:00 
module load R
cd ${RUNDIR}

cp $USEDIR1/$VAR1.$PREF.ods  u1.ods
cp $USEDIR1/$VAR1.$PREF.evn  u1.evn
cp $USEDIR1/$PREF.$VAR1.Rqq  u1.Rqq
cp $USEDIR1/fil2p.ods  u1p.ods
cp $USEDIR1/fil2p.evn  u1p.evn

cp $USEDIR1/$VAR2.$PREF.ods  u2.ods
cp $USEDIR1/$VAR2.$PREF.evn  u2.evn
cp $USEDIR1/$PREF.$VAR2.Rqq  u2.Rqq
cp $USEDIR1/fil2p.ods  u2p.ods
cp $USEDIR1/fil2p.evn  u2p.evn

cp $USEDIR1/$VAR3.$PREF.ods  u3.ods
cp $USEDIR1/$VAR3.$PREF.evn  u3.evn
cp $USEDIR1/$PREF.$VAR3.Rqq  u3.Rqq
cp $USEDIR1/fil2p.ods  u3p.ods
cp $USEDIR1/fil2p.evn  u3p.evn


gnuplot $PLTFILE.plt

gnuplot $PDTFILE.plt

mv -f $PLTFILE.ps ../$PLTFILE.mh.ps
mv -f $PDFFILE.pdf ../$PLTFILE.mh.pdf

R CMD BATCH $PLTFILE.R

mv -f $PLTFILE.png ../$PLTFILE.qq.png
cd ..
rm -f -r $VAR

__EOT__

qsub $JOBRUN


exit


# job7:  making lcz files

LCZDIR=${PREF}_lcz
mkdir -p $LCZDIR
cd $LCZDIR

USEDIR1=/working/qdu/$RUN/$DIRUSE/${PDTDIR}
USEDIR2=/working/qdu/$RUN/$DIRUSE/${OUTDIR}

awk '{ print $2 }' $INFDIR/$PREF.inf >var.lst

# choose top 3 CHRs!!!

TOPFILE=topsnp.imp2.txt
>$TOPFILE
for VAR in `cat var.lst` 
do
sort -nr -k3,3 $USEDIR1/$PREF.$VAR.top | awk -v var=$VAR '{ print $0, var }' >>$TOPFILE

mkdir -p $VAR
cd $VAR

RUNDIR=/working/qdu/$RUN/$DIRUSE/$LCZDIR/$VAR

sort -nr -k3,3 $USEDIR1/$PREF.$VAR.top  >top.$VAR.txt
awk 'NR<4 { print $2 }' top.$VAR.txt >chr.$VAR.lst

JOBSUB=lcz.$VAR.sub

cat << __EOT__  >$JOBSUB
#PBS -l ncpus=1
#PBS -l mem=1gb
#PBS -l walltime=40:00:00 
cd $RUNDIR
for CHR in \`cat chr.$VAR.lst \`
do
zless $USEDIR2/$PREF.$VAR.use.gz | awk -v chr=\$CHR 'BEGIN {print "SNP P CHR BP LGP ALE"
                      }
                      { if(\$6==chr && \$10>=0.01 && \$10<=0.99 ) print \$2, \$(NF-1), \$6, \$7, -log(\$(NF-1))*0.434294482, \$8"_"\$9"_"\$10
}' >$PREF.$VAR.\$CHR.lcz
mv -f $PREF.$VAR.\$CHR.lcz ../

#chr
done
__EOT__


qsub $JOBSUB

cd ..

#var
done


exit




# job6: plot_use miami GEMMA & VEGAS2

mkdir -p plot_use
cd plot_use

REFDIR=~/working/ref
PMTFILE=miami.gwas.tmp
PQQFILE=miami.qq.R.tmp

JOB=miami_vgs2
mkdir -p $JOB
cd $JOB

VGSDIR=${PREF}_vgs2

# miami for gwas & gene-base

RUNDIR=/working/qdu/$RUN/$DIRUSE/plot_use/$JOB

awk '{ print $2 }' $INFDIR/$PREF.inf >var.lst
for VAR in `cat var.lst` 
do

mkdir -p $VAR
cd $VAR

RUNDIR=/working/qdu/$RUN/$DIRUSE/plot_use/$JOB/$VAR
USEDIR1=/working/qdu/$RUN/$DIRUSE/${PDTDIR}
USEDIR2=/working/qdu/$RUN/$DIRUSE/${VGSDIR}/vgs2_${VAR}

NAME=`awk -v var=$VAR '$2==var { print $4 }' $INFDIR/$PREF.inf`

NUM1=`awk -v var=$VAR '$2==var { print $3 }' $INFDIR/$PREF.inf`
NUM2=`wc -l $USEDIR2/vgs2.$VAR.Rqq | awk '{ print $1 }'`

NAME1="$NAME (GWAS Sample-N=${NUM1})"
NAME2="$NAME (GENE-BASE Gene-N=${NUM2})"

PROG1=GEMMA
PROG2=VEGAS2

USE1="SNP N = "
USE2="GENE N = "

PLTFILE=miami_${VAR}

sort -nr -k3,3 $USEDIR1/$PREF.$VAR.Rqq.tail  >u1.tail
cp $USEDIR1/$PREF.$VAR.qq.lmd  u1.lmd
awk 'NR==1 {print }' u1.tail >u1.tail.1

sort -nr -k3,3 $USEDIR2/vgs2.$VAR.Rqq.head  >u2.tail
cp $USEDIR2/vgs2.$VAR.qq.lmd  u2.lmd
awk 'NR==1 {print }' u2.tail >u2.tail.1

# read u2.Rqq same vars (6)!! file with 9 (lgpSNP snpP SNP!!!)
awk '{ print $1, $2, $3, $4, $5, $6 }' $USEDIR2/vgs2.$VAR.Rqq  >u2.Rqq

TOPSNP1=`awk 'NR==1 {print $5 }' u1.tail `
SNPLNG1=`awk 'NR==1 {print $1 }' u1.tail `
SNPLGP1=`awk 'NR==1 {print $3 }' u1.tail `
SNPNUM1=`tail -1 u1.lmd | awk '{print $2 }'`

TOPSNP2=`awk 'NR==1 {print $5 }' u2.tail `
SNPLNG2=`awk 'NR==1 {print $1 }' u2.tail `
SNPLGP2=`awk 'NR==1 {print $3 }' u2.tail `
SNPNUM2=`tail -1 u2.lmd | awk '{print $2 }'`

sed -e "s/REPFILE/$PLTFILE/g" \
    -e "s/REPMAXY/20/g" \
    -e "s/REPPROG1/$PROG1/g" \
    -e "s/REPPROG2/$PROG2/g" \
    -e "s/REPNAME1/$NAME1/g" \
    -e "s/REPNAME2/$NAME2/g" \
    -e "s/REPNUM1/$NUM1/g" \
    -e "s/REPNUM2/$NUM2/g" \
    -e "s/REPSNPN1/$SNPNUM1/g"  \
    -e "s/REPSNPN2/$SNPNUM2/g"  \
    -e "s/REPSNP1/$TOPSNP1/g"  \
    -e "s/REPSNP2/$TOPSNP2/g"  \
    -e "s/REPPOS1/$SNPLNG1/g" \
    -e "s/REPPOS2/$SNPLNG2/g" \
    -e "s/REPTOP1/$SNPLGP1/g"  \
    -e "s/REPTOP2/$SNPLGP2/g"  $REFDIR/$PMTFILE >$PLTFILE.plt

sed -e "s/REPFILE/$PLTFILE/g" \
    -e "s/REPPROG1/$PROG1/g" \
    -e "s/REPPROG2/$PROG2/g" \
    -e "s/REPNAME1/$NAME1/g" \
    -e "s/REPNAME2/$NAME2/g" \
    -e "s/REPUSE1/$USE1/g" \
    -e "s/REPUSE2/$USE2/g" \
    -e "s/RPRQQ1/u1/g" \
    -e "s/RPRQQ2/u2/g"  $REFDIR/$PQQFILE >$PLTFILE.R


JOBRUN=$VAR.sub

cat << __EOT__  >$JOBRUN
#PBS -l ncpus=1
#PBS -l mem=4gb
#PBS -l walltime=40:00:00 
module load R
cd ${RUNDIR}

cp $USEDIR1/$VAR.$PREF.ods  u1.ods
cp $USEDIR1/$VAR.$PREF.evn  u1.evn
cp $USEDIR1/fil2p.ods  u1p.ods
cp $USEDIR1/fil2p.evn  u1p.evn
cp $USEDIR1/$PREF.$VAR.Rqq  u1.Rqq

cp $USEDIR2/$VAR.vgs2.ods  u2.ods
cp $USEDIR2/$VAR.vgs2.evn  u2.evn

gnuplot $PLTFILE.plt

R CMD BATCH $PLTFILE.R
#rm -f .RData
#rm -f u*
mv -f $PLTFILE.pdf ../
mv -f $PLTFILE.ps ../
mv -f $PLTFILE.png ../
cd ..
rm -f -r $VAR

__EOT__

qsub $JOBRUN

cd ..

done



exit




# job5a: psc out & plot

JOBDIR=${PREF}_psc
cd $JOBDIR

USE=psc
REFDIR=~/working/ref

RUNDIR=/working/qdu/$RUN/$DIRUSE/${JOBDIR}

# dir for plots!!!

PSCOPT=out_plot
mkdir -p $PSCOPT

PROG=PASCAL-ASN

awk '{ print $2 }' $INFDIR/$PREF.inf >var.lst
for VAR in `cat var.lst` 
do

NAME=`awk -v var=$VAR '$2==var { print $4 }' $INFDIR/$PREF.inf`

cd ASN_${VAR}

RUNDIR=/working/qdu/$RUN/$DIRUSE/${JOBDIR}/ASN_${VAR}

OUT=${USE}_${VAR}
JOBSUB=$USE.$VAR.sub
JOBRUN=$USE.$VAR.run
JOBOUT=${USE}_${VAR}.log

cat << __EOT__  >$JOBSUB
#PBS -l ncpus=1
#PBS -l mem=1gb
#PBS -l walltime=40:00:00 
module load R
cd $RUNDIR || exit 
(
./$JOBRUN > $JOBOUT  2>&1 
) &
wait
__EOT__

sed -e "s/REPRUN/$USE/g" \
    -e "s/REPVAR/$VAR/g" \
    -e "s/REPOUT/$OUT/g" \
    -e "s/REPLAB/$NAME/g" \
    -e "s/REPDIR/$PSCOPT/g" \
    -e "s/REPPROG/$PROG/g"  $REFDIR/psc_opdt.tmp  >$JOBRUN

chmod +x $JOBRUN
qsub $JOBSUB


cd ..

#var
done


cd ..


exit




# job4a: vgs out & plot

VGSDIR=${PREF}_vgs2
cd $VGSDIR

REFDIR=/working/qdu/u1_dbsnp/dbsnp_vgs2

USE=vgs2

RUNDIR=/working/qdu/$RUN/$DIRUSE/${VGSDIR}

rm -f *.easia  *.easia.e* *.easia.o*

VGSOPT=out_plot
mkdir -p $VGSOPT

PROG=vegas2-ASN

for k in `cat k.lst`
do

VAR=`awk -v kk=$k '$1==kk { print $2 }' $INFDIR/$PREF.inf`
NAME=`awk -v kk=$k '$1==kk { print $4 }' $INFDIR/$PREF.inf`

OUT=${USE}_${VAR}
JOBSUB=$USE.$VAR.sub
JOBRUN=$USE.$VAR.run
JOBOUT=${USE}_${VAR}.log

cat << __EOT__  >$JOBSUB
#PBS -l ncpus=1
#PBS -l mem=1gb
#PBS -l walltime=40:00:00 
module load R
cd $RUNDIR || exit 
(
./$JOBRUN > $JOBOUT  2>&1 
) &
wait
__EOT__

sed -e "s/REPRUN/$USE/g" \
    -e "s/REPVAR/$VAR/g" \
    -e "s/REPOUT/$OUT/g" \
    -e "s/REPLAB/$NAME/g" \
    -e "s/REPDIR/$VGSOPT/g" \
    -e "s/REPFILE/$ASSOC/g" \
    -e "s/REPPROG/$PROG/g"  $REFDIR/vgs_opdt.tmp  >$JOBRUN

chmod +x $JOBRUN
qsub $JOBSUB


# k
done

cd ..



exit



#  job5:  pathway PASCAL

# for VAR in fev fvc evr

JOBDIR=${PREF}_psc
mkdir -p $JOBDIR
cd $JOBDIR

PASCAL=/working/qdu/bin/PASCAL/Pascal
PSCDIR=/working/qdu/bin/PASCAL

# run ASN pathway

#awk '{ print $2 }' $INFDIR/$PREF.inf >var.lst

for VAR in fev fvc evr
do
mkdir -p ASN_${VAR}
done

USEDIR=/working/qdu/$RUN/$DIRUSE/${OUTDIR}
USE=psc
RUNDIR=/working/qdu/$RUN/$DIRUSE/${JOBDIR}
JOBRUN=$USE.ASN.run

cat << __EOT__  >$JOBRUN
#PBS -l ncpus=1
#PBS -l mem=4gb
#PBS -l walltime=40:00:00 
module load vegas/2
cd $RUNDIR
for VAR in fev fvc evr
do

OUTDIR=/working/qdu/$RUN/$DIRUSE/$JOBDIR/ASN_\${VAR}
cd ${PSCDIR} 
for ((CHR=1;CHR<=22;CHR++))
do

# only rs,  prb
zless ${USEDIR}/$PREF.\$VAR.use.gz | awk -v chr=\$CHR '\$6==chr && substr(\$2,1,2)== "rs" {print \$2, \$(NF-1) }'  >\$CHR.tmp
sed "s/ /\t/g" \$CHR.tmp >\$VAR.\$CHR.txt
rm -f \$CHR.tmp 

./Pascal --pval=\$VAR.\$CHR.txt --chr=\$CHR --runpathway=on --custom=ASN --customdir=resources/ASN/ --outsuffix=ASN

cd output
mv -f \${VAR}.\${CHR}.*  \$OUTDIR/
cd ..
rm -f \$VAR.\$CHR.txt

#chr
done

cd $RUNDIR

#var
done

__EOT__

qsub $JOBRUN

# above run one by one !!!


exit





# job4:  gene-based vegas2

VGSDIR=${PREF}_vgs2
mkdir -p $VGSDIR
cd $VGSDIR

mkdir -p out

VEGAS=vegas2
RUNDIR=/working/qdu/$RUN/$DIRUSE/${VGSDIR}
USEDIR=/working/qdu/$RUN/$DIRUSE/${OUTDIR}

awk '{ print $1 }' $INFDIR/$PREF.inf >k.lst

for k in `cat k.lst`
do

VAR=`awk -v kk=$k '$1==kk { print $2 }' $INFDIR/$PREF.inf`

for ((CHR=1;CHR<=23;CHR++))
do

JOBRUN=vgs.$VAR.$CHR.easia

cat << __EOT__  >$JOBRUN
#PBS -l ncpus=1
#PBS -l mem=1gb
#PBS -l walltime=40:00:00 
module load vegas/2
cd ${RUNDIR} 
mkdir -p ${VAR}_${CHR}_EASIA
cd ${VAR}_${CHR}_EASIA

# only rs,  prb
zless ${USEDIR}/$PREF.$VAR.use.gz | awk -v chr=$CHR '\$6==chr && substr(\$2,1,2)== "rs" {print \$2, \$(NF-1) }'  >$CHR.tmp
sed "s/ /\t/g" $CHR.tmp >$VAR.$CHR.txt
rm -f $CHR.tmp 
wc -l $VAR.$CHR.txt 

$VEGAS $VAR.$CHR.txt  -pop 1000GASN -subpop ASN -chr $CHR -out ${VAR}-${CHR}.out

mv -f ${VAR}-${CHR}.out ../out/
cd ..
rm -f -r ${VAR}_${CHR}_EASIA
__EOT__

qsub $JOBRUN

#chr
done

#var
done


exit




#  job3: plots

REFDIR=~/working/ref
DIRUSE=fev_${PREF}
cd $DIRUSE

INFDIR=/working/qdu/$RUN/$DIRUSE

DIRRESULT=out_${PREF}
OUTDIR=${PREF}_out

PDTDIR=${PREF}_pdt
mkdir -p $PDTDIR
cd $PDTDIR

RUNDIR=/working/qdu/$RUN/$DIRUSE/$PDTDIR
DIR="..\/..\/${OUTDIR}"

awk '{ print $1 }' $INFDIR/$PREF.inf >k.lst

for k in `cat k.lst`
do

VAR=`awk -v kk=$k '$1==kk { print $2 }' $INFDIR/$PREF.inf`
NAME=`awk -v kk=$k '$1==kk { print $4 }' $INFDIR/$PREF.inf`

XYFILE=$PREF.$VAR.XYFILE
cp ../$OUTDIR/$XYFILE .

ASSOC=$PREF.$VAR.use.gz
PROG=GEMMA

JOBSUB=$PREF.$VAR.sub
JOBRUN=$PREF.$VAR.run
JOBOUT=${PREF}_${VAR}.log

cat << __EOT__  >$JOBSUB
#PBS -l ncpus=1
#PBS -l mem=4gb
#PBS -l walltime=40:00:00 
module load R
cd $RUNDIR || exit 
(
./$JOBRUN > $JOBOUT  2>&1 
) &
wait
__EOT__

sed -e "s/REPRUN/$PREF/g" \
    -e "s/REPVAR/$VAR/g" \
    -e "s/REPLAB/$NAME/g" \
    -e "s/REPDIR/$DIR/g" \
    -e "s/REPFILE/$ASSOC/g" \
    -e "s/REPNAME/$NAME/g" \
    -e "s/REPPROG/$PROG/g" \
    -e "s/REPXYF/$XYFILE/g"  $REFDIR/gem.$GVER.23.pdt.use  >$JOBRUN

chmod +x $JOBRUN
qsub $JOBSUB

# k
done


exit




# 2: get outs

REFDIR=/working/qdu/qd1

DIRUSE=fev_${PREF}
cd $DIRUSE

DIRRESULT=out_${PREF}

cp $PLKDIR/$KLST .
rm -f *.cb *.cb.e* *.cb.o* 

# sample N=278
#  fev: VAR: fev fvc evr
#  cov9: [sex,age,ht,sm3 & 5PCs]

cat << __EOT__  >$PREF.inf
1 fev 278 FEV1
2 fvc 278 FVC
3 evr 278 FEV1FVC
__EOT__


OUTDIR=${PREF}_out
mkdir -p $OUTDIR
cd $OUTDIR

awk '{ print $1 }' ../$PREF.inf >k.lst

for k in `cat k.lst`
do

VAR=`awk -v kk=$k '$1==kk { print $2 }' ../$PREF.inf`
NUM=`awk -v kk=$k '$1==kk { print $3 }' ../$PREF.inf`

FILE=$PREF.$k
OUT=$PREF.$VAR

# keep a few files
cp ../$DIRRESULT/$FILE.1.1.log.txt $OUT.1.1.log.txt
cp ../$DIRRESULT/$FILE.1.1.assoc.txt.gz $OUT.1.1.assoc.txt.gz
cp ../$DIRRESULT/$FILE.23.1.log.txt $OUT.23.1.log.txt

JOBSUB=$USE.$VAR.sub
JOBRUN=$USE.$VAR.run
JOBOUT=${USE}_${VAR}.out

RUNDIR=/working/qdu/$RUN/$DIRUSE/${PREF}_out
DIR="\/working\/qdu\/$RUN\/$DIRUSE\/${DIRRESULT}"

cat << __EOT__  >$JOBSUB
#PBS -l ncpus=1
#PBS -l mem=1gb
#PBS -l walltime=40:00:00 
cd $RUNDIR || exit 
(
./$JOBRUN > $JOBOUT  2>&1 
) &
wait
__EOT__

sed -e "s/REPPREF/$USE/g" \
    -e "s/REPVAR/$VAR/g" \
    -e "s/REPDIR/$DIR/g" \
    -e "s/REPOUT/$OUT/g" \
    -e "s/REPKLST/$KLST/g" \
    -e "s/REPNUM/$NUM/g" \
    -e "s/REPFILE/$FILE/g"  $REFDIR/gem.$GVER.info.use  >$JOBRUN

chmod +x $JOBRUN
qsub $JOBSUB

# k
done



exit



# 1: run gwas

DIRUSE=fev_${PREF}
mkdir -p $DIRUSE
cd $DIRUSE

DIRRESULT=out_${PREF}
mkdir -p $DIRRESULT

PEDDIR=/working/qdu/$RUN
GRMDIR=/working/qdu/bmsw/grm
PLKDIR=/working/qdu/impute2/imp2/$GVER
FAMDIR=/working/qdu/impute2/imp2/merge/frq_hwe

RUNDIR=/working/qdu/$RUN/$DIRUSE

paste -d" " $FAMDIR/u2.1.fam $PEDDIR/$USE.fam  >fam.chk
awk '$2!=$8 { print }' fam.chk | wc -l

KLST=$GVER.60k.lst

for ((CHR=1;CHR<=23;CHR++))
do
awk -v chr=$CHR '$1==chr {print $2 }' $PLKDIR/$KLST >$CHR.jobs.lst
  for j in `cat $CHR.jobs.lst`
  do

#  fev: fev fvc evr

for ((k=1;k<=3;k++))
do

JOB=${PREF}_${k}_${CHR}_${j}
JOBRUN=$JOB.cb

cat << __EOT__  >$JOBRUN
#PBS -l ncpus=1
#PBS -l mem=1gb
#PBS -l walltime=40:00:00 
cd $RUNDIR
mkdir -p $JOB
cd $JOB
cp $PLKDIR/set.$CHR.$j.bed .
cp $PLKDIR/set.$CHR.$j.bim .
cp $PEDDIR/$USE.fam  set.$CHR.$j.fam
cp $PEDDIR/$USE.cov .

$GEMMA -bfile set.$CHR.$j -k $GRMDIR/grm.bmsw.sXX.txt.gz  -lmm -n $k -miss 0.8 -c $USE.cov -o $PREF.$k.$CHR.$j

cd output
gzip $PREF.$k.$CHR.$j.assoc.txt
mv -f $PREF.$k.$CHR.$j.assoc.txt.gz $RUNDIR/$DIRRESULT/
mv -f $PREF.$k.$CHR.$j.log.txt $RUNDIR/$DIRRESULT/
cd ../..
rm -f -r $JOB
__EOT__

qsub $JOBRUN

#k
  done

#j
 done

rm -f $CHR.jobs.lst

#chr
done

#
## 60K files ready, then run gwas
#
## END 
#################


