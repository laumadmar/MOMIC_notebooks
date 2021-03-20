#!/bin/bash
if [ $2 == "true" ]; then
for i in {1..22}
do
java -jar /mnt/data/GWAS/tools/SnpSift.jar extractFields $1/chr$i.dose.vcf.gz ID REF ALT "GEN[*].DS" | tail -n +2 | gzip -c > $1/chr$i.dose.rsq.0.3.DS.vcf.gz
done
else
for i in {1..22}
do
java -jar /mnt/data/GWAS/tools/SnpSift.jar filter "(R2>0.3)" $1/chr$i.dose.vcf.gz | gzip -c > $1/chr$i.dose.rsq.0.3.vcf.gz
java -jar /mnt/data/GWAS/tools/SnpSift.jar extractFields $1/chr$i.dose.rsq.0.3.vcf.gz ID REF ALT "GEN[*].DS" | tail -n +2 | gzip -c > $1/chr$i.dose.rsq.0.3.DS.vcf.gz
rm $1/chr$i.dose.rsq.0.3.vcf.gz
done
fi
