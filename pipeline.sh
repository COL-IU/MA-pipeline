#!/bin/bash
. $1
echo ${MAHOME}

if [[ "$2" == LT* ]]
then
    echo "==================== MAKE LIST ========================"
    java -jar ${MAHOME}/bin/MakeList.jar ${DATADIR} ${STRAIN_NAME} ${JOB} ${MAHOME}/util
    cd ${DATADIR};cat ${LIST} | xargs -P${NUM_CPU_OUTER} -I{} sh -c "cd {};mv {}_1.fq.gz {}_1.fq.orig.gz;mv {}_2.fq.gz {}_2.fq.orig.gz;java -jar ${TRIMJAR} PE -threads 1 -phred33 -trimlog {}.trimlog {}_1.fq.orig.gz {}_2.fq.orig.gz {}_1.fq.gz {}_1.fq.singleton.gz {}_2.fq.gz {}_2.fq.singleton.gz ILLUMINACLIP:${TRIMHOME}adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:70 &> {}.trimstat"
    cd ${MAHOME}
elif [[ "$2" == L* ]]
then
    java -jar ${MAHOME}/bin/MakeList.jar ${DATADIR} ${STRAIN_NAME} ${JOB} ${MAHOME}/util
fi

if [[ "$2" == *1* ]]
then
    	mkdir ${OUTDIR}
	mkdir ${OUTDIR}/scripts
	#[9/14/13 by heewlee] ADDED this part to accomodate different read lenths across the lines in a single MA set
	cat ${LIST} | xargs -I{} sh -c "zcat ${DATADIR}/{}/{}_1.fq.gz | head -n2 | tail -n1 " > ${OUTDIR}/fr.tmp
	paste ${LIST} ${OUTDIR}/fr.tmp > ${OUTDIR}/firstReads.list
	rm ${OUTDIR}/fr.tmp
	#[9/14/13 by heewlee] END OF ADDITION
	cat ${LIST} | xargs -I{} -P10 sh -c "java -jar ${MAHOME}/bin/GenerateScript.jar -B {} -R ${MAHOME}/util -r ${REF} -D ${DATADIR} -o ${OUTDIR} -l ${READ_LEN} -N ${OUTDIR}/scripts/${JOB}{} -n ${NUM_CPU_BWA} -m ${MAX_NUM_MISMATCH} -S F -W F -U T -K F -FS .fq.gz -M ${MAHOME} -H ${HEADERLIST} -MQ ${MIN_MAPQUAL}"

	chmod a+x ${OUTDIR}/scripts/*

	cat ${LIST} | xargs -I{} mkdir ${OUTDIR}/{}

	ls ${OUTDIR}/scripts/ | xargs -I{} -P${NUM_CPU_OUTER} sh -c "${OUTDIR}/scripts/{}"
fi

if [[ "$2" == *2* ]]
then
	#[9/10/13] ADDED SelectLinesAndUPdateListFile to accmodate different sequencing covereage vals across lines. 
	cd ${OUTDIR};java -jar ${MAHOME}/bin/SelectLinesAndUpdateListFile.jar ${LIST} ${OUTDIR} ${MAINCHROMOSOME}
fi

if [[ "$2" == *3* ]]
then
	cd ${OUTDIR};cat ${LIST} | xargs -I{} -P20 sh -c "cd {};perl ${MAHOME}/bin/printBases.pl {}_${REFBASE}.qc.PE.pileup > {}.printBases"
fi

if [[ "$2" == *F ]]
then
	cd ${OUTDIR};cat ${LIST} | xargs -I{} -P20 sh -c "cd {};python ${MAHOME}/util/alleleFreq.py {}.printBases {};${MAHOME}/util/typifyAlleles.py {}.printBases.allele > {}.types;${MAHOME}/util/alleleToPutations.py {}.printBases.allele > {}.allele.putations;echo -e '1\t'{} >> dummy.wNum;java -jar ${MAHOME}/bin/AnalyzeMutations.jar {}.allele.putations ${MAHOME}/util/standard_codon.txt ${MAHOME}/util/${REFBASE}.fna ${MAHOME}/util/BLOSUM62.txt 0 Y dummy.wNum N {}_snp ${PTTLIST};${MAHOME}/util/augmentAnnotation.py {}.printBases.allele {}.allele.putations.${MAINCHROMOSOME}.detail.wBases > {}.allele.putations.${MAINCHROMOSOME}.detail.wBases.wAlleles"
fi

if [[ "$2" == *I ]]
then
	#cat ${LIST} | xargs -I{} -P${NUM_CPU_OUTER} sh -c "java -jar ${MAHOME}/bin/ProcessPileup_indel.jar ${OUTDIR}/{}/{}_${REFBASE}.qc.PE.pileup > ${OUTDIR}/{}/{}.indels"
  cat ${LIST} | xargs -I{} -P${NUM_CPU_OUTER} sh -c "python ${MAHOME}/util/indelAlleleFreq.py ${OUTDIR}/{}/{}.indels"
	cat ${LIST} | xargs -I{} -P${NUM_CPU_OUTER} sh -c "python ${MAHOME}/util/prepIndelsForAnalysis.py ${OUTDIR}/{}/{}.indels > ${OUTDIR}/{}/{}.prepped.indels; java -jar ${MAHOME}/bin/AnalyzeIndel.jar ${OUTDIR}/{}/{}.prepped.indels ${OUTDIR}/{}/dummy.wNum ${MAHOME}/util/${REFBASE}.fna -O:${OUTDIR}/{}/{} ${PTTLIST}"
	cat ${LIST} | xargs -I{} -P${NUM_CPU_OUTER} sh -c "python ${MAHOME}/util/mergeIndelAllele.py ${OUTDIR}/{}/{}.indels.allele ${OUTDIR}/{}/{}.pindels > ${OUTDIR}/{}/{}.merged.pindels"
fi

if [[ "$2" == *4* ]]
then
  java -jar ${MAHOME}/bin/MergePrintBases.jar ${HEADERLIST} ${LIST} ${JOB}.alignPos
fi

if [[ "$2" == *5* ]]
then
	java -jar ${MAHOME}/bin/CheckScars.jar ${OUTDIR}/${JOB}.alignPos ${MAHOME}/util/scars > ${OUTDIR}/${JOB}.scars
fi

if [[ "$2" == *6* ]]
then
	#[9/10/13] Ported the perl code to java to speed up and support deviating coverage vals across lines in a MA set
	java -jar ${MAHOME}/bin/PrintConsensus.jar ${OUTDIR}/${JOB}.alignPos ${OUTDIR}/${JOB}.consensus ${LIST} ${OUTDIR} ${MAINCHROMOSOME}
#perl ${MAHOME}/bin/printConsensus.pl ${OUTDIR}/${JOB}.alignPos ${OUTDIR}/${JOB}.consensus 3

	perl ${MAHOME}/bin/filterSig.pl ${OUTDIR}/${JOB}.consensus > ${OUTDIR}/${JOB}.putations

	perl ${MAHOME}/bin/filterSig2.pl ${OUTDIR}/${JOB}.consensus > ${OUTDIR}/${JOB}.shared.putations

	java -jar ${MAHOME}/bin/AnalyzeMutations.jar ${OUTDIR}/${JOB}.putations ${MAHOME}/util/standard_codon.txt ${MAHOME}/util/${REFBASE}.fna  ${MAHOME}/util/BLOSUM62.txt 0 Y ${NUMLIST} N ${OUTDIR}/${JOB}_snp ${PTTLIST}

	#java -jar ${MAHOME}/bin/AppendBasesAround.jar ${MAHOME}/util/${REFBASE}.fna ${OUTDIR}/${JOB}.putations.detail > ${OUTDIR}/${JOB}.putations.detail.wBases

	cp ${OUTDIR}/${JOB}.putations ${OUTDIR}/${JOB}.fake

	java -jar ${MAHOME}/bin/GetPolymorphisms.jar ${OUTDIR}/${JOB}.consensus ${OUTDIR}/${JOB}.fake > ${OUTDIR}/${JOB}.polymorphisms

	java -jar ${MAHOME}/bin/AnalyzeMutations.jar ${OUTDIR}/${JOB}.fake ${MAHOME}/util/standard_codon.txt ${MAHOME}/util/${REFBASE}.fna ${MAHOME}/util/BLOSUM62.txt 0 Y ${MAHOME}/util/poly_num2sam.txt N ${OUTDIR}/${JOB}_poly ${PTTLIST}

	java -jar ${MAHOME}/bin/AnalyzeMutationsShared.jar ${OUTDIR}/${JOB}.shared.putations ${MAHOME}/util/standard_codon.txt ${MAHOME}/util/${REFBASE}.fna ${MAHOME}/util/BLOSUM62.txt 0 Y ${NUMLIST} N ${OUTDIR}/${JOB}_shared_snp ${PTTLIST}
	
	#java -jar ${MAHOME}/bin/AppendBasesAround.jar ${MAHOME}/util/${REFBASE}.fna ${OUTDIR}/${JOB}.shared.putations.detail > ${OUTDIR}/${JOB}.shared.putations.detail.wBases
fi

if [[ "$2" == *7* ]]
then
	cat ${LIST} | xargs -I{} -P${NUM_CPU_OUTER} sh -c "java -jar ${MAHOME}/bin/ProcessPileup_indel.jar ${OUTDIR}/{}/{}_${REFBASE}.qc.PE.pileup > ${OUTDIR}/{}/{}.indels"

	cat ${LIST} | xargs -I{} echo -n "${OUTDIR}/{}/{}.indels " > ${OUTDIR}/${JOB}.args
	
	echo ${OUTDIR}/${JOB}.indels >> ${OUTDIR}/${JOB}.args

	cat ${OUTDIR}/${JOB}.args | xargs java -jar ${MAHOME}/bin/MergePileups.jar ${HEADERLIST}

	java -jar ${MAHOME}/bin/AnalyzeIndel.jar ${OUTDIR}/${JOB}.indels ${NUMLIST} ${MAHOME}/util/${REFBASE}.fna -O:${OUTDIR}/${JOB} ${PTTLIST}
fi

if [[ "$2" == *8* ]]
then
	cat ${LIST} | xargs -I{} -P${NUM_CPU_OUTER} sh -c "bzip2 -f -5 ${OUTDIR}/{}/{}_${REFBASE}.qc.PE.pileup ${OUTDIR}/{}/{}.printBases"
	
	bzip2 -f ${OUTDIR}/${JOB}.alignPos ${OUTDIR}/${JOB}.consensus 
fi

