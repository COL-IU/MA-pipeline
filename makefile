all: create generatescript processsam perlscripts mergeprintbases checkscars mutationsimulator analyzemutations appendbasesaround getpolymorphisms analyzemutationsshared processpileupindel mergepileups analyzeindel selectlines selectlinesandupdatelistfile printconsensus mutationsimulator makelist

clean:
	rm -rf src/*/classes
	rm bin/*.jar bin/*.pl

create:
	mkdir -p bin

generatescript:
#	rm -rf src/GenerateScript/classes
	mkdir -p src/GenerateScript/classes
	javac src/GenerateScript/*.java -d src/GenerateScript/classes
	jar cfe bin/GenerateScript.jar GenerateScript -C src/GenerateScript/classes .

processsam:
#	rm -rf src/ProcessSam/classes
	mkdir -p -p src/ProcessSam/classes
	javac -cp src/ProcessSam src/ProcessSam/MappingStat.java -d src/ProcessSam/classes
	javac -cp src/ProcessSam src/ProcessSam/ProcessAllSAMs.java -d src/ProcessSam/classes
	jar cfe bin/ProcessAllSAMs.jar ProcessAllSAMs -C src/ProcessSam/classes .

perlscripts:
	cp src/printBases.pl bin/
	cp src/printConsensus.pl bin/
	cp src/filterSig.pl bin/
	cp src/filterSig2.pl bin/

selectlinesandupdatelistfile:
	mkdir -p src/SelectLinesAndUpdateListFile/classes
	javac src/SelectLinesAndUpdateListFile/SelectLinesAndUpdateListFile.java -d src/SelectLinesAndUpdateListFile/classes
	jar cfe bin/SelectLinesAndUpdateListFile.jar SelectLinesAndUpdateListFile -C src/SelectLinesAndUpdateListFile/classes .

printconsensus:
	mkdir -p src/PrintConsensus/classes
	javac src/PrintConsensus/PrintConsensus.java -d src/PrintConsensus/classes
	jar cfe bin/PrintConsensus.jar PrintConsensus -C src/PrintConsensus/classes .

mergeprintbases:
#	rm -rf src/MergePrintBases/classes
	mkdir -p src/MergePrintBases/classes
	javac -cp src/MergePrintBases/ src/MergePrintBases/MergePrintBases.java -d src/MergePrintBases/classes
	jar cfe bin/MergePrintBases.jar MergePrintBases -C src/MergePrintBases/classes .

checkscars:
	mkdir -p src/MergePrintBases/classes
	javac -cp src/MergePrintBases/ src/MergePrintBases/CheckScars.java -d src/MergePrintBases/classes
	jar cfe bin/CheckScars.jar CheckScars -C src/MergePrintBases/classes .

mutationsimulator:
	mkdir -p src/AnalyzeMutations/classes
	javac src/AnalyzeMutations/SubstitutionMatrixParser.java
	javac -cp src/AnalyzeMutations src/AnalyzeMutations/MutationSimulatorDriver.java
	mv src/AnalyzeMutations/*.class src/AnalyzeMutations/classes
	jar cfe bin/MutationSimulator.jar MutationSimulatorDriver -C src/AnalyzeMutations/classes .

analyzemutations:
#	rm -rf src/AnalyzeMutations/classes
	mkdir -p src/AnalyzeMutations/classes
	javac src/AnalyzeMutations/SubstitutionMatrixParser.java
	javac -cp src/AnalyzeMutations src/AnalyzeMutations/AnalyzeMutations.java
	mv src/AnalyzeMutations/*.class src/AnalyzeMutations/classes
	jar cfe bin/AnalyzeMutations.jar AnalyzeMutations -C src/AnalyzeMutations/classes .

appendbasesaround:
#	rm -rf src/AppendBasesAround/classes
	mkdir -p src/AppendBasesAround/classes
	javac src/AppendBasesAround/SubstitutionMatrixParser.java
	javac -cp src/AppendBasesAround src/AppendBasesAround/AppendBasesAround.java
	mv src/AppendBasesAround/*.class src/AppendBasesAround/classes
	jar cfe bin/AppendBasesAround.jar AppendBasesAround -C src/AppendBasesAround/classes .

getpolymorphisms:
#	rm -rf src/GetPolymorphisms/classes
	mkdir -p src/GetPolymorphisms/classes
	javac -cp src/GetPolymorphisms src/GetPolymorphisms/GetPolymorphisms.java
	mv src/GetPolymorphisms/*.class src/GetPolymorphisms/classes
	jar cfe bin/GetPolymorphisms.jar GetPolymorphisms -C src/GetPolymorphisms/classes .

analyzemutationsshared:
#	rm -rf src/AnalyzeMutationsShared/classes
	mkdir -p src/AnalyzeMutationsShared/classes
	javac src/AnalyzeMutationsShared/SubstitutionMatrixParser.java
	javac -cp src/AnalyzeMutationsShared src/AnalyzeMutationsShared/AnalyzeMutationsShared.java
	mv src/AnalyzeMutationsShared/*.class src/AnalyzeMutationsShared/classes
	jar cfe bin/AnalyzeMutationsShared.jar AnalyzeMutationsShared -C src/AnalyzeMutationsShared/classes .

processpileupindel:
#	rm -rf src/ProcessPileup_indel/classes
	mkdir -p src/ProcessPileup_indel/classes
	javac -cp src/ProcessPileup_indel src/ProcessPileup_indel/ProcessPileup.java
	mv src/ProcessPileup_indel/*.class src/ProcessPileup_indel/classes
	jar cfe bin/ProcessPileup_indel.jar ProcessPileup -C src/ProcessPileup_indel/classes .

mergepileups:
#	rm -rf src/ProcessPileup_indel/classes
	mkdir -p src/ProcessPileup_indel/classes
	javac -cp src/ProcessPileup_indel src/ProcessPileup_indel/MergePileups.java
	mv src/ProcessPileup_indel/*.class src/ProcessPileup_indel/classes
	jar cfe bin/MergePileups.jar MergePileups -C src/ProcessPileup_indel/classes .

analyzeindel:
#	rm -rf src/ProcessPileup_indel/classes
	mkdir -p src/ProcessPileup_indel/classes
	javac src/ProcessPileup_indel/SubstitutionMatrixParser.java
	javac -cp src/ProcessPileup_indel src/ProcessPileup_indel/AnalyzeIndel.java
	mv src/ProcessPileup_indel/*.class src/ProcessPileup_indel/classes
	jar cfe bin/AnalyzeIndel.jar AnalyzeIndel -C src/ProcessPileup_indel/classes .

selectlines:
#	rm -rf src/SelectLines/classes
	mkdir -p src/SelectLines/classes
	javac src/SelectLines/ExtComposition.java
	javac -cp src/SelectLines src/SelectLines/PBCount.java
	javac -cp src/SelectLines src/SelectLines/SelectLines.java
	mv src/SelectLines/*.class src/SelectLines/classes
	jar cfe bin/SelectLines.jar SelectLines -C src/SelectLines/classes .

makelist:
	mkdir -p src/MakeList/classes
	javac -cp src/MakeList src/MakeList/MakeList.java
	mv src/MakeList/*.class src/MakeList/classes
	jar cfe bin/MakeList.jar MakeList -C src/MakeList/classes .
