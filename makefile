SHELL := /bin/bash
CPP_FLAG = -lpthread
CPP_FOLDER=cpp/src
CPP_FILE=${CPP_FOLDER}/*.cpp
OBJECT_FILE=${CPP_FOLDER}/*.o


all: hawk.out  preProcess.out log_reg_case.out log_reg_control.out bonf_fasta.out\
     kmersearch.out kmersummary.out convertToFasta_bf_correction.out convertToFasta_bh_correction.out\
	 kmerStats.out

hawk.out: hawk.o ${OBJECT_FILE}
	g++ $^ ${CPP_FLAG} -o $@

preProcess.out: preProcess.o
	g++ $^ -o $@

log_reg_case.out: log_reg_case.o ${OBJECT_FILE}
	g++ $^  ${CPP_FLAG} -o $@

log_reg_control.out: log_reg_control.o ${OBJECT_FILE}
	g++ $^ ${CPP_FLAG} -o $@

bonf_fasta.out: bonf_fasta.o
	g++ $^ -o $@

kmersearch.out: kmersearch.o
	g++ $^ -o $@

kmersummary.out: kmersummary.o
	g++ $^ -o $@

convertToFasta_bf_correction.out: convertToFasta_bf_correction.o
	g++ $^ -o $@

convertToFasta_bh_correction.out: convertToFasta_bh_correction.o
	g++ $^ -o $@

kmerStats.out: kmerStats.o
	g++ $^ -o $@

hawk.o: hawk.cpp
	g++ $^ -c -o $@

preProcess.o: preProcess.cpp
	g++ $^ -c -o $@

log_reg_case.o: log_reg_case.cpp
	g++ $^ -c -o $@

log_reg_control.o: log_reg_control.cpp
	g++ $^  -c -o $@

bonf_fasta.o: bonf_fasta.cpp
	g++ $^ -c -o $@

kmersearch.o: kmersearch.cpp
	g++ $^ -c -o $@

kmersummary.o: kmersummary.cpp
	g++ $^ -c -o $@

convertToFasta_bf_correction.o: convertToFasta_bf_correction.cpp
	g++ $^ -c -o $@

convertToFasta_bh_correction.o: convertToFasta_bh_correction.cpp
	g++ $^ -c -o $@

kmerStats.o: kmerStats.cpp
	g++ $^ -c -o $@

${OBJECT_FILE}: ${CPP_FILE}
	for f in `ls ${CPP_FILE}`;do echo $$f;g++ $$f -c -o $$f.o; done

clean:
	rm *.o
	rm ${CPP_FOLDER}/*.o

#all:
#	g++ hawk.cpp cpp/src/*.cpp -lm -lpthread -o hawk
#	g++ bonf_fasta.cpp -o bonf_fasta
#	g++ kmersearch.cpp -o kmersearch
#	g++ kmersummary.cpp -o kmersummary
#	g++ preProcess.cpp -o preProcess
#	g++ convertToFasta.cpp -o convertToFasta
