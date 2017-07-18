all:
	g++ hawk.cpp cpp/src/*.cpp -lm -lpthread -o hawk
	g++ bonf_fasta.cpp -o bonf_fasta
	g++ kmersearch.cpp -o kmersearch
	g++ kmersummary.cpp -o kmersummary