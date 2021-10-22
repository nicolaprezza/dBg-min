// Copyright (c) 2021, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

#include <cassert>
#include <chrono>

#include "internal/dBg.hpp"

using namespace std;
using namespace dbg;

uint64_t nlines = 0;

format_t format = fastq;

bool D = false;
bool pause_ = false;
bool prune = false;

void help(){
	cout << "dBg-min: builds the de Bruijn graph and compute the minimum equivalent Wheeler DFA." << endl << endl;
	cout << "Usage: dBg-min [options] <input> <k>" << endl;
	cout << "   Options:"<<endl;
	cout << "   -l <nlines>         Use only the first nlines sequences from the input file to build the graph. If set to 0, use all lines. Default: 0."<<endl;
	cout << "   -a                  The input file is fasta. If not specified, it is assumed that the input file is fastq."<<endl;
	cout << "   -p                  Optimize space by pruning unnecessary dummy nodes (relevant if number of sequences is large)."<<endl;
	cout << "   <input>             Input fasta/fastq file (see option -a). Mandatory."<<endl;
	cout << "   <k>                 Order of the de Bruijn graph in [1,28]. Mandatory."<<endl;
	exit(0);
}

void parse_args(char** argv, int argc, int &ptr){

	assert(ptr<argc);

	string s(argv[ptr]);
	ptr++;

	if(s.compare("-l")==0){

		if(ptr>=argc-2){
			cout << "Error: missing parameter after -l option." << endl;
			help();
		}

		nlines = atoi(argv[ptr++]);

	}else if(s.compare("-a")==0){

		format = fasta;

	}else if(s.compare("-p")==0){

		prune = true;

	}else{
		cout << "Error: unrecognized '" << s << "' option." << endl;
		help();
	}

}

int main(int argc, char** argv){

	//parse options

	int ptr = 1;

	if(argc<3) help();

	while(ptr<argc-2)
		parse_args(argv, argc, ptr);

	auto input_file = string(argv[ptr++]);
	uint8_t k = atoi(argv[ptr]);

	if(k>41 or k==0){
		cout << "Error: k must be in [1,41]" << endl;
		help();
	}

	cout << "Building de Bruijn graph of input file " << input_file << endl;
	cout << "called as: dBg-min " << (format==fasta?"-a ":"") << (prune?"-p ":"") << "-l " << nlines << " " << input_file << " " << int(k) << endl;

	auto t1 = std::chrono::high_resolution_clock::now();

	dBg G(input_file, format, nlines, k, true);

	auto t2 = std::chrono::high_resolution_clock::now();

	uint64_t elapsed = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
	cout << "\nDone. Elapsed time (hh:mm:ss): " << elapsed/3600 << ":" << (elapsed%3600)/60 << ":" << (elapsed%3600)%60 << endl;

	cout << "Number of kmers " << G.number_of_distinct_kmers() << endl;
	cout << "Number of dummy nodes " << G.number_of_padded_kmers() << endl;
	cout << "Total number of nodes (kmers + dummy nodes) " << G.number_of_nodes() << endl;
	cout << "Number of edges " << G.number_of_edges() << endl;

	auto t3 = std::chrono::high_resolution_clock::now();

	if(prune){

		cout << "\nRemoving unnecessary padded nodes ... " << endl;
		G.prune();
		t3 = std::chrono::high_resolution_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::seconds>(t3 - t2).count();
		cout << "\nDone. Elapsed time (hh:mm:ss): " << elapsed/3600 << ":" << (elapsed%3600)/60 << ":" << (elapsed%3600)%60 << endl;

		cout << "Number of kmers " << G.number_of_distinct_kmers() << endl;
		cout << "Number of dummy nodes " << G.number_of_padded_kmers() << endl;
		cout << "Total number of nodes (kmers + dummy nodes) " << G.number_of_nodes() << endl;
		cout << "Number of edges " << G.number_of_edges() << endl;

	}

	cout << "\nMinimizing dBg ... " << endl;
	G.minimize();
	auto t4 = std::chrono::high_resolution_clock::now();
	elapsed = std::chrono::duration_cast<std::chrono::seconds>(t4 - t3).count();
	cout << "\nDone. Elapsed time (hh:mm:ss): " << elapsed/3600 << ":" << (elapsed%3600)/60 << ":" << (elapsed%3600)%60 << endl;

	//cout << "Number of kmers " << G.number_of_distinct_kmers() << endl;
	//cout << "Number of dummy nodes " << G.number_of_padded_kmers() << endl;
	//cout << "Total number of nodes (kmers + dummy nodes) " << G.number_of_nodes() << endl;
	//cout << "Number of edges " << G.number_of_edges() << endl;


}
