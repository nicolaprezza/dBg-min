// Copyright (c) 2021, Nicola Prezza.  All rights reserved.
// Use of this source code is governed
// by a MIT license that can be found in the LICENSE file.

/*
 * dBg.hpp
 *
 *  Created on: Aug 17, 2021
 *      Author: nico
 *
 *
 */

#ifndef INCLUDED_DBG
#define INCLUDED_DBG

#include <vector>
#include <algorithm>
#include <stack>
#include <queue>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;

namespace dbg{

enum format_t {fasta, fastq};
typedef pair<uint64_t,uint8_t> edge_t;

/*
 * convert DNA alphabet {$,A,C,G,T} to integers in [0,4] (3 bits per int)
 */
uint8_t toINT(char c){

	switch(c){
		case '$': return 0; break;
		case 'A': case 'a': return 1; break;
		case 'C': case 'c': return 2; break;
		case 'G': case 'g': return 3; break;
		case 'T': case 't': return 4; break;
		default:break;

	}

	return 1;//this includes 'N' characters.
}

//inverse of the above
char toCHAR(uint8_t x){

	switch(x){
		case 0: return '$'; break;
		case 1: return 'A'; break;
		case 2: return 'C'; break;
		case 3: return 'G'; break;
		case 4: return 'T'; break;
		default:break;

	}

	return '$';

}

//format example: if k=3 and we have a kmer ACG followed by T, then an edge is represented as an integer rev(ACG)T = GCAT

/*
 * input: edge (XYZ,W) stored as integer of 128 bits (see "format example" above), character c stored in 3 bits (see function toINT), and order k
 * output: edge (YZW,c)
 */
__uint128_t edge(__uint128_t kmer, uint8_t c, uint8_t k){

	return (((kmer >> 3) | ((kmer & __uint128_t(7))<<(3*k))) & (~__uint128_t(7))) | c;

}

/*
 * input: edge (XYZ,W) stored as integer of 128 bits (see "format example" above)
 * output: character Z
 */
uint8_t last_char(__uint128_t kmer, uint8_t k){

	return (kmer >> (3*k)) & __uint128_t(7);

}

/*
 * input: edge (XYZ,W) stored as integer of 128 bits (see "format example" above)
 * output: character W on the edge
 */
char get_edge(__uint128_t kmer){

	return toCHAR(kmer & __uint128_t(7));

}

string kmer_to_str_(__uint128_t kmer, int k){

	string km;

	for(int i=0;i<k;++i){

		km += toCHAR(kmer & __uint128_t(7));
		kmer = kmer >> 3;

	}

	return km;

}

string kmer_to_str(__uint128_t kmer, int k){

	char edge = toCHAR(kmer & __uint128_t(7));

	string km;

	kmer = kmer >> 3;

	for(int i=0;i<k;++i){

		km += toCHAR(kmer & __uint128_t(7));
		kmer = kmer >> 3;

	}

	km += "|";
	km += edge;

	return km;

}

//has $ symbols in the kmer part?
bool has_dollars(__uint128_t kmer){

	return ((kmer >> 3) &  __uint128_t(7)) == 0;

}



class dBg{

public:

	dBg();

	dBg( const string & filename,
			format_t format,
			int nlines = 0,
			uint8_t k = 28,
			bool do_not_optimize = false,
			bool verbose = true) : k(k){

		assert(k>0 and k<=41);

		if(verbose)
			cout << "Computing how much memory I need to allocate ..." << endl;

		uint64_t pre_allocation = 0;
		uint64_t tot_bases = 0;

		{
			//count how many kmers we will generate (one per base)

			ifstream file(filename);
			int read_lines=0;

			string str;
			while (std::getline(file, str) and (nlines == 0 or read_lines < nlines)) { //getline reads header

				getline(file, str);//getline reads DNA

				pre_allocation += str.length()+1;
				tot_bases += str.length();

				if(format == fastq){
					getline(file, str);//getline reads +
					getline(file, str);//getline reads quality
				}

				read_lines++;

				if(read_lines%1000000==0 and verbose)
					cout << "read " << read_lines << " sequences" << endl;

			}

		}

		ifstream file(filename);

		if(verbose){
			cout << "Number of bases: " << tot_bases << endl;
			cout << "Trying to allocate " << pre_allocation*16 << " Bytes ... " << endl;
		}

		int read_lines=0;

		/*
		 *  vector storing all (k+1)-mers (i.e. edges) using 3 bits per char
		 *  format example: if k=3 and we have a kmer ACG followed by T, then we store an integer rev(ACG)T = GCAT
		 *
		 */
		vector<__uint128_t> kmers;
		kmers.reserve(pre_allocation);

		if(verbose)
			cout << "Extracting k-mers from dataset ..." << endl;

		string str;
		while (std::getline(file, str) and (nlines == 0 or read_lines < nlines)) { //getline reads header

			getline(file, str);//getline reads DNA

			// start processing DNA fragment

			//first kmer: ($$$,C), where C is the first letter of str

			uint8_t first_char = toINT(str[0]);

			//if(first_char==0) cout << str<<endl;
			assert(first_char!=0);
			__uint128_t kmer = first_char;

			kmers.push_back(kmer);

			//push the other kmers
			for(int i=1;i<str.length();++i){

				kmer = edge(kmer, toINT(str[i]),k);
				kmers.push_back(kmer);

			}

			//last kmer: (ACG,$), where ACG was the last kmer in the DNA fragment
			kmer = edge(kmer, toINT('$'),k);
			kmers.push_back(kmer);

			if(format == fastq){
				getline(file, str);//getline reads +
				getline(file, str);//getline reads quality
			}

			read_lines++;

			if(read_lines%1000000==0 and verbose)
				cout << "read " << read_lines << " sequences" << endl;

		}

		if(verbose)
			cout << "Sorting and de-duplicating k-mers ..." << endl;

		sort(kmers.begin(),kmers.end());
		auto it = unique(kmers.begin(), kmers.end());
		auto new_size = distance(kmers.begin(), it);
		kmers.resize(new_size);
		kmers.shrink_to_fit();

		nr_of_edges = 0;

		if(verbose)
			cout << "Packing dBg graph into a BWT ..." << endl;

		//previous kmer read from kmers.
		__uint128_t prev_kmer = kmers[0];
		string out;
		out += get_edge(kmers[0]);

		for(uint64_t i = 1; i<=kmers.size();++i){

			//cout << "out: " << out << endl;

			//if kmer changes or i goes beyond the end of kmers
			if(i == kmers.size() || (kmers[i]>>3) != (prev_kmer>>3)){

				assert(i == kmers.size() or out.length() > 0);

				cout << "out = " << out << endl;

				if(out.length() == 1){

					BWT += out[0];
					OUT.push_back(true);

				}else{

					for(int j=0;j<out.size();++j){

						if(out[j]!='$'){

							BWT += out[j];
							OUT.push_back(false);

						}

					}

					OUT.back() = true;

				}

				for(auto c:out) nr_of_edges += c!='$';


				nr_of_nodes++;

				out = "";
				if(i < kmers.size()) out += get_edge(kmers[i]);

				if(has_dollars(prev_kmer)) padded_kmers++;

			}else{

				out += get_edge(kmers[i]);

			}

			if(i < kmers.size()) prev_kmer = kmers[i];

		}

		if(verbose)
			cout << "Computing incoming edges ..." << endl;

		//compute successor of every kmer: edge (ACG,T) becomes (CGT,$). Ignore edges labeled with $: (ACG,$)

		new_size = 0;
		for(uint64_t i=0;i<kmers.size();++i){

			if(get_edge(kmers[i]) != '$'){

				kmers[new_size++] = edge(kmers[i],0,k);

			}

		}

		sort(kmers.begin(),kmers.begin() + new_size);
		kmers.resize(new_size);
		kmers.shrink_to_fit();

		//for each letter, count how many incoming edges with that letter we have
		C[toINT('$')] = 1; // kmer $$$ (with incoming dummy label $) is not in kmers but we must take it into account
		C[toINT('A')] = 0;
		C[toINT('C')] = 0;
		C[toINT('G')] = 0;
		C[toINT('T')] = 0;

		assert(last_char(kmers[0],k) != toINT('$'));

		C[last_char(kmers[0],k)]++;
		IN.push_back(true); // kmer $$$ is not in kmers but we must take it into account

		prev_kmer = kmers[0];

		for(uint64_t i = 1; i<=kmers.size();++i){

			if(i<kmers.size()){

				assert(last_char(kmers[i],k) != toINT('$'));
				C[last_char(kmers[i],k)]++;

			}

			//if kmer changes or i goes beyond the end of kmers
			if(i == kmers.size() || (kmers[i]>>3) != (prev_kmer>>3)){

				IN.push_back(true);

			}else{

				IN.push_back(false);

			}

			if(i < kmers.size()) prev_kmer = kmers[i];

		}

		//cumulate C

		for(int i=4;i>0;--i) C[i] = C[i-1];
		C[0] = 0;
		for(int i=1;i<5;++i) C[i] += C[i-1];



		cout << "BWT/OUT/IN:\n" << BWT << endl;
		for(auto b : OUT) cout << uint(b);
		cout << endl;
		for(auto b : IN) cout << uint(b);
		cout << endl;

		cout << "C[$] = " << C[toINT('$')] << endl;
		cout << "C[A] = " << C[toINT('A')] << endl;
		cout << "C[C] = " << C[toINT('C')] << endl;
		cout << "C[G] = " << C[toINT('G')] << endl;
		cout << "C[T] = " << C[toINT('T')] << endl;

	}

	uint64_t number_of_nodes(){

		return nr_of_nodes;

	}

	/*
	 * number of nodes (distinct kmers) in the de Bruijn graph. Note: also padded kmers are counted.
	 */
	uint64_t number_of_padded_kmers(){

		return padded_kmers;
	}

	uint64_t number_of_distinct_kmers(){

		return nr_of_nodes - padded_kmers;

	}

	/*
	 * number of edges in the de Bruijn graph.
	 */
	uint64_t number_of_edges(){

		return nr_of_edges;

	}

	void prune(){}

	void minimize(){};

private:

	uint8_t k = 0; //order
	uint64_t nr_of_nodes = 0; //number of nodes
	uint64_t nr_of_edges = 0; //number of edges
	uint64_t padded_kmers = 0;//number of nodes corresponding to a k-mers padded with $

	//packed Wheeler graph representation
	string BWT;
	vector<bool> OUT; //marks (bit 1) the last outgoing edge of each node in BWT order
	vector<bool> IN;  //marks (bit 1) the last incoming edge of each node in BWT order
	uint64_t C[5];//starting positions of each letter in column F. Accessed as C[toINT(char)]

};



}

#endif
