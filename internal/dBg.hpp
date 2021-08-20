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

				//cout << "out = " << out << endl;

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

				if(has_dollars(prev_kmer)){

					padded_node.push_back(true);
					padded_kmers++;

				}else{

					padded_node.push_back(false);

				}

			}else{

				out += get_edge(kmers[i]);

			}

			if(i < kmers.size()) prev_kmer = kmers[i];

		}

		assert(BWT.length() == OUT.size());

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

		//for each letter, count how many nodes with that incoming letter we have
		C_node[toINT('$')] = 1;
		C_node[toINT('A')] = 0;
		C_node[toINT('C')] = 0;
		C_node[toINT('G')] = 0;
		C_node[toINT('T')] = 0;

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

				assert(last_char(prev_kmer,k) != toINT('$'));
				C_node[last_char(prev_kmer,k)]++;

			}else{

				IN.push_back(false);

			}

			if(i < kmers.size()) prev_kmer = kmers[i];

		}

		//cumulate C

		for(int i=4;i>0;--i) C[i] = C[i-1];
		C[0] = 0;
		for(int i=1;i<5;++i) C[i] += C[i-1];

		for(int i=4;i>0;--i) C_node[i] = C_node[i-1];
		C_node[0] = 0;
		for(int i=1;i<5;++i) C_node[i] += C_node[i-1];


		/*cout << "BWT/OUT/IN/padded:\n" << BWT << endl;
		for(auto b : OUT) cout << uint(b);
		cout << endl;
		for(auto b : IN) cout << uint(b);
		cout << endl;
		for(auto b : padded_node) cout << uint(b);
		cout << endl;

		cout << "C[$] = " << C[toINT('$')] << endl;
		cout << "C[A] = " << C[toINT('A')] << endl;
		cout << "C[C] = " << C[toINT('C')] << endl;
		cout << "C[G] = " << C[toINT('G')] << endl;
		cout << "C[T] = " << C[toINT('T')] << endl;

		cout << "C_node[$] = " << C_node[toINT('$')] << endl;
		cout << "C_node[A] = " << C_node[toINT('A')] << endl;
		cout << "C_node[C] = " << C_node[toINT('C')] << endl;
		cout << "C_node[G] = " << C_node[toINT('G')] << endl;
		cout << "C_node[T] = " << C_node[toINT('T')] << endl;*/

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

	void prune(){

		//reverse edges of the dBg
		//for each node, the starting position of its edges in vector edges
		auto nodes_start = vector<uint64_t>(nr_of_nodes,0);
		vector<uint64_t> edges;

		uint64_t F_pos[5];
		uint64_t curr_node_in[5];

		for(int i=0;i<5;++i){
			F_pos[i] = C[i];
			curr_node_in[i] = C_node[i];
 		}

		uint64_t curr_node_out = 0;
		uint64_t m = 0; //number of edges

		for(uint64_t bwt_pos=0;bwt_pos<BWT.length();++bwt_pos){

			char c = BWT[bwt_pos];

			if(c != '$'){

				//the dBg has an edge: (u,v)
				uint64_t u = curr_node_out;
				uint64_t v = curr_node_in[toINT(c)];

				//in the reverse dBg, node v has an outgoing edge
				nodes_start[v]++;
				m++;

				curr_node_in[toINT(c)] += IN[F_pos[toINT(c)]++];

			}

			curr_node_out += OUT[bwt_pos];

		}

		//cumulate nodes_start
		for(uint64_t i=nodes_start.size()-1;i>0;--i) nodes_start[i] = nodes_start[i-1];
		nodes_start[0] = 0;
		for(int i=1;i<nodes_start.size();++i) nodes_start[i] += nodes_start[i-1];

		edges = vector<uint64_t>(m);

		for(int i=0;i<5;++i){
			F_pos[i] = C[i];
			curr_node_in[i] = C_node[i];
 		}

		curr_node_out = 0;

		for(uint64_t bwt_pos=0;bwt_pos<BWT.length();++bwt_pos){

			char c = BWT[bwt_pos];

			if(c != '$'){

				//the dBg has an edge: (u,v)
				uint64_t u = curr_node_out;
				uint64_t v = curr_node_in[toINT(c)];

				//in the reverse dBg, node v has an outgoing edge towards u
				edges[nodes_start[v]++] = u;

				curr_node_in[toINT(c)] += IN[F_pos[toINT(c)]++];

			}

			curr_node_out += OUT[bwt_pos];

		}

		for(uint64_t i=nodes_start.size()-1;i>0;--i) nodes_start[i] = nodes_start[i-1];
		nodes_start[0] = 0;

		nodes_start.push_back(m);

		/*cout << "graph: " << endl;
		for(auto x : nodes_start) cout << x << " ";
		cout << endl;
		for(auto x : edges) cout << x << " ";
		cout << endl;*/

		//find unnecessary nodes
		auto necessary = vector<bool>(nr_of_nodes);
		for(uint64_t i=0;i<nr_of_nodes;++i) necessary[i] = not padded_node[i];

		for(uint64_t i=0;i<nr_of_nodes;++i){

			if(not padded_node[i]){

				//pick first successor of i
				uint64_t next = nodes_start[i+1]-nodes_start[i] > 0 ? edges[nodes_start[i]] : nr_of_nodes;

				// if this non-padded node has only one successor, which is a padded node,
				// then the successor is necessary.
				if(nodes_start[i+1]-nodes_start[i] == 1 && padded_node[next]){

					//padded nodes must have exactly 1 successor (except node 0 = $$$)
					assert(next == 0 or nodes_start[next+1]-nodes_start[next] == 1);

					while(next != 0 and not necessary[next]){

						necessary[next] = true;
						next = edges[nodes_start[next]];

					}

					if(next == 0) necessary[next] = true;

				}

			}

		}


		uint64_t new_nr_of_nodes = 0;
		for(uint64_t i=0;i<necessary.size();++i) new_nr_of_nodes += necessary[i];

		vector<bool> new_padded_node;
		for(uint64_t i=0;i<necessary.size();++i)
			if(necessary[i])
				new_padded_node.push_back(padded_node[i]);

		vector<bool> new_IN;
		uint64_t new_C[5];
		uint64_t new_C_node[5];

		for(int k = 0;k<5;++k){
			new_C[k]=0;
			new_C_node[k]=0;
		}

		if(necessary[0]){
			new_C[toINT('$')]++;
		}

		uint8_t F_char = 0;//current char on F column


		for(uint64_t i=0;i<nr_of_nodes;++i){

			//incoming letter of node i
			while(F_char < 4 && not (i >= C_node[F_char] and i <  C_node[F_char+1]) ) F_char++;

			if(necessary[i]){

				new_C_node[F_char]++;

				uint64_t in_deg = i == 0;

				assert(i == 0 or nodes_start[i+1]-nodes_start[i] > 0);

				for(uint64_t j = nodes_start[i]; j < nodes_start[i+1]; ++j){

					if(necessary[edges[j]])
						in_deg++;

				}

				assert(in_deg>0);

				for(uint64_t j = 0;j< in_deg-1;++j){

					new_IN.push_back(false);
					new_C[F_char]++;

				}

				new_C[F_char]++;
				new_IN.push_back(true);

			}//else: skip this node

		}

		for(int i=4;i>0;--i) new_C[i] = new_C[i-1];
		new_C[0] = 0;
		for(int i=1;i<5;++i) new_C[i] += new_C[i-1];

		for(int i=4;i>0;--i) new_C_node[i] = new_C_node[i-1];
		new_C_node[0] = 0;
		for(int i=1;i<5;++i) new_C_node[i] += new_C_node[i-1];


		string new_BWT;
		vector<bool> new_OUT;

		for(int i=0;i<5;++i){
			F_pos[i] = C[i];
			curr_node_in[i] = C_node[i];
 		}

		curr_node_out = 0;

		int new_out_deg = 0;

		for(uint64_t bwt_pos=0;bwt_pos<BWT.length();++bwt_pos){

			char c = BWT[bwt_pos];

			if(necessary[curr_node_out]){

				if(c != '$'){

					//the dBg has an edge labeled c: (u,v)
					uint64_t u = curr_node_out;
					uint64_t v = curr_node_in[toINT(c)];

					new_out_deg++;
					new_BWT += c;

				}else{

					new_out_deg++;
					new_BWT += c;

				}

				if(OUT[bwt_pos]){

					assert(new_out_deg > 0);

					for(int k=0;k<new_out_deg-1;++k) new_OUT.push_back(false);
					new_OUT.push_back(true);

				}

			}

			if(c != '$') curr_node_in[toINT(c)] += IN[F_pos[toINT(c)]++];

			if(OUT[bwt_pos]){
				new_out_deg = 0;
				curr_node_out++;
			}

		}

		nr_of_nodes = new_nr_of_nodes;
		padded_kmers = 0;

		padded_node = new_padded_node;
		for(auto b:padded_node) padded_kmers += b;

		BWT = new_BWT;

		nr_of_edges = 0;
		for(auto c:BWT) nr_of_edges += c!='$';

		OUT = new_OUT;
		IN = new_IN;

		for(int k=0;k<5;++k){
			C[k] = new_C[k];
			C_node[k] = new_C_node[k];
		}

		/*
		cout << "BWT/OUT/IN/padded:\n" << BWT << endl;
		for(auto b : OUT) cout << uint(b);
		cout << endl;
		for(auto b : IN) cout << uint(b);
		cout << endl;
		for(auto b : padded_node) cout << uint(b);
		cout << endl;

		cout << "C[$] = " << C[toINT('$')] << endl;
		cout << "C[A] = " << C[toINT('A')] << endl;
		cout << "C[C] = " << C[toINT('C')] << endl;
		cout << "C[G] = " << C[toINT('G')] << endl;
		cout << "C[T] = " << C[toINT('T')] << endl;

		cout << "C_node[$] = " << C_node[toINT('$')] << endl;
		cout << "C_node[A] = " << C_node[toINT('A')] << endl;
		cout << "C_node[C] = " << C_node[toINT('C')] << endl;
		cout << "C_node[G] = " << C_node[toINT('G')] << endl;
		cout << "C_node[T] = " << C_node[toINT('T')] << endl;
		*/

	}

	void minimize(){};

private:

	uint8_t k = 0; //order
	uint64_t nr_of_nodes = 0; //number of nodes (= kmers + dummy-padded nodes)
	uint64_t nr_of_edges = 0; //number of edges
	uint64_t padded_kmers = 0;//number of nodes corresponding to a k-mers padded with $

	//packed Wheeler graph representation
	string BWT;
	vector<bool> OUT; //marks (bit 1) the last outgoing edge of each node in BWT order
	vector<bool> IN;  //marks (bit 1) the last incoming edge of each node in BWT order
	vector<bool> padded_node; // marks padded (dummy) kmers, i.e. kmers containing a $
	uint64_t C[5];//starting positions of each letter in column F. Accessed as C[toINT(char)]

	uint64_t C_node[5]; //C_node[c] = starting positions of nodes with incoming letter c in the list of all nodes.

};



}

#endif
