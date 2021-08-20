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

using namespace std;

namespace dbg{

typedef pair<uint64_t,uint8_t> edge_t;

enum format_t {fasta, fastq};

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

/*
 * encode small integers as small positive integers (negative -> odd, positive -> even)
 * output: integer >= 0
 */
uint64_t int_to_positive(int w){

	return w<0?2*(-w)-1:2*w;

}

//reverse of the above
int positive_to_int(uint64_t x){

	return x%2 ? -((int(x)+1)/2) : x/2;

}

//number of bits of x>=1
uint64_t n_bits(uint64_t x){

	assert(x>=1);
	return 64 - __builtin_clzll(x);

}

//gamma-encoding length of the integer w >= 1
uint64_t cost_gamma(uint64_t w){

	assert(w>=1);
	return 2*n_bits(w)-1;

}

//gamma-encoding length of vector of integers >= 0
uint64_t cost_gamma(vector<uint64_t>& V){

	uint64_t cost = 0;

	for(auto w : V) cost += cost_gamma(w+1);

	return cost;

}


//delta-encoding length of the integer w >= 1
uint64_t cost_delta(uint64_t w){

	assert(w>=1);
	return (n_bits(w)-1) + cost_gamma(n_bits(w));

}

//cost function (bits to represent integer). Input: positive integer >= 0
uint64_t cost_of_int(uint64_t w){

	assert(w>=0);
	w += 1; //because gamma and delta can only encode integers >= 1
	return cost_gamma(w);
	//return cost_delta(w);

}


//the cost of a (possibly negative) weight, i.e. the number of bits used to represent it
uint64_t cost_of_weight(int w){

	//return cost_of_int(abs(w));
	return cost_of_int(int_to_positive(w));

}

/*
 * input: edge (XYZ,W) stored as integer of 128 bits (see "format example" below), character c stored in 3 bits (see function toINT), and order k
 * output: edge (YZW,c)
 */
__uint128_t edge(__uint128_t kmer, uint8_t c, uint8_t k){

	return (((kmer >> 3) | ((kmer & __uint128_t(7))<<(3*k))) & (~__uint128_t(7))) | c;

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

	dBg(){}

	/*
	 * constructor from fastq/fasta file
	 *
	 * filename: input file path
	 * format: either fasta or fastq
	 * nlines: if 0, process all DNA fragments. Otherwise, only the first nlines
	 * k: de Bruijn graph order (limited to 41)
	 * do_not_optimize: turn off space optimization (does not prune dBg)
	 */
	cw_dBg( const string & filename,
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
			cout << "Sorting k-mers ..." << endl;

		sort(kmers.begin(),kmers.end());





















		if(verbose)
			cout << "Computing in/out-degrees, edge labels, and weights ..." << endl;

		//previous kmer read from kmers.
		__uint128_t prev_kmer = kmers[0];

		start_positions_out_.push_back(0);
		out_labels_.push_back(toCHAR(kmers[0] & __uint128_t(7)));
		char c = out_labels_[out_labels_.size()-1];
		assert(c=='$' or c=='A' or c=='C' or c=='G' or c=='T');

		uint32_t count = 1;

		for(uint64_t i = 1; i<kmers.size();++i){

			//if kmer changes
			if((kmers[i]>>3) != (prev_kmer>>3)){

				//we set to 0 the counters of kmers that contain $
				count = has_dollars(prev_kmer)?0:count;
				weights_.push_back(count); //I'm pushing the weight of the previous kmer

				if(not has_dollars(prev_kmer)){
					MAX_WEIGHT = count>MAX_WEIGHT?count:MAX_WEIGHT;
					MEAN_WEIGHT += count;
				}else{
					padded_kmers++;
				}

				count = 1; //start counting weight of this new kmer

				start_positions_out_.push_back(out_labels_.size());
				out_labels_.push_back(toCHAR(kmers[i] & __uint128_t(7)));//append to BWT first outgoing edge of this new kmer
				char c = out_labels_[out_labels_.size()-1];
				assert(c=='$' or c=='A' or c=='C' or c=='G' or c=='T');

			}else{//same kmer

				count++;

				uint8_t curr_char = kmers[i] & __uint128_t(7);
				uint8_t prev_char = prev_kmer & __uint128_t(7);

				//if char of outgoing edge has changed
				if( curr_char != prev_char ){

					out_labels_.push_back(toCHAR(curr_char));
					char c = out_labels_[out_labels_.size()-1];
					assert(c=='$' or c=='A' or c=='C' or c=='G' or c=='T');

				}

			}

			prev_kmer = kmers[i];

		}

		if(not has_dollars(prev_kmer)){
			MAX_WEIGHT = count>MAX_WEIGHT?count:MAX_WEIGHT;
			MEAN_WEIGHT += count;
		}else{
			padded_kmers++;
		}

		weights_.push_back(count);//push back weight of the last kmer

		nr_of_nodes = start_positions_out_.size();

		OUT_ = vector<uint64_t>(out_labels_.size());

		//delete char labeling outgoing edge from each (k+1)-mer, obtaining the k-mers
		for(auto & k : kmers)
			k = k>>3;

		//if(verbose)
			//cout << "Resizing k-mers ..." << endl;

		//compact kmers by removing duplicates
		auto it = unique(kmers.begin(), kmers.end());
		auto new_size = distance(kmers.begin(), it);
		kmers.resize(new_size);
		kmers.shrink_to_fit();

		assert(kmers.size() == nr_of_nodes);

		start_positions_in_ = vector<uint64_t>(nr_of_nodes,0);

		for(uint64_t i=0;i<nr_of_nodes;++i){

			uint64_t out_deg = out_degree_(i);

			for(uint8_t off=0;off<out_deg;++off){

				//outgoing label
				assert(start_positions_out_[i]+off<out_labels_.size());
				char c = out_labels_[start_positions_out_[i]+off];

				assert(c=='$' or c=='A' or c=='C' or c=='G' or c=='T');

				if(c!='$'){

					__uint128_t succ_kmer = (kmers[i]>>3) | (__uint128_t(toINT(c))<<(3*(k-1)));

					//cout << "----" << kmer_to_str_(kmers[i],k) << " " << c << " " << kmer_to_str_(succ_kmer,k)  << endl;

					auto it = lower_bound(kmers.begin(), kmers.end(), succ_kmer);

					assert(it != kmers.end());//destination kmer must be present

					uint64_t successor = distance(kmers.begin(), it);

					//cout << kmer_to_str_(succ_kmer,k) << " " << kmer_to_str_(kmers[successor],k) << endl;
					//cout << successor << " / " << nr_of_nodes << endl;

					assert( kmers[successor]==succ_kmer );

					OUT_[start_positions_out_[i]+off] = successor;

					//count how many edges enter in successor
					start_positions_in_[successor]++;

					//if(start_positions_in_[successor]>5)
						//cout << "hmmmm: " << kmer_to_str_(succ_kmer,k) << endl;

				}

			}

		}

		kmers.clear();
		kmers.shrink_to_fit();

		assert(start_positions_in_[0]==0);

		//add one incoming edge (the only one labeled $) to the source
		start_positions_in_[0]=1;

		for(auto x : start_positions_in_){
			assert(x>0); //every node must have >0 and <=5 incoming edges
			assert(x<=5);
		}

		/*
		 *     1  3  2  1  2 --> cumulate -->
		 *     1  4  6  7  9 --> right-shift -->
		 *     0  1  4  6  7
		 */

		//cumulate in start_positions_in_
		for(uint64_t i=1;i<nr_of_nodes;++i){

			start_positions_in_[i] += start_positions_in_[i-1];

		}

		//right-shift start_positions_in_

		for(uint64_t i=nr_of_nodes-1; i>0; --i){

			start_positions_in_[i] = start_positions_in_[i-1];

		}

		//the $ of the source starts at position 0 in IN
		start_positions_in_[0]=0;

		C = vector<uint64_t>(5);

		uint64_t F_length = 1;

		for(auto c : out_labels_){

			C[toINT(c)]++;
			F_length += c!='$';

		}

		C[0] = 1;//there is one $ in the F column (it's the only incoming edge of the root kmer $$$)
		C[1] += C[0];
		C[2] += C[1];
		C[3] += C[2];
		C[4] += C[3];

		C[4] = C[3];
		C[3] = C[2];
		C[2] = C[1];
		C[1] = C[0];
		C[0] = 0; //$ starts at position 0 in the F column

		IN_ = vector<uint64_t>(F_length,nr_of_nodes);

		assert(start_positions_in_[nr_of_nodes-1] < IN_.size());

		//fill IN_

		for(uint64_t i=0;i<nr_of_nodes;++i){

			uint64_t out_deg = out_degree_(i);

			for(uint8_t off=0;off<out_deg;++off){

				//outgoing label
				char c = out_labels_[start_positions_out_[i]+off];
				assert(c=='$' or c=='A' or c=='C' or c=='G' or c=='T');

				if(c!='$'){

					uint64_t successor = OUT_[start_positions_out_[i]+off];

					IN_[start_positions_in_[successor]++] = i;


				}

			}

		}

		assert(start_positions_in_[nr_of_nodes-1] == IN_.size());

		//right-shift start_positions_in_

		for(uint64_t i=nr_of_nodes-1; i>0; --i){

			start_positions_in_[i] = start_positions_in_[i-1];

		}

		//the $ of the source starts at position 0 in IN
		start_positions_in_[0]=0;
		start_positions_in_[1]=1;

		MEAN_WEIGHT /= (weights_.size()-padded_kmers);

		if((not do_not_optimize) and (not XBWT))
			prune(verbose);

		//compute MST

		if(verbose)
			cout << "Computing MST ... " << endl;

		auto w = MST_(verbose);

		if(verbose)
			cout << " Weight of MST: " << w << " bits (" << double(w)/number_of_distinct_kmers() << " bits/kmer if stored with Elias Gamma)" << endl;

		if(XBWT)
			make_forest(verbose);

		if((not do_not_optimize) and XBWT)
			prune(verbose);


		if(D){

			//Count distinct abundances
			cout << "Counting distinct abundances ..." << endl;
			vector<uint32_t> W(weights_);
			sort( W.begin(), W.end() );
			W.erase( unique( W.begin(), W.end() ), W.end() );
			cout << "done. " << W.size() << " distinct abundances" << endl;


		}


		if(verbose)
			cout << "Computing deltas on the MST edges ... " << endl;

		build_deltas();

		if(verbose)
			cout << "Computing MST tree decomposition (samples) ... " << endl;

		macro_tree_decomposition(verbose);

		if(verbose)
			cout << "Compressing the data structures ... " << endl;

		compress_structures();

		if(verbose)
			cout << "Freeing memory ... " << endl;

		free_uncompressed_structures();

	}

	uint8_t get_order(){
		return k;
	}

	/*
	 * save the structure to the path specified.
	 */
	void save_to_file(string path){

		std::ofstream out(path);
		serialize(out);
		out.close();

	}

	uint64_t serialize(ostream& out){

		DBG_SIZE = 0;//size in bytes of the dBg
		DELTAS_SIZE = 0;//size in bytes of the deltas
		MST_SIZE = 0;//size in bytes of the MST + samples

		uint64_t written_bytes = 0;//total size in bytes

		out.write((char*)&k,sizeof(k));
		DBG_SIZE += sizeof(k);

		out.write((char*)&srate,sizeof(srate));
		DBG_SIZE += sizeof(srate);

		out.write((char*)&nr_of_nodes,sizeof(nr_of_nodes));
		DBG_SIZE += sizeof(nr_of_nodes);

		out.write((char*)&MAX_WEIGHT,sizeof(MAX_WEIGHT));
		DBG_SIZE += sizeof(MAX_WEIGHT);

		out.write((char*)&MEAN_WEIGHT,sizeof(MEAN_WEIGHT));
		DBG_SIZE += sizeof(MEAN_WEIGHT);

		out.write((char*)&padded_kmers,sizeof(padded_kmers));
		DBG_SIZE += sizeof(padded_kmers);

		char XBWT_c = XBWT;
		out.write(&XBWT_c,sizeof(XBWT_c));
		DBG_SIZE += sizeof(XBWT_c);

		DBG_SIZE += BWT.serialize(out);
		DBG_SIZE += IN.serialize(out);
		DBG_SIZE += IN_rank.serialize(out);
		DBG_SIZE += IN_sel.serialize(out);
		DBG_SIZE += OUT.serialize(out);
		DBG_SIZE += OUT_rank.serialize(out);
		DBG_SIZE += OUT_sel.serialize(out);
		out.write((char*)C.data(),5*sizeof(uint64_t));
		DBG_SIZE += 5*sizeof(uint64_t);

		written_bytes += DBG_SIZE;

		DELTAS_SIZE += deltas->serialize(out);

		written_bytes += DELTAS_SIZE;

		MST_SIZE += mst.serialize(out);
		MST_SIZE += sampled.serialize(out);
		MST_SIZE += sampled_rank.serialize(out);
		MST_SIZE += samples->serialize(out);

		written_bytes += MST_SIZE;

		return written_bytes;

	}


	/*
	 * load the structure from the path specified.
	 */
	void load_from_file(string path){

		std::ifstream in(path);
		load(in);
		in.close();

	}


	/* load the structure from the istream
	 */
	void load(std::istream& in) {

		in.read((char*)&k,sizeof(k));
		in.read((char*)&srate,sizeof(srate));
		in.read((char*)&nr_of_nodes,sizeof(nr_of_nodes));
		in.read((char*)&MAX_WEIGHT,sizeof(MAX_WEIGHT));
		in.read((char*)&MEAN_WEIGHT,sizeof(MEAN_WEIGHT));
		in.read((char*)&padded_kmers,sizeof(padded_kmers));

		char XBWT_c;
		in.read(&XBWT_c,sizeof(XBWT_c));
		XBWT = XBWT_c;

		BWT.load(in);

		IN.load(in);
		IN_rank.load(in);//TODO might have to rebuild all rank/select DS
		IN_rank = typename bitv_type::rank_1_type(&IN);
		IN_sel.load(in);
		IN_sel = typename bitv_type::select_1_type(&IN);

		OUT.load(in);
		OUT_rank.load(in);
		OUT_rank = typename bitv_type::rank_1_type(&OUT);
		OUT_sel.load(in);
		OUT_sel = typename bitv_type::select_1_type(&OUT);

		C = vector<uint64_t>(5,0);
		in.read((char*)C.data(),5*sizeof(uint64_t));

		deltas = new cint_vector();
		deltas->load(in);

		mst.load(in);

		sampled.load(in);
		sampled_rank.load(in);
		sampled_rank = rrr_vector<>::rank_1_type(&sampled);

		samples = new cint_vector();
		samples->load(in);

	}

	/*
	 * size of the de Bruijn graph
	 */
	uint64_t dbg_size_in_bits(){

		return	DBG_SIZE*8;

	}

	/*
	 * size (bits) of just the compressed deltas
	 */
	uint64_t deltas_size_in_bits(){

		return	DELTAS_SIZE*8;

	}

	/*
	 * size (bits) of the MST topology and samples
	 */
	uint64_t mst_size_in_bits(){

		return	MST_SIZE*8;

	}

	/*
	 * total size (Bytes) of the data structure
	 */
	uint64_t size_in_bits(){

		return dbg_size_in_bits() + deltas_size_in_bits() + mst_size_in_bits();

	}

	/*
	 * number of nodes (distinct kmers + padded kmers) in the de Bruijn graph. Note: also padded kmers are counted.
	 */
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

		return BWT.size();

	}

	uint64_t max_weight(){

		return MAX_WEIGHT;

	}

	double mean_weight(){

		return MEAN_WEIGHT;

	}

	/*
	 * first column of the BWT matrix
	 */
	char F(uint64_t i){

		return toCHAR(F_int(i));

	}

	/*
	 * returns node corresponding to input kmer (as string)
	 *
	 * if kmer does not exist, returns number_of_nodes()
	 *
	 */
	uint64_t find_kmer(string& kmer){

		assert(kmer.length()==k);

		auto range = full_range();

		for(auto c : kmer){

			//cout << "range: " << range.first << " " << range.second << endl;
			range = LF(range,c);

		}

		//cout << "range: " << range.first << " " << range.second << endl;

		if(range.second != range.first+1)
			return number_of_nodes();

		return range.first;

	}

	/*
	 * input: kmer
	 * output: abundance
	 */
	uint64_t operator[](string& kmer){

		return abundance(kmer);

	}

	/*
	 * input: kmer
	 * output: abundance
	 */
	uint64_t abundance(string& kmer){

		//auto idx = find_kmer(kmer);
		//return idx == number_of_nodes()?0:weights_[idx]; //TODO: use the compressed weights once 'weights' has been freed

		uint64_t n = find_kmer(kmer);

		if(n==number_of_nodes()) return 0;

		//cumulate deltas on the path from n to its nearest sampled ancestor
		int cumulated_deltas = 0;

		while(not sampled[n]){

			assert(n<deltas->size());

			cumulated_deltas -= positive_to_int(deltas->operator[](n));
			n = mst_parent(n);

		}
		assert(sampled[n]);

		auto sum = int64_t(samples->operator[](sampled_rank(n))) + cumulated_deltas;

		assert(sum>=0);

		return sum;

	}

	/*
	 * input: packed kmer
	 * output: abundance
	 */
	uint64_t abundance(__uint128_t& kmer, int k){

		string km = kmer_to_str_(kmer, k);

		return abundance(km);

	}

	/*
	 * edge is represented as pair <node, rank> where rank is the rank among outgoing edges of node
	 */
	int weight_of_edge_(edge_t e){

		assert(e.first<number_of_nodes());
		assert(e.second < out_degree_(e.first));

		assert(out_label_(e.first, e.second) != '$');

		auto target = OUT_[start_positions_out_[e.first] + e.second];

		assert(target<number_of_nodes());

		return abundance_(e.first) - abundance_(target);

	}

private:


	void free_uncompressed_structures(){

		//labels of outgoing edges. Corresponds to the BWT, except that there can be unnecessary $ here (will be removed in BWT)
		out_labels_.clear();
		out_labels_.shrink_to_fit();

		OUT_.clear(); //outgoing edges
		OUT_.shrink_to_fit();
		IN_.clear(); //incoming edges
		IN_.shrink_to_fit();

		start_positions_in_.clear(); //for each node, its starting position in IN_
		start_positions_in_.shrink_to_fit();
		start_positions_out_.clear(); //for each node, its starting position in OUT_ and out_labels_
		start_positions_out_.shrink_to_fit();

		weights_.clear(); //one weight per node
		weights_.shrink_to_fit();

		parent_in_mst_.clear(); //one per node. For the roots we write nr_of_nodes
		parent_in_mst_.shrink_to_fit();
		mst_out_edges_.clear(); //marks outgoing edges of the MST, in OUT_ order
		mst_out_edges_.shrink_to_fit();

		//deltas on MST input edges
		//for each node n, deltas_[n] is the weight of the edge leading to n, or 0 if n is a root.
		deltas_.clear();
		deltas_.shrink_to_fit();

		sampled_ = bit_vector(0); 	//marks sampled nodes on the dBg

		samples_.clear(); 	//weight samples
		samples_.shrink_to_fit();

	}


	/*
	 * returns parent of n in the MST
	 */
	uint64_t mst_parent(uint64_t n){

		auto in_deg = in_degree(n);

		assert(n>0);

		uint64_t first_pos = n==0?0:IN_sel__(n)+1;
		uint64_t last_pos = first_pos+in_deg-1;

		//there can be only 1 incoming edge in the MST

		while(first_pos <= last_pos && not mst__(first_pos)) first_pos++;

		assert(first_pos > last_pos || mst__(first_pos));

		return first_pos <= last_pos ? OUT_rank(FL(first_pos)) : nr_of_nodes;

	}

	/*
	 * returns the in-degree of the node
	 * since this is a dBg and there is only one source node ($$$),
	 * this number is always between 1 and 4
	 *
	 */
	uint8_t in_degree(uint64_t node){
		assert(node<number_of_nodes());
		return IN_sel__(node+1) - (node == 0 ? 0 : IN_sel__(node)+1) +1;
	}
	/*
	 * returns the out-degree of the node
	 * this number is always between 0 and 5 (count also $ if present)
	 *
	 */
	uint8_t out_degree(uint64_t node){
		assert(node<number_of_nodes());
		return OUT_sel(node+1) - (node == 0 ? 0 : OUT_sel(node)+1) +1;
	}

	/*
	 * input: node (represented as an integer) and character
	 * returns: node reached
	 *
	 * if node has no out-edge labeled c, the function returns number_of_nodes()
	 *
	 * c cannot be $
	 *
	 */
	uint64_t move_forward_by_char_(uint64_t node, char c){

		assert(c != '$');

		assert(node <= nr_of_nodes);

		if(node == nr_of_nodes) return node;

		uint64_t idx = start_positions_out_[node];
		uint8_t k = 0;

		while(out_labels_[idx+k] != c and k < out_degree_(node)) ++k;

		return k == out_degree_(node) ? nr_of_nodes : OUT_[idx+k];

	}

	/*
	 * input: node (represented as an integer) and rank between 0 and out_degree_(node)-1
	 * returns: node reached by following the k-th outgoing edge of the node
	 *
	 */
	uint64_t move_forward_by_rank_(uint64_t node, uint8_t k){

		if(node == nr_of_nodes) return node;

		assert(k<out_degree_(node));

		return OUT_[start_positions_out_[node]+k];

	}

	/*
	 * input: edge (represented as pair node, rank of outgoing edge)
	 * returns: position in OUT_ of the edge
	 *
	 */
	uint64_t edge_pos_in_OUT_(edge_t e){

		assert(e.second < out_degree_(e.first));
		return start_positions_out_[e.first] + e.second;

	}

	/*
	 * input: edge (represented as pair node, rank of outgoing edge)
	 * returns: position in array IN (equivalently, in F column) of the edge
	 *
	 */
	uint64_t edge_pos_in_IN_(edge_t e){

		assert(e.second < in_degree_(e.first));

		return start_positions_in_[e.first] + e.second;


	}

	/*
	 * input: node (represented as an integer) and rank between 0 and out_degree_(node)-1
	 * returns: label of k-th outgoing edge of the node
	 *
	 */
	char out_label_(uint64_t node, uint8_t k){

		assert(k<out_degree_(node));

		return out_labels_[start_positions_out_[node] + k];

	}

	/*
	 * input: node (represented as an integer) and rank between 0 and in_degree_(node)-1
	 * returns: node reached by following the k-th incoming edge of the node
	 *
	 * node cannot be the root 0 ($$$)
	 */
	uint64_t move_backward_(uint64_t node, uint8_t k){

		assert(k<in_degree_(node));
		assert(in_degree_(node)>0);

		assert(node > 0);
		assert(node < number_of_nodes());

		return IN_[start_positions_in_[node]+k];

	}

	/*
	 * returns the in-degree of the node
	 * since this is a dBg and there is only one source node ($$$),
	 * this number is always between 1 and 4
	 *
	 */
	uint8_t in_degree_(uint64_t node){

		assert(node<number_of_nodes());

		return node == nr_of_nodes-1 ? IN_.size() - start_positions_in_[node] : start_positions_in_[node+1]-start_positions_in_[node];

	}

	/*
	 * returns the out-degree of the node
	 * this number is always between 0 and 5 (count also $ if present)
	 *
	 */
	uint8_t out_degree_(uint64_t node){

		assert(node<number_of_nodes());

		return node == nr_of_nodes-1 ? OUT_.size() - start_positions_out_[node] : start_positions_out_[node+1]-start_positions_out_[node];

	}

	/*
	 * removes all the unnecessary padded nodes. The problem with dBgs over read sets is that we add k padded nodes
	 * for each read. Many of these however are not necessary if the dBg is well connected. In particular, only sources
	 * in the dBg (entry points in the connected components) need padded nodes, and those are few.
	 */
	void prune(bool verbose){


		if(verbose){
			cout << "Pruning ..." << endl;
		}

		uint64_t initial_nr_of_nodes = number_of_nodes();
		uint64_t initial_nr_of_edges = out_labels_.size();

		//mark all padded nodes
		vector<bool> padded(number_of_nodes(),false);

		//nodes that are necessary (not to be deleted)
		vector<bool> necessary_node(number_of_nodes(),false);

		//mark elements in OUT_ and IN_
		vector<bool> remove_out_edge(OUT_.size(),false);
		vector<bool> remove_in_edge(IN_.size(),false);

		{
			//find all nodes that are within distance k-1 from root (i.e. all padded nodes)

			//pairs <node, distance from root>
			stack<pair<uint64_t, uint8_t> > S;
			S.push({0,0}); //push the root (at distance 0)

			while(not S.empty()){

				auto p = S.top();
				S.pop();

				uint64_t n = p.first;
				uint8_t d = p.second;

				assert(d<k);

				padded[n] = true;

				for(uint8_t i = 0; i < out_degree_(n) && d<k-1 ;++i){

					if(out_label_(n,i) != '$'){

						S.push({move_forward_by_rank_(n,i),d+1});

						assert(in_degree_(move_forward_by_rank_(n,i))==1);
						assert(move_backward_(move_forward_by_rank_(n,i),0) == n);

					}

				}

			}

		}

		//now, for each non-padded kmer X that is preceded only by a padded kmer (i.e. a source of the dBg),
		//mark as necessary all padded kmers that lead from the root to X

		for(uint64_t n = 0; n<number_of_nodes();++n){

			if(not padded[n]){

				necessary_node[n] = true;

			}

			uint64_t node;
			if(not padded[n] && n>0 && in_degree_(n)==1 && padded[node = move_backward_(n,0)]){

				while(node != 0){

					//cout << "necessary padded node: " << node << endl;

					assert(in_degree_(node) == 1); //all padded nodes have in-degree 1
					assert(padded[node]);
					necessary_node[node] = true;
					node = move_backward_(node,0);

				}

			}

		}

		necessary_node[0] = true;

		//update padded_kmers
		padded_kmers = 0;
		for(uint64_t i=0;i<number_of_nodes();++i)
			padded_kmers += padded[i] and necessary_node[i];

		{

			string new_out_labels;
			vector<uint64_t> new_start_positions_in; //for each node, its starting position in IN_
			vector<uint64_t> new_start_positions_out; //for each node, its starting position in OUT_ and out_labels_

			uint64_t start_IN = 0;
			uint64_t start_OUT = 0;

			for(uint64_t n = 0;n<number_of_nodes();++n){

				if(necessary_node[n]){

					new_start_positions_in.push_back(start_IN);
					new_start_positions_out.push_back(start_OUT);

					//in- and out-degree of the node
					auto in_deg = in_degree_(n);
					auto out_deg = out_degree_(n);

					//compute new in-degree
					uint8_t new_in_deg = 0;
					for(uint8_t k = 0;k<in_deg;++k){

						if(n == 0 || necessary_node[move_backward_(n,k)]){

							new_in_deg++;

						}

					}

					start_IN += new_in_deg;

					assert(new_in_deg>0);

					//compute new out-degree
					uint8_t new_out_deg = 0;

					assert(out_deg>0);

					for(uint8_t k = 0;k<out_deg;++k){

						char c = out_label_(n,k);

						if(c == '$'){

							//insert $ in new_out_labels only if the node has out-degree = 1 (i.e. has only $ as outgoing edge)
							if(out_deg==1){
								new_out_deg++;
								new_out_labels.push_back(c);
							}

						}else{

							//cout << "edge " << char(c) << " node " << n << " out node: " << int(necessary_node[move_forward_by_rank_(n,k)]) << endl;

							if(necessary_node[move_forward_by_rank_(n,k)]){
								new_out_labels.push_back(c);
								new_out_deg++;
							}

						}

					}

					//the root is allowed to have out degree 1 because we mark it as necessary even if it might not be
					//we fix this in the next line
					assert(n == 0 or new_out_deg>0);

					//if root is not actually necessary, we add to it just 1 outgoing edge labeled $
					if(n==0 and new_out_deg == 0){
						new_out_labels.push_back('$');
						new_out_deg=1;
					}

					start_OUT += new_out_deg;

				}

			}

			out_labels_ = new_out_labels;
			start_positions_in_ = new_start_positions_in; //for each node, its starting position in IN_
			start_positions_out_ = new_start_positions_out; //for each node, its starting position in OUT_ and out_labels_

			assert(start_positions_in_.size() == start_positions_out_.size());

			nr_of_nodes = start_positions_in_.size();

		}

		//overwrite C

		C = vector<uint64_t>(5,0);

		uint64_t F_length = 1;

		for(auto c : out_labels_){

			C[toINT(c)]++;
			F_length += c!='$';

		}

		C[0] = 1;//there is one $ in the F column (it's the only incoming edge of the root kmer $$$)
		C[1] += C[0];
		C[2] += C[1];
		C[3] += C[2];
		C[4] += C[3];

		C[4] = C[3];
		C[3] = C[2];
		C[2] = C[1];
		C[1] = C[0];
		C[0] = 0; //$ starts at position 0 in the F column

		IN_ = vector<uint64_t>(F_length);
		IN_[0] = nr_of_nodes;
		OUT_ = vector<uint64_t>(out_labels_.size(),nr_of_nodes);

		//now fill IN_ and OUT_

		//current node in IN_, for each letter
		vector<uint64_t> curr_node_in(5);

		for(auto c : {'A','C','G','T'}){

			curr_node_in[toINT(c)] = distance(start_positions_in_.begin(),lower_bound(start_positions_in_.begin(),start_positions_in_.end(),C[toINT(c)]));
			assert(start_positions_in_[curr_node_in[toINT(c)]] == C[toINT(c)]);

		}

		uint64_t curr_node_out = 0;

		vector<uint64_t> curr_pos_in(5);

		for(auto c : {'A','C','G','T'}){

			curr_pos_in[toINT(c)] = C[toINT(c)];

		}

		for(uint64_t curr_pos_out=0;curr_pos_out<out_labels_.size();++curr_pos_out){

			//update curr_node_out if we reach the position of the next node
			if(curr_node_out < nr_of_nodes-1 && curr_pos_out == start_positions_out_[curr_node_out+1]) curr_node_out++;

			char c  = out_labels_[curr_pos_out];

			if(c != '$'){

				//update curr_node_in if we reach the position of the next node
				if(curr_node_in[toINT(c)] < nr_of_nodes-1 && curr_pos_in[toINT(c)] == start_positions_in_[curr_node_in[toINT(c)]+1]) curr_node_in[toINT(c)]++;

				OUT_[curr_pos_out] = curr_node_in[toINT(c)];

				assert(curr_pos_in[toINT(c)] != 0);
				IN_[curr_pos_in[toINT(c)]] = curr_node_out;

				curr_pos_in[toINT(c)]++;

			}

		}


		//prune weights

		uint64_t idx=0;

		for(uint64_t i=0;i<weights_.size();++i)
			if(necessary_node[i]) weights_[idx++] = weights_[i];

		assert(idx == number_of_nodes());

		weights_.resize(number_of_nodes());

		for(uint64_t n=0;n<nr_of_nodes;++n){

			assert(in_degree_(n)>0);
			assert(in_degree_(n)<=5);
			assert(out_degree_(n)<=5);
			assert(out_degree_(n)>0);

		}

		if(verbose){
			cout << " Removed " << initial_nr_of_nodes-number_of_nodes() << " nodes and " << initial_nr_of_edges-out_labels_.size() << " edges." << endl;
		}

	}



	/*
	 * turns the graph into a forest corresponding to the MST forest previously computed
	 */
	void make_forest(bool verbose){

		if(verbose){
			cout << "Turning graph into a forest ..." << endl;
		}

		uint64_t initial_nr_of_nodes = number_of_nodes();
		uint64_t initial_nr_of_edges = out_labels_.size();

		{

			string new_out_labels;
			vector<uint64_t> new_start_positions_in; //for each node, its starting position in IN_
			vector<uint64_t> new_start_positions_out; //for each node, its starting position in OUT_ and out_labels_
			vector<bool> new_mst_out_edges_;

			uint64_t start_IN = 0;
			uint64_t start_OUT = 0;

			for(uint64_t n = 0;n<number_of_nodes();++n){

				new_start_positions_in.push_back(start_IN);
				new_start_positions_out.push_back(start_OUT);

				//in- and out-degree of the node
				auto in_deg = in_degree_(n);
				auto out_deg = out_degree_(n);

				//it's a tree: in-degree always equal to 1 (for the roots we have an incoming edge labeled $)
				uint8_t new_in_deg = 1;

				start_IN += new_in_deg;

				//compute new out-degree
				uint8_t new_out_deg = 0;

				assert(out_deg>0);

				for(uint8_t k = 0;k<out_deg;++k){

					char c = out_label_(n,k);

					if(c == '$'){

						//previously we allowed $ only if it was the only outgoing label
						assert(out_deg==1);
						new_out_deg++;
						new_out_labels.push_back(c);
						new_mst_out_edges_.push_back(false);

					//keep outgoing edge only if it's part of the MST
					}else if(mst_out_edges_[start_positions_out_[n]+k]){

						new_out_labels.push_back(c);
						new_out_deg++;
						new_mst_out_edges_.push_back(true);

					}

				}

				//in this case it means that all children are already "captured" in the
				//MST by some other node. Add just one outgoing label $ to this node.
				if(new_out_deg==0){

					new_out_deg++;
					new_out_labels.push_back('$');
					new_mst_out_edges_.push_back(false);

				}

				start_OUT += new_out_deg;

			}

			out_labels_ = new_out_labels;
			start_positions_in_ = new_start_positions_in; //for each node, its starting position in IN_
			start_positions_out_ = new_start_positions_out; //for each node, its starting position in OUT_ and out_labels_
			mst_out_edges_ = new_mst_out_edges_;

			assert(start_positions_in_.size() == start_positions_out_.size());

			nr_of_nodes = start_positions_in_.size();

		}

		//overwrite C

		C = vector<uint64_t>(5,0);

		uint64_t F_length = 1;

		for(auto c : out_labels_){

			C[toINT(c)]++;
			F_length += c!='$';

		}

		C[0] = 1;//there is one $ in the F column (it's the only incoming edge of the root kmer $$$)
		C[1] += C[0];
		C[2] += C[1];
		C[3] += C[2];
		C[4] += C[3];

		C[4] = C[3];
		C[3] = C[2];
		C[2] = C[1];
		C[1] = C[0];
		C[0] = 0; //$ starts at position 0 in the F column

		IN_ = vector<uint64_t>(F_length);
		IN_[0] = nr_of_nodes;
		OUT_ = vector<uint64_t>(out_labels_.size(),nr_of_nodes);

		//now fill IN_ and OUT_

		//current node in IN_, for each letter
		vector<uint64_t> curr_node_in(5);

		for(auto c : {'A','C','G','T'}){

			curr_node_in[toINT(c)] = distance(start_positions_in_.begin(),lower_bound(start_positions_in_.begin(),start_positions_in_.end(),C[toINT(c)]));
			assert(start_positions_in_[curr_node_in[toINT(c)]] == C[toINT(c)]);

		}

		uint64_t curr_node_out = 0;

		vector<uint64_t> curr_pos_in(5);

		for(auto c : {'A','C','G','T'}){

			curr_pos_in[toINT(c)] = C[toINT(c)];

		}

		for(uint64_t curr_pos_out=0;curr_pos_out<out_labels_.size();++curr_pos_out){

			//update curr_node_out if we reach the position of the next node
			if(curr_node_out < nr_of_nodes-1 && curr_pos_out == start_positions_out_[curr_node_out+1]) curr_node_out++;

			char c  = out_labels_[curr_pos_out];

			if(c != '$'){

				//update curr_node_in if we reach the position of the next node
				if(curr_node_in[toINT(c)] < nr_of_nodes-1 && curr_pos_in[toINT(c)] == start_positions_in_[curr_node_in[toINT(c)]+1]) curr_node_in[toINT(c)]++;

				OUT_[curr_pos_out] = curr_node_in[toINT(c)];

				assert(curr_pos_in[toINT(c)] != 0);
				IN_[curr_pos_in[toINT(c)]] = curr_node_out;

				curr_pos_in[toINT(c)]++;

			}

		}

		//the number of nodes cannot change as this is a MST
		assert(initial_nr_of_nodes == number_of_nodes());

		for(uint64_t n=0;n<nr_of_nodes;++n){

			assert(in_degree_(n)>0);
			assert(in_degree_(n)<=5);
			assert(out_degree_(n)<=5);
			assert(out_degree_(n)>0);

		}

		assert(mst_out_edges_.size() == OUT_.size());

		if(verbose){
			cout << " Removed " << initial_nr_of_edges-out_labels_.size() << " edges." << endl;
		}

	}


	/*
	 * abundance of node n (represented as integer, i.e. its rank among all nodes)
	 */
	int abundance_(uint64_t n){

		assert(n<number_of_nodes());
		return int(weights_[n]);

	}

	/*
	 * returns character in position i of column F
	 */
	uint8_t F_int(uint64_t i){

		uint8_t  c = i>=C[1] and i<C[2] ? 1 :
				 i>=C[2] and i<C[3] ? 2 :
				 i>=C[3] and i<C[4] ? 3 :
				 i>=C[4]            ? 4 : 0;

		return c;

	}

	uint64_t LF(uint64_t i){

		char c = BWT[i];

		assert(c!='$');

		return C[toINT(c)] + BWT.rank(i,c);

	}

	pair<uint64_t, uint64_t> full_range(){
		return {0,number_of_nodes()};
	}

	/*
	 * LF function. Input: range OF NODES [begin, end) of nodes reached by path P, and char c
	 * output: nodes reached by path Pc.
	 *
	 * output is an empty range (end <= begin) if input is an empty range or if Pc does not occur
	 *
	 */
	pair<uint64_t, uint64_t> LF(pair<uint64_t, uint64_t> range, char c){

		assert(c!='$');

		uint64_t begin_on_BWT = range.first == 0 ? 0 : OUT_sel(range.first)+1; //inclusive
		uint64_t end_on_BWT = OUT_sel(range.second)+1; //exclusive

		//cout << "begin on BWT: " << begin_on_BWT << " " << end_on_BWT << endl;

		uint64_t c_before_begin = BWT.rank(begin_on_BWT,c);
		uint64_t c_before_end = BWT.rank(end_on_BWT,c);

		//cout << "c before: " << c_before_begin << " " << c_before_end << endl;

		uint64_t first_on_F = C[toINT(c)] + c_before_begin;
		uint64_t last_on_F = C[toINT(c)] + c_before_end;

		//cout << "on F " << first_on_F << " " << last_on_F << endl;

		return {IN_rank__(first_on_F), IN_rank__(last_on_F)};

	}

	uint64_t FL(uint64_t i){

		assert(i>0);
		assert(i<IN.size());

		uint8_t  c = F_int(i);

		assert(i>=C[c]);

		uint64_t sel = (i-C[c])+1;
		char ch = toCHAR(c);

		assert(sel<=BWT.rank(BWT.size(),ch));
		return BWT.select(sel,ch);

	}

	/*
	 * computes MST forest.
	 * returns weight of the MST
	 */
	uint64_t MST_(bool verbose){

		parent_in_mst_ = vector<uint64_t>(nr_of_nodes,nr_of_nodes);//roots i will have parent_in_mst_[i] = nr_of_nodes
		mst_out_edges_ = vector<bool>(OUT_.size(),false);

		uint64_t weight = 0; //weight of the MST

		set<uint64_t> not_in_mst;//nodes not yet in the MST

		int n_MST = 0; //number of trees in the MST forest

		//insert all nodes but the root in MST
		for(uint64_t u = 0; u<number_of_nodes(); ++u) not_in_mst.insert(u);

		while(not_in_mst.size()>0){

			n_MST++;

			uint64_t u = *not_in_mst.begin();

			not_in_mst.erase(not_in_mst.find(u));

			//pq contains edges (u,k). Let v = move_forward_by_rank_(u,k). Then edge (u,v) is not in the MST
			priority_queue<edge_t, vector<edge_t>, comp_edge<decltype(*this)> > pq(*this);

			for(uint8_t k = 0; k<out_degree_(u);++k)
				if(out_label_(u,k)!='$')
					pq.push({u,k});

			while(not pq.empty()){

				edge_t e;
				uint64_t v = 0;

				//find a minimum weight edge on the frontier
				do{

					e = pq.top();
					pq.pop();

					v = move_forward_by_rank_(e.first,e.second);

				}while((not pq.empty()) && not_in_mst.find(v) == not_in_mst.end());

				//frontier edge found
				if(not_in_mst.find(v) != not_in_mst.end()){

					weight += cost_of_weight(weight_of_edge_(e));//weight of MST

					not_in_mst.erase(not_in_mst.find(v));
					for(uint8_t k = 0; k<out_degree_(v);++k)
						if(out_label_(v,k)!='$')
							pq.push({v,k});

					parent_in_mst_[v] = e.first;
					mst_out_edges_[start_positions_out_[e.first] + e.second] = true;

				}

			}

		}

		if(verbose)
			cout << " Number of connected components : " << n_MST << endl;

		assert(not_in_mst.size()==0);

		return weight;

	}

	void build_deltas(){

		//for each node
		for(uint64_t i = 0; i<nr_of_nodes; ++i){

			//if the node is not the root of a tree in the MST forest
			if(parent_in_mst_[i] != nr_of_nodes){

				int w1 = abundance_(parent_in_mst_[i]);
				int w2 = abundance_(i);

				auto encoded_diff = int_to_positive(w1-w2);

				deltas_.push_back(encoded_diff);

			}else{//push 0

				deltas_.push_back(0);

			}

		}

		assert(deltas_.size() == number_of_nodes());

	}

	/*
	 * is n a leaf of the MST?
	 */
	bool is_leaf_(uint64_t n){

		uint8_t out_deg = out_degree_(n);

		assert(out_deg>0);

		bool out_edge_found = false;

		uint64_t start = start_positions_out_[n];

		for(uint8_t k=0; k < out_deg;++k){

			out_edge_found = out_edge_found or mst_out_edges_[start+k];

		}

		return not out_edge_found;

	}

	/*
	 * true iff n is a root in the mst forest
	 */
	bool is_root_in_mst_(uint64_t n){

		assert(n < nr_of_nodes);

		return parent_in_mst_[n] == nr_of_nodes;

	}

	/*
	 * marks sampled nodes on the dBg using a tree decomposition that
	 * decomposes the tree in subtrees of size Theta(srate). In the end, marks
	 * the roots of the subtrees using bitvector 'sampled'
	 */
	void macro_tree_decomposition(bool verbose){

		//rrr_vector<> sampled;
		//rrr_vector<>::rank_1_type sampled_rank;

		//nodes in post-order
		vector<uint64_t> post_order(nr_of_nodes);
		int64_t idx = nr_of_nodes-1;

		{

			stack<uint64_t> nodes;

			//push all MST roots in the stack
			for(uint64_t i=0;i<nr_of_nodes;++i)
				if(is_root_in_mst_(i))
					nodes.push(i);

			while(not nodes.empty()){

				auto n = nodes.top();
				nodes.pop();

				assert(idx < nr_of_nodes);
				assert(idx >= 0);
				post_order[idx--] = n;

				for(uint8_t k = 0; k < out_degree_(n);++k){

					uint64_t child_pos = start_positions_out_[n]+k;

					if(mst_out_edges_[child_pos])
						nodes.push(OUT_[child_pos]);

				}

			}

			assert(idx==-1);

		}

		sampled_ = bit_vector(nr_of_nodes,false);

		//size of each component of the tree decomposition
		//at the beginning, each component contains just the node itself.
		vector<uint64_t> component_size(nr_of_nodes,1);

		for(uint64_t i=0;i<nr_of_nodes;++i){

			auto n = post_order[i];

			if(is_root_in_mst_(n)){

				sampled_[n] = 1;//roots are always sampled

			}else{

				if(component_size[n] < srate){

					component_size[parent_in_mst_[n]] += component_size[n];

				}else{

					sampled_[n] = 1;

				}

			}

		}

		assert(sampled_[0] == 1);

		for(uint64_t i=0;i<number_of_nodes();++i){

			if(sampled_[i])
				samples_.push_back(weights_[i]);

		}

		if(verbose)
			cout << " " << samples_.size() << " sampled weights" << endl;

	}

	void compress_structures(){

		//str_type BWT;

		construct_im(BWT,out_labels_,1);

		//vector<uint64_t> C; //C array (encoding F column of BWT matrix)

		//bitv_type IN;
		//typename bitv_type::rank_1_type IN_rank;
		//typename bitv_type::select_1_type IN_sel;

		if(not XBWT){

			bit_vector IN_bv(IN_.size());

			uint64_t idx=0;

			for(uint64_t n=0;n<number_of_nodes();++n){

				auto in_deg = in_degree_(n);

				assert(in_deg>0);

				for(uint8_t off = 0; off < in_deg-1; ++off){

					IN_bv[idx++] = 0;

				}

				IN_bv[idx++] = 1;

			}

			assert(idx == IN_.size());

			IN = bitv_type(IN_bv);
			IN_rank = typename bitv_type::rank_1_type(&IN);
			IN_sel = typename bitv_type::select_1_type(&IN);

		}

		//bitv_type OUT;
		//typename bitv_type::rank_1_type OUT_rank;
		//typename bitv_type::select_1_type OUT_sel;

		{

			bit_vector OUT_bv(OUT_.size());

			uint64_t idx=0;

			for(uint64_t n=0;n<number_of_nodes();++n){

				auto out_deg = out_degree_(n);

				assert(out_deg>0);

				for(uint8_t off = 0; off < out_deg-1; ++off){

					OUT_bv[idx++] = 0;

				}

				OUT_bv[idx++] = 1;

			}

			assert(idx == OUT_.size());

			OUT = bitv_type(OUT_bv);
			OUT_rank = typename bitv_type::rank_1_type(&OUT);
			OUT_sel = typename bitv_type::select_1_type(&OUT);

		}

		//cint_vector deltas; //compressed deltas on the edges of the MST, in IN order (i.e. on the F column)

		deltas = new cint_vector(deltas_);

		//cout << "******* cost gamma = " << cost_gamma(deltas_) << endl;

		if(not XBWT){

			//marks edges (in IN order) that belong to the MST
			//rrr_vector<> mst;

			bit_vector mst_bv(IN_.size(),0);

			uint64_t idx=0;

			for(uint64_t n=0;n<number_of_nodes();++n){

				auto in_deg = in_degree_(n);

				assert(in_deg>0);

				for(uint8_t off = 0; off < in_deg; ++off){

					if(parent_in_mst_[n] != nr_of_nodes and parent_in_mst_[n] == IN_[start_positions_in_[n]+off]){

						mst_bv[idx++] = 1;

					}else{

						mst_bv[idx++] = 0;

					}

				}

			}

			assert(idx == IN_.size());

			mst = rrr_vector<>(mst_bv);

		}

		//marks sampled nodes on the dBg
		//rrr_vector<> sampled;
		//rrr_vector<>::rank_1_type sampled_rank;

		sampled = rrr_vector<>(sampled_);
		sampled_rank = rrr_vector<>::rank_1_type(&sampled);

		//sampled weights
		//cint_vector samples;

		samples = new cint_vector(samples_);

	}

	/*
	 * these functions are an interface over IN, IN_rank, IN_sel, mst: if XBWT is true, these structures are
	 * not actually necessary
	 */

	bool IN__(uint64_t i){

		return XBWT ? true : IN[i];

	}

	uint64_t IN_rank__(uint64_t i){

		return XBWT ? i : IN_rank(i);

	}

	uint64_t IN_sel__(uint64_t i){

		assert(i>0);

		return XBWT ? i-1 : IN_sel(i);

	}

	bool mst__(uint64_t i){

		return XBWT ? true : mst[i];

	}


	//parameters:

	uint8_t k; //order
	uint64_t nr_of_nodes; //number of nodes
	uint64_t padded_kmers = 0;//number of nodes corresponding to a k-mers padded with $

	/*
	 * temporary fast data structures used during construction
	 * convention: we end names by underscore (as well as functions that use them)
	 */

	//labels of outgoing edges. Corresponds to the BWT, except that there can be unnecessary $ here (will be removed in BWT)
	string out_labels_;
	vector<uint64_t> OUT_; //outgoing edges
	vector<uint64_t> IN_; //incoming edges

	vector<uint64_t> start_positions_in_; //for each node, its starting position in IN_
	vector<uint64_t> start_positions_out_; //for each node, its starting position in OUT_ and out_labels_

};

}

#endif
