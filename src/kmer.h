/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * Description:
 *   A 64-bit implementation can handle k-mer with size up tp 32 bases.
 * ============================================================================
 */

#pragma once
#ifndef KMER_H
#define KMER_H

#include "base.h"

namespace DNA{
  // const char bases[4]={'A', 'C', 'G', 'T'};
  //Map bases to 2 bit representation: A|a:00; C|c:01; G|g:10; T|t:11; N|n:arbitrary.
  const uint8_t bits2RCbits[4] = {0b11, 0b10, 0b01, 0b00};
  const uint8_t base2bits[]={0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x0F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x1F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x2F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x3F
													0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x4F
													0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x5F
													0x00, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x6F
													0x00, 0x00, 0x00, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x7F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x8F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0x9F
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xAF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xBF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xCF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xDF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,//0xEF
													0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};//0xFF
  const char bits2base[]={'A', 'C', 'G', 'T', 'C', ' ', ' ', ' ', 'G', ' ', ' ', ' ', 'T', ' ', ' ', ' ',//0x0F
  												'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x1F
  												'G', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x2F
  												'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x3F
  												'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x4F
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x5F
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x6F
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x7F
  												'G', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x8F
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x9F
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xAF
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xBF
  												'T', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xCF
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xDF
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xEF
  												' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};//0xFF
  const char bits2RCbase[]={'T', 'G', 'C', 'A', 'G', ' ', ' ', ' ', 'C', ' ', ' ', ' ', 'A', ' ', ' ', ' ',//0x0F
                					'G', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x1F
                					'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x2F
                					'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x3F
                					'G', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x4F
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x5F
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x6F
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x7F
                					'C', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x8F
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x9F
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xAF
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xBF
                					'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xCF
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xDF
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xEF
                					' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};//0xFF
}

typedef uint64_t kmer_t;

kmer_t encode(const char* seq, int k=20){
	kmer_t res=0;
	for(int x =0; x < k; x++, seq++){
		res <<= 2;
		switch(*seq){
			case 'A':
			case 'a':
				break;
			case 'T':
			case 't':
				res |= 0b11;
				break;
			case 'C':
			case 'c':
				res |= 0b01;
				break;
			case 'G':
			case 'g':
				res |= 0b10;
				break;
			default:
				std::cerr<<"[error] unsupported base: "<<*seq<<" in "<<x<<" base of the sequence "<<(*seq)<<std::endl;
				exit(0);
		}
	}
	return res;
}

string decode(kmer_t kmer, int k=20){
	const char idx2base[4]= {'A', 'C', 'G', 'T'};
	string res(k, ' ');
	for(int x = k-1; x >= 0; x--){
		int idx = (kmer & 0b11);
		res[x] = idx2base[idx];
		kmer >>= 2;
	}
	return res;
}

kmer_t RC(kmer_t kmer, int k){
	kmer_t rc=0;
	for(int x = 0; x < k; x++){
		rc <<= 2;
		switch(kmer&0b11){
			case 0:
				rc |= 0b11;
				break;
			case 1:
				rc |= 0b10;
				break;
			case 2:
				rc |= 0b01;
				break;
			case 3:
				break;
			default:
				std::cerr<<"[error] unexpected in RC(kmer_t)"<<endl;
				break;
		}
		kmer >>= 2;
	}
	return rc;
}

//user can choose to save gRNA in plain ATGC bases, or encoded forms.
void save(const unordered_map<kmer_t, vector<GLoc>>& kmer2loci, const string out_file, const int k, const bool plain = false, const uint32_t freq_min = 0, const uint32_t freq_max = 0xFFFFFFFF){
	ofstream fout;
	fout.open(out_file, ios::out);
	if(!fout.is_open()){
		std::cerr<<"[Error] failed to open "<<out_file<<endl;
	}
	fout<<"gRNA\tfreq.\tloci\n";
	for(auto record : kmer2loci){
    auto loci_num = record.second.size();
    if(loci_num < freq_min or loci_num > freq_max){
      continue;
    }
		fout<<(plain?decode(record.first, k):std::to_string(record.first))<<"\t"<<loci_num<<"\t";
		fout<<record.second[0].ref<<(record.second[0].strand==DNA::STRAND::PLUS?"+":"-")<<record.second[0].loc;
		for(int idx = 1; idx < record.second.size(); idx++){
			fout<<"|"<<record.second[idx].ref<<(record.second[idx].strand==DNA::STRAND::PLUS?"+":"-")<<record.second[idx].loc;
		}
		fout<<endl;
	}
	fout.close();
	return ;
}

void load(unordered_map<kmer_t, vector<GLoc>>& kmer2loci, const string kmer_file, const int k, const bool plain = false, const uint32_t freq_min = 0, const uint32_t freq_max = 0xFFFFFFFF){
	ifstream fin;
	fin.open(kmer_file, ios::in);
	if(!fin.is_open()){
		std::cerr<<"[Error] failed to open "<<kmer_file<<endl;
	}
	string gRNA_seq, loci_info, line;
	uint32_t freq;
	getline(fin, line);
	while(getline(fin, line)){
		stringstream ss(line);
		ss>>gRNA_seq>>freq>>loci_info;
    if(freq < freq_min or freq > freq_max){
      continue;
    }
		vector<GLoc> loci;
		int pos;
		string chr_name;
		DNA::STRAND sense;
		int tmp1, tmp2;
		tmp1 = 0;
		tmp2 = loci_info.find('|', tmp1);
		while(tmp2 != string::npos){
			int tmp = tmp1;
			while(loci_info[tmp] != '+' && loci_info[tmp] != '-'){
				tmp++;
			}
			if(loci_info[tmp] == '+'){
				sense = DNA::STRAND::PLUS;
			}else{
				sense = DNA::STRAND::MINUS;
			}
			chr_name = loci_info.substr(tmp1, tmp-tmp1);
			pos = std::stoi(loci_info.substr(tmp+1, tmp2 - tmp - 1));
			loci.push_back(GLoc(chr_name, pos, sense));

			tmp1 = tmp2 + 1;
			tmp2 = loci_info.find('|', tmp1);
		}
		//handle the last loc
		tmp2 = loci_info.length();
		int tmp = tmp1;
		while(loci_info[tmp] != '+' && loci_info[tmp] != '-'){
			tmp++;
		}
		if(loci_info[tmp] == '+'){
			sense = DNA::STRAND::PLUS;
		}else{
			sense = DNA::STRAND::MINUS;
		}
		chr_name = loci_info.substr(tmp1, tmp-tmp1);
		pos = std::stoi(loci_info.substr(tmp+1, tmp2 - tmp - 1));
		loci.push_back(GLoc(chr_name, pos, sense));

		kmer_t kmer;
		if(plain){
			kmer = encode(gRNA_seq.c_str(), k);
		}else{
			kmer = std::stoull(gRNA_seq);
		}
		kmer2loci[kmer] = loci;
	}
	fin.close();
	return ;
}

#endif
