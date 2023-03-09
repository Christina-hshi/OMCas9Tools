/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * ============================================================================
 */

#pragma once
#ifndef FIND_20MERNGG_H
#define FIND_20MERNGG_H

#include "SeqIO.h"
#include "Params.h"
#include "kmer.h"

void find_20merNGG(Params_find20merNGG& params){
  const int k=20;
  unordered_map<kmer_t, vector<GLoc>> gRNA2loci;

  DNASeqIO seqIn(params.seqfiles);
  seqIn.start();
  DNASeq dnaseq = seqIn.read_next();
  while(not dnaseq.is_empty()){
    int seq_len = dnaseq.length();
    string& seq = dnaseq.seq;
    // check forward strand
		for(int x = k+1, y = 0; x < seq_len-1; x++, y++){
			if(seq[x] == 'G' && seq[x+1] == 'G'){
				//check if any 'N'
				bool noN = true;
				for(int pos = x-1; pos >= y; pos--){
					if(seq[pos] == 'N'){
						y = pos;
						x = pos + k + 1;
						noN = false;
						break;
					}
				}
				if(noN){
					kmer_t gRNA = encode(seq.c_str()+y, k);
					gRNA2loci[gRNA].push_back(GLoc(dnaseq.name, y, DNA::STRAND::PLUS));
				}
			}
		}
    // check reverse strand
		for(int x = 0; x <= seq_len - k - 3; x++){
			if(seq[x] == 'C' && seq[x+1] == 'C'){
				//check if any 'N'
				bool noN = true;
				for(int pos = x+k+2; pos > x+1; pos--){
					if(seq[pos] == 'N'){
						x = pos;
						noN = false;
						break;
					}
				}
				if(noN){
					string gRNA_seq = seq.substr(x+3, k);
					gRNA_seq = DNA::RC(gRNA_seq);
					kmer_t gRNA = encode(gRNA_seq.c_str(), k);
					gRNA2loci[gRNA].push_back(GLoc(dnaseq.name, x+k+2, DNA::STRAND::MINUS));
				}
			}
		}
    // next read
    dnaseq = seqIn.read_next();
  }
  // save kmerNGG to file
  save(gRNA2loci, params.output, k);
}

//find all
void find_20merNGG(GRange& region, unordered_map<kmer_t, vector<GLoc>>& gRNA2loci){
  const int k=20;
  if(region.seq == ""){
    MSG::error("Seq was not set for region to find 20merNGG.");
  }

  int seq_len = region.seq.length();
  string& seq = region.seq;
  // check forward strand
  for(int x = k+1, y = 0; x < seq_len-1; x++, y++){
    if(seq[x] == 'G' && seq[x+1] == 'G'){
      //check if any 'N'
      bool noN = true;
      for(int pos = x-1; pos >= y; pos--){
        if(seq[pos] == 'N'){
          y = pos;
          x = pos + k + 1;
          noN = false;
          break;
        }
      }
      if(noN){
        kmer_t gRNA = encode(seq.c_str()+y, k);
        gRNA2loci[gRNA].push_back(GLoc(region.ref, region.start+y, DNA::STRAND::PLUS));
      }
    }
  }
  // check reverse strand
  for(int x = 0; x <= seq_len - k - 3; x++){
    if(seq[x] == 'C' && seq[x+1] == 'C'){
      //check if any 'N'
      bool noN = true;
      for(int pos = x+k+2; pos > x+1; pos--){
        if(seq[pos] == 'N'){
          x = pos;
          noN = false;
          break;
        }
      }
      if(noN){
        string gRNA_seq = seq.substr(x+3, k);
        gRNA_seq = DNA::RC(gRNA_seq);
        kmer_t gRNA = encode(gRNA_seq.c_str(), k);
        gRNA2loci[gRNA].push_back(GLoc(region.ref, region.start+x+k+2, DNA::STRAND::MINUS));
      }
    }
  }
}

#endif
