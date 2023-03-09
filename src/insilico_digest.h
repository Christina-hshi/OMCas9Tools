/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * ============================================================================
 */

#pragma once
#ifndef INSILICO_DIGEST_H
#define INSILICO_DIGEST_H

#include "base.h"
#include "SeqIO.h"
#include "Params.h"

// We don't accept non-zero mismatch allowed for restriction enzyme at this stage
void insilico_digest_restriction(DNASeq& dnaseq, Enzyme& enzyme, vector<OMLabel>& labels){
  if(enzyme.mismatch_max > 0){
    MSG::error("We don't accept non-zero mismatch allowed for restriction enzyme at this stage.");
  }
  const int elen = enzyme.length();
  const int dlen = dnaseq.length();
  const char* enz_ptr = enzyme.seq.c_str();
  const char* seq_ptr = dnaseq.seq.c_str();
  int seqidx = 0;
  while(seqidx <= dlen - elen){
    if (memcmp(seq_ptr, enz_ptr, elen) == 0){
      labels.push_back(OMLabel(dnaseq.name, seqidx, DNA::STRAND::PLUS, enzyme.channel));
    }
    seqidx++;
    seq_ptr++;
  }
  //find target in RC seq
  string enzymeRC = DNA::RC(enzyme.seq);
  if (enzyme.seq == enzymeRC){ return; } //skip RC search if the recognition sequence is panlindrome
  const char* enzRC_ptr = enzymeRC.c_str();
  seq_ptr = dnaseq.seq.c_str();
  seqidx = 0;
  while(seqidx <= dlen - elen){
    if (memcmp(seq_ptr, enzRC_ptr, elen) == 0){
      labels.push_back(OMLabel(dnaseq.name, seqidx+elen-1, DNA::STRAND::MINUS, enzyme.channel));
    }
    seqidx++;
    seq_ptr++;
  }
}

// PAM NGG is a must for CRISPR Cas9 digestion
void insilico_digest_Cas9(DNASeq& dnaseq, Enzyme& enzyme, vector<OMLabel>& labels){
  //find NGG and then compare with the enzyme
  int plen=enzyme.length();
  int seqidx = plen+2;
  while(seqidx < dnaseq.length()){
    if(dnaseq.seq[seqidx-1] != 'G' or dnaseq.seq[seqidx] != 'G'){
      seqidx++;
      continue;
    }
    //compare with the enzyme
    int mismatch=0;
    for(int x = seqidx-plen-2, y=0; y < plen; x++, y++){
      if(dnaseq.seq[x] != enzyme.seq[y]){
        mismatch++;
        if(mismatch > enzyme.mismatch_max) break;
      }
    }
    if(mismatch <= enzyme.mismatch_max){
      labels.push_back(OMLabel(dnaseq.name, seqidx-plen-2+enzyme.digest_loc, DNA::STRAND::PLUS, enzyme.channel));
    }
    seqidx++;
  }
  //find target in RC seq
  string enzymeRC = DNA::RC(enzyme.seq);
  if (enzyme.seq == enzymeRC){ return; } //skip RC search if the recognition sequence is panlindrome
  seqidx = 0;
  while(seqidx < dnaseq.length()-plen-2){
    if(dnaseq.seq[seqidx] != 'C' or dnaseq.seq[seqidx+1] != 'C'){
      seqidx++;
      continue;
    }
    //compare enzyme
    int mismatch=0;
    for(int x = seqidx+3, y=0; y < plen; x++, y++){
      if(dnaseq.seq[x] != enzymeRC[y]){
        mismatch++;
        if(mismatch > enzyme.mismatch_max) break;
      }
    }
    if(mismatch <= enzyme.mismatch_max){
      labels.push_back(OMLabel(dnaseq.name, seqidx+plen+2-enzyme.digest_loc, DNA::STRAND::MINUS, enzyme.channel));
    }
    seqidx++;
  }
  return ;
}

void insilico_digest(DNASeq& dnaseq, Enzyme& enzyme, vector<OMLabel>& labels){
  switch (enzyme.type) {
    case ENZYME_TYPE::RESTRICTION:
      insilico_digest_restriction(dnaseq, enzyme, labels);
      break;
    case ENZYME_TYPE::CRISPR_CAS9:
      insilico_digest_Cas9(dnaseq, enzyme, labels);
      break;
    default:
      MSG::error("Unrecognized enzyme type: " + enzyme.type);
      break;
  }
}

void insilico_digest(Params_insilicoDigest& params){
  uint64_t labelNumTotal = 0;
  //Initialize the output cmap file
  ofstream fout;
  fout.open(params.mapfile, std::ios::out);
  if(not fout.is_open()){
    MSG::error("Failed to open file "+params.mapfile);
  }
  //TODO: aggregate cross channels
  fout<<"# CMAP File Version:\t0.1\n";
    // <<"# Label Channels:\t1\n";
    //<<"# Nickase Recognition Site 1:\t"<<endl;
  int channel_count = 1;
  for(auto enzyme : params.enzymes){
    fout<<"# Label channel "<<channel_count<<" : "<<enzyme.seq<<endl;
    channel_count++;
  }
  fout<<"# Number of Consensus Nanomaps:\tUNKNOWN"<<endl
    <<boost::algorithm::join(vector<string>{"#h CMapId", "ContigLength", "NumSites",
                  "SiteID", "LabelChannel", "Position",
                  "StdDev", "Coverage", "Occurrence"}, "\t")<<endl
    <<boost::algorithm::join(vector<string>{"#f int", "float", "int", "int", "int", "float", "float", "int", "int"}, "\t")<<endl;

  DNASeqIO seqIn(params.seqfiles);
  seqIn.start();
  DNASeq dnaseq = seqIn.read_next();
  int cMapId = 1; //only count the ones with labels.
  while(not dnaseq.is_empty()){
    vector<OMLabel> labels;
    //digest by all gRNAs
    for(auto enzyme : params.enzymes){
      insilico_digest(dnaseq, enzyme, labels);
    }
    //sort labels
    std::sort(labels.begin(), labels.end(), [](const OMLabel& x1, const OMLabel& x2){
      return x1.loc < x2.loc;
    });
    //merge labels being too close to each other
    vector<OMLabel> labels_merged;
    vector<int> cluster_sizes; // only count in clusters with at least two labels

    int cluster_start=0, cluster_end=0;
    while(cluster_end < labels.size()){
      cluster_end++;
      if ((cluster_end == labels.size()) or (labels[cluster_end].loc - labels[cluster_end-1].loc > params.merge_dis)){
        int cluster_center=0;
        if(cluster_end - cluster_start > 1){
          for(int x = cluster_start; x < cluster_end; x++){
            cluster_center += labels[x].loc;
          }
          cluster_center /= (cluster_end - cluster_start);
          cluster_sizes.push_back(labels[cluster_end-1].loc-labels[cluster_start].loc);
        }else{
          cluster_center = labels[cluster_start].loc;
        }
        labels_merged.push_back(OMLabel(labels[cluster_start].ref, cluster_center,
          labels[cluster_start].strand, labels[cluster_start].channel));
        cluster_start = cluster_end;
      }
    }

    cout<<  "[Info] On "<<dnaseq.name<<", merged "<<cluster_sizes.size()<<" clusters of labels <= "<<params.merge_dis<<" bp. "<<endl;
    if(cluster_sizes.size() > 0){
      cout<<"       Span of clusters: min. "<<*min_element(cluster_sizes.begin(), cluster_sizes.end())
          <<"       max. "<<*max_element(cluster_sizes.begin(), cluster_sizes.end())
          <<"       ave. "<<float(accumulate(cluster_sizes.begin(), cluster_sizes.end(), 0))/cluster_sizes.size()<<endl;
    }

    int numSites = labels_merged.size();
    if(numSites == 0){
      dnaseq = seqIn.read_next();
      continue;
    }
    labelNumTotal += numSites;

    int siteID, stdDev, coverage, occurrence;
    siteID = stdDev = coverage = occurrence = 1;
    //label.loc converts to 1-based coordinate.
    for(OMLabel& label : labels_merged){
      fout<<cMapId<<"\t"<<dnaseq.length()<<"\t"<<numSites<<"\t"
        <<siteID<<"\t"<<label.channel<<"\t"<<label.loc+1<<"\t"
        <<stdDev<<"\t"<<coverage<<"\t"<<occurrence<<endl;
      siteID++;
    }
    fout<<cMapId<<"\t"<<dnaseq.length()<<"\t"<<numSites<<"\t"
      <<siteID<<"\t"<<0<<"\t"<<dnaseq.length()<<"\t"
      <<0<<"\t"<<1<<"\t"<<0<<endl;
    cMapId++;
    dnaseq = seqIn.read_next();
  }
  seqIn.close();
  fout.close();
  cout<<labelNumTotal<<" labels in total."<<endl;
}

#endif
