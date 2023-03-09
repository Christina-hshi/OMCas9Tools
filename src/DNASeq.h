/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * ============================================================================
 */

#pragma once

#ifndef DNASEQ_H
#define DNASEQ_H

#include "base.h"

namespace DNA{
  const char bases[5]={'A', 'C', 'G', 'T', 'N'}; //allow N
  enum STRAND{PLUS, MINUS};
  const char base2RCbase[]={' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x0F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x1F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x2F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x3F
													' ', 'T', ' ', 'G', ' ', ' ', ' ', 'C', ' ', ' ', ' ', ' ', ' ', ' ', 'N', ' ',//0x4F
													' ', ' ', ' ', ' ', 'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x5F
													' ', 'T', ' ', 'G', ' ', ' ', ' ', 'C', ' ', ' ', ' ', ' ', ' ', ' ', 'N', ' ',//0x6F
													' ', ' ', ' ', ' ', 'A', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x7F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x8F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0x9F
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xAF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xBF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xCF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xDF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',//0xEF
													' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};//0xFF
  inline string RC(const string& seq){
    size_t seqlen = seq.length();
    string seqRC(seqlen, ' ');
    for(int x = 0, y = seqlen-1; x < seqlen; x++, y--){
      seqRC[y] = DNA::base2RCbase[seq[x]];
    }
    return seqRC;
  }
}

class DNASeq{
public:
  string name;
  string seq;
  string qual;//quality string, maybe empty for fasta file
  DNASeq(string n="", string s="", string q=""):name(n), seq(s), qual(q){}
  inline string RC_seq(){ return DNA::RC(this->seq); }
  inline bool is_empty(){
    return seq==""?true:false;
  }
  inline size_t length(){ return seq.length(); }
  void toBinary();
};

enum ENZYME_TYPE{RESTRICTION=1, CRISPR_CAS9=2};

class Enzyme{
public:
  string name;
  string seq;
  uint32_t digest_loc;
  uint32_t mismatch_max;
  ENZYME_TYPE type;
  int channel;
  Enzyme(string n, string s, uint32_t loc=0, uint32_t m=0, ENZYME_TYPE t=ENZYME_TYPE::RESTRICTION, int c=1):
    name(n), seq(s), digest_loc(loc), mismatch_max(m), type(t), channel(c) {}
  inline size_t length() { return seq.length(); }
  static vector<Enzyme> load_from_tsv(string filename){
    vector<Enzyme> enzymes;
    fstream fin;
    fin.open(filename, ios::in);
    if(not fin.is_open()){ MSG::error("Failed to open file "+filename); }
    string enzyme_name, enzyme_seq;
    int digest_site, mis_max, enzyme_type, chan;
    string line;
    getline(fin, line);//escape the header line
    while(fin>>enzyme_name>>enzyme_seq>>digest_site>>mis_max>>enzyme_type>>chan){
      enzymes.push_back(Enzyme(enzyme_name, enzyme_seq, digest_site, mis_max, static_cast<ENZYME_TYPE>(enzyme_type), chan));
    }
    fin.close();
    return enzymes;
  }
  string to_string(){
    return "Enzyme[]";
  }
};

// class Cas9Probe{
// public:
//   string name;
//   string seq;//not include PAM
//   uint32_t digest_loc;//nicking location relative to the start of the seq
//   int channel;
//   Cas9Probe(string n, string s, uint32_t loc=0, int c=1):name(n), seq(s), digest_loc(loc), channel(c) {}
//   inline size_t length() { return seq.length(); }
// };

// genomic location
class GLoc{
public:
  string ref;
  uint32_t loc;
  DNA::STRAND strand;
  GLoc(string r, uint32_t l, DNA::STRAND s): ref(r), loc(l), strand(s){}
};

class OMLabel: public GLoc{
public:
  int channel;
  OMLabel(string r, uint32_t l, DNA::STRAND s, int c):GLoc(r, l, s), channel(c){}
};

// genomic range
class GRange{
public:
  string ref;
  uint32_t start;
  uint32_t end;
  DNA::STRAND strand;
  string seq; //optional for space efficiency
  GRange(string r, uint32_t st, uint32_t en, DNA::STRAND s=DNA::STRAND::PLUS): ref(r), start(st), end(en), strand(s), seq(""){}
};

#endif
