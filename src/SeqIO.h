/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * ============================================================================
 */

#pragma once

#ifndef SEQIO_H
#define SEQIO_H

#include "base.h"
#include "DNASeq.h"

enum SEQ_FILE_TYPE{FASTA, FASTQ};
//enum SEQ_FILE_MODE{TEXT, GZIP, BZIP2};

class SeqIO{
protected:
  uint32_t fileidx;//index of current opened file.
  ifstream filehandler;//file stream
  vector<SEQ_FILE_TYPE> filetypes;
public:
  vector<string> seqfiles;
  // The sequence files are lines in the file <fname>
  SeqIO(string fname){
    ifstream fin;
    fin.open(fname);
    if (not fin.is_open()){
      MSG::error("Failed to open file "+fname);
    }
    string line;
    while(std::getline(fin, line)){
      seqfiles.push_back(line);
    }
    fin.close();
    // inferFileType();
  }
  SeqIO(vector<string> files):seqfiles(files){
    // inferFileType();
  }
  // virtual void inferFileType(){ MSG::message("SeqIO inferFileType"); };
  void start(){
    if(seqfiles.size() == 0){
      MSG::error("No sequence files.");
    }
    fileidx=0;
    filehandler.open(seqfiles[fileidx]);
    if(not filehandler.is_open()){
      MSG::error("Failed to open file "+seqfiles[fileidx]);
    }
  }
  void close(){
    if(filehandler.is_open()){
      filehandler.close();
    }
  }
  void reset(){
    if(filehandler.is_open()){
      filehandler.close();
    }
    start();
  }
  ~SeqIO(){
    if(filehandler.is_open()){
      filehandler.close();
    }
  }
};

class DNASeqIO : public SeqIO{
public:
  DNASeqIO(vector<string> files): SeqIO(files){ inferFileType(); }
  //return NULL when reaching end.
  DNASeq read_next();
  void inferFileType();
};

#endif
