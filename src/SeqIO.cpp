/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * ============================================================================
 */

 #include "SeqIO.h"

 DNASeq DNASeqIO::read_next(){
   if(not filehandler.is_open()){
     return DNASeq();
   }
   string seqname, seq, quality="";
   if(filetypes[fileidx] == SEQ_FILE_TYPE::FASTA){
     string line;
     while(getline(filehandler, line) and line[0] != '>'){}
     size_t idx = 0;
     while(idx < line.length() and line[idx] != ' ' and line[idx] != '\t'){
       idx++;
     }
     seqname = line.substr(1, idx-1);
     seq="";
     while(getline(filehandler, line)){
       seq = seq + line;
       if(filehandler.peek() == '>'){
         break;
       }
     }
   }else if(filetypes[fileidx] == SEQ_FILE_TYPE::FASTQ){
     string line;
     while(getline(filehandler, line) and line[0] != '@'){}
     size_t idx = 0;
     while(idx < line.length() and line[idx] != ' ' and line[idx] != '\t'){
       idx++;
     }
     seqname = line.substr(1, idx-1);
     getline(filehandler, seq);
     getline(filehandler, line);
     getline(filehandler, quality);
   }else{
     MSG::error("Unsupported file type "+seqfiles[fileidx]);
   }
   if(filehandler.peek() == EOF){
     filehandler.close();
     fileidx++;
     if(fileidx < seqfiles.size()){
       filehandler.open(seqfiles[fileidx]);
       if(not filehandler.is_open()){
         MSG::error("Failed to open file "+seqfiles[fileidx]);
       }
     }
   }
   std::transform(seq.begin(), seq.end(),seq.begin(), ::toupper);
   return DNASeq(seqname, seq, quality);
 }

 void DNASeqIO::inferFileType(){
   MSG::message("DNASeqIO: infer read file types.");
   filetypes.resize(0);
   for(string fname : seqfiles){
       transform(fname.begin(), fname.end(), fname.begin(), ::toupper);
       if(fname.find(".FA") == fname.length()-3 or fname.find(".FASTA") == fname.length()-6){
         filetypes.push_back(SEQ_FILE_TYPE::FASTA);
         MSG::message("Please make sure no multiple line FASTA files.");
       }else if(fname.find(".FQ") == fname.length()-3 or fname.find(".FASTQ") == fname.length()-6){
         filetypes.push_back(SEQ_FILE_TYPE::FASTQ);
       }else{
         MSG::error("Failed to infer file type of "+fname);
       }
   }
 }
