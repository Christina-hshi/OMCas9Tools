/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * ============================================================================
 */

#include "insilico_digest.h"
#include "find_20merNGG.h"

void analyze(Params_analyze& params);

int main(int argc, const char* argv[]){
  string cmd;
  Params params = parse_opts(argc, argv, cmd);
  if(cmd == "insilicoDigest"){
    MSG::message(currentDateTime()+" Start in silico digestion...");
    insilico_digest(std::get<Params_insilicoDigest>(params));
    MSG::message(currentDateTime()+" Finished digestion.");
  }else if(cmd == "find20merNGG"){
    MSG::message(currentDateTime()+" Start finding 20merNGG...");
    find_20merNGG(std::get<Params_find20merNGG>(params));
    MSG::message(currentDateTime()+" Finished finding.");
  }else if(cmd == "analyze"){
    MSG::message(currentDateTime()+" Start analyzing...");
    analyze(std::get<Params_analyze>(params));
    MSG::message(currentDateTime()+" Finished analysis.");
  }else{
    MSG::error(cmd+" unsupported.");
  }
  return 0;
}

void analyze(Params_analyze& params){
  // Load reference sequences
  map<string, string> refname2seq;
  DNASeqIO seqIn(params.seqfiles);
  seqIn.start();
  DNASeq dnaseq = seqIn.read_next();
  while(not dnaseq.is_empty()){
    refname2seq[dnaseq.name] = dnaseq.seq;
    dnaseq = seqIn.read_next();
  }
  vector<size_t> refseqLens;
  // size_t total=0;
  for(auto ref : refname2seq){
    refseqLens.push_back(ref.second.length());
    // cout<<ref.first<<" "<<ref.second.length()<<endl;
    // total += ref.second.length();
  }
  size_t refseqLenTotal = std::accumulate(refseqLens.begin(), refseqLens.end(), (size_t)0);
  cout<<"[Info] "<<refname2seq.size()<<" reference sequences with "<<refseqLenTotal<<" bp in total."<<endl;
  // cout<<total<<endl;

  // Load only singleton gRNAs
  unordered_map<kmer_t, vector<GLoc>> kmer2loci;
  const int k=20;
  load(kmer2loci, params.kmerfile, k, false, 1, 1);

  map<int, string> cMapId2name;
  for(int x = 1; x <= 22; x++){
    cMapId2name[x] = "chr"+std::to_string(x);
  }
  cMapId2name[23] = "chrX";
  cMapId2name[24] = "chrY";
  cMapId2name[25] = "chrM";

  // Load cmap info
  ifstream fin;
	fin.open(params.cmapfile, ios::in);
	if(!fin.is_open()){
		std::cerr<<"[Error] failed to open "<<params.cmapfile<<endl;
	}
  map<int, vector<int>> cMap2labels;
  map<int, int> cMap2Lens;
	int cMapId, contigLength, numSites, siteID, labelChannel, position;
  string line;
	while(getline(fin, line) and line[0]=='#') continue;
	do{
		stringstream ss(line);
		ss>>cMapId>>contigLength>>numSites>>siteID>>labelChannel>>position;
    cMap2Lens[cMapId] = contigLength;
    //start form the first line
    vector<int> labels;
    for(int x = 0; x < numSites; x++){
      labels.push_back(position);
      getline(fin, line);
      stringstream ss_tmp(line);
      ss_tmp>>cMapId>>contigLength>>numSites>>siteID>>labelChannel>>position;
    }
    cMap2labels[cMapId] = labels;
    // cout<<cMapId<<" "<<labels.size()<<endl;
  }while(getline(fin, line));
  fin.close();

  cout<<"GapSizeMin,GapNum,GapLenTotal,GapMin,GapMax,GapAve,GapFilledNum,GapFilledLen,gRNAnum"<<endl;
  for(int gapsize_min = 10000; gapsize_min <= 50000; gapsize_min+=10000){
    params.gapsize_min = gapsize_min;
    cout<<gapsize_min<<",";
  //check gaps
  vector<GRange> gaps;
  for (auto& cMapInfo : cMap2labels){
    auto cMapId = cMapInfo.first;
    auto cMapName = cMapId2name[cMapId];
    auto cMapLen = cMap2Lens[cMapId];
    auto labels = cMapInfo.second;
    if(labels[0] >= params.gapsize_min){
      gaps.push_back(GRange(cMapName, 0, labels[0]));
    }
    for(int x = 1; x < labels.size(); x++){
      if(labels[x] - labels[x-1] >= params.gapsize_min){
        gaps.push_back(GRange(cMapName, labels[x-1], labels[x]));
      }
    }
    if(cMapLen - labels.back() >= params.gapsize_min){
      gaps.push_back(GRange(cMapName, labels.back(), cMapLen));
    }
  }
  if(false){ // write the gap info to file
    ofstream fout;
    fout.open("gapAtLeast"+std::to_string(params.gapsize_min)+".csv", std::ios_base::out);
    if(not fout.is_open()){
      MSG::error("Failed to open file.");
    }
    fout<<"Chr,Start,End"<<endl;
    for(auto gap : gaps){
      fout<<gap.ref<<","<<gap.start<<","<<gap.end<<endl;
    }
    fout.close();
  }
  // statistics of the gaps
  if(true){
    int gap_num = gaps.size();
    vector<int> gapLens;
    for(auto gap : gaps){
      gapLens.push_back(gap.end - gap.start);
    }
    int total_len = std::accumulate(gapLens.begin(), gapLens.end(), 0);
    int min = *std::min_element(gapLens.begin(), gapLens.end());
    int max = *std::max_element(gapLens.begin(), gapLens.end());
    cout<<gap_num<<","<<total_len<<","<<min<<","<<max<<","<<total_len/gap_num<<",";
    // cout<<"[Info] "<<gap_num<<" gaps with "<<total_len<<" bp in total. "
    //   <<"Min:"<<min<<"  Max:"<<max<<"  Ave:"<<total_len/gap_num<<endl;
  }

  // For each gap, identify all unique gRNAs inside,
  // and see whether we add smallest of them, such that all neighboring label distances are < gapsize_min
  vector<pair<GRange, int>> gaps_filled;
  for(auto gap : gaps){
    gap.seq = refname2seq[gap.ref].substr(gap.start, gap.end - gap.start);
    unordered_map<kmer_t, vector<GLoc>> reg_kmer2loci;
    find_20merNGG(gap, reg_kmer2loci);
    // gap.seq = "";
    // use only singletons
    vector<int> singles;
    for(auto kmer : reg_kmer2loci){
      if(kmer2loci.find(kmer.first) != kmer2loci.end()){
        singles.push_back(kmer.second[0].loc);
      }
    }
    if(singles.size() == 0) continue;
    std::sort(singles.begin(), singles.end());
    // check how many singletons are needed to fill the gaps
    // BFS
    bool reachEnd = false;
    int step = 1;
    set<int> nexts;
    for(int x = 0; x < singles.size(); x++){
      if(singles[x] - gap.start < params.gapsize_min){
        nexts.insert(x);
        if(gap.end - singles[x] < params.gapsize_min){
          reachEnd = true;
        }
      }
    }
    while( (not reachEnd) and (not nexts.empty()) ){
      step++;
      set<int> nexts_new;
      for(auto now : nexts){
        if(now == singles.size()-1) continue;
        for(int x = now+1; x < singles.size(); x++){
          if(singles[x] - singles[now] < params.gapsize_min){
            nexts_new.insert(x);
            if(gap.end - singles[x] < params.gapsize_min){
              reachEnd = true;
              break;
            }
          }
        }
      }
      nexts = nexts_new;
    }
    if(reachEnd){
      gaps_filled.push_back(pair<GRange, int>(gap, step));
    }
  }
  // filled gaps statistics
  vector<int> filledGapLens, steps;
  for(auto gap : gaps_filled){
    filledGapLens.push_back(gap.first.end - gap.first.start);
    steps.push_back(gap.second);
  }
  cout<<gaps_filled.size()<<","<<std::accumulate(filledGapLens.begin(), filledGapLens.end(), 0)<<","
    <<std::accumulate(steps.begin(), steps.end(), 0)<<endl;
  // cout<<"[Info] "<<gaps_filled.size()<<" gaps filled with "
  //   <<std::accumulate(filledGapLens.begin(), filledGapLens.end(), 0)<< " bp in total. "
  //   <<std::accumulate(steps.begin(), steps.end(), 0)<<" singleton gRNAs were used."<<endl;

  }//end loop
}
