/*
 * ============================================================================
 *
 *        Authors: Christina SHI <Christina.hshi@gmail.com>
 *
 * ============================================================================
 */

#pragma once

#ifndef PARAMS_H
#define PARAMS_H

#include "base.h"

class Params_generic{
public:
  size_t thread_num;
  //extra thread to open for running the program
  Params_generic(int th=0): thread_num(th){}
  friend ostream& operator<<( ostream &output, const Params_generic& params);
};

ostream& operator<<( ostream &output, const Params_generic& params){
  output<<"Params:"<<endl
  // <<"\tkmer size:                "<<params.K<<endl
  // <<"\tkmer min. abundance:      "<<params.kmer_abundance_min<<endl
  // <<"\tsolid kmer min. abundance:"<<params.solid_kmer_abundance_min<<endl
  // <<"\tsolid kmer max. abundance:"<<params.solid_kmer_abundance_max<<endl
  <<"\tthreads:                  "<<params.thread_num<<endl;
  return output;
}

class Params_insilicoDigest : public Params_generic{
public:
  vector<Enzyme> enzymes;
  vector<string> seqfiles;
  string mapfile;//digested sequences will be saved in <mapfile> in CMAP format.
  uint32_t mismatch_max; // maximum allowed mismatches when finding
  uint32_t merge_dis;
  Params_insilicoDigest(){}
  Params_insilicoDigest(vector<Enzyme> ens, vector<string> seqs, string outfile,
    uint32_t mis, uint32_t merge, size_t th):
      enzymes(ens), seqfiles(seqs), mapfile(outfile),
      mismatch_max(mis), merge_dis(merge), Params_generic(th){}
  static Params_insilicoDigest parse_opts(const string& program, const string& cmd,  vector<string> opts_remained){
    namespace po = boost::program_options;
    po::options_description desc(program+" "+cmd+" <options>\nOptions");
    desc.add_options()
        ("help,h", "print help messages")
        ("enzymes,e", po::value<string>()->required(), "file containing enzyme info with one enzyme per line in TSV format: <Name> <Seq> <Digest_site> <Mismatch_max> <Type(1:Restriction Enzyme 2:CRISPR/Cas9)> <Channel>")
        ("seqfiles,f", po::value<string>()->required(), "file containing sequence files with one file per line.")
        ("mismatch,m", po::value<int>()->default_value(0), "maximum number of mismatches allowed when searching for targets during in silico digestion. [Deprecated! please set for each enzyme in enzyme file.]")
        ("mergeDis,d", po::value<int>()->default_value(1000), "merge nearby labels <= <mergeDis>")
        ("output,o", po::value<string>()->default_value("digested.cmap"), "file to save the digested map in CMAP format")
        ("thread,t", po::value<int>()->default_value(1), "number of threads to run program.");

    po::variables_map vm;
    po::store(po::command_line_parser(opts_remained).options(desc).run(), vm);
    if(opts_remained.size()==0 or vm.count("help")){
      std::cerr<<std::endl<<desc<<std::endl;
      exit(0);
    }
    po::notify(vm);

    string enzymefile = vm["enzymes"].as<string>();
    vector<Enzyme> enzymes = Enzyme::load_from_tsv(enzymefile);
    string seqfilesfile = vm["seqfiles"].as<string>();
    vector<string> seqfiles;
    ifstream fin;
    string line;
    fin.open(seqfilesfile, ios::in);
    while(getline(fin, line)){
      seqfiles.push_back(line);
    }
    fin.close();

    Params_insilicoDigest params;
    params.enzymes = enzymes;
    params.seqfiles = seqfiles;
    params.mapfile = vm["output"].as<string>();
    params.mismatch_max = vm["mismatch"].as<int>();
    params.merge_dis = vm["mergeDis"].as<int>();
    params.thread_num = vm["thread"].as<int>();
    return params;
  }
};

class Params_find20merNGG : public Params_generic{
public:
  vector<string> seqfiles;
  string output;
  bool skiploc;
  Params_find20merNGG(){}
  Params_find20merNGG(vector<string> seqs, string out, bool skip): seqfiles(seqs), output(out), skiploc(skip){}
  static Params_find20merNGG parse_opts(const string& program, const string& cmd, vector<string> opts_remained){
    namespace po = boost::program_options;
    po::options_description desc(program+" "+cmd+" <options>\nOptions");
    desc.add_options()
        ("help,h", "print help messages.")
        ("seqfiles,f", po::value<string>()->required(), "file containing sequence files with one file per line.")
        ("skiploc", po::bool_switch()->default_value(false), "not output the loci of 20merNGG (only counts) to save space.")
        ("output,o", po::value<string>()->default_value("20merNGG.tsv"), "output file.");

    po::variables_map vm;
    po::store(po::command_line_parser(opts_remained).options(desc).run(), vm);
    if(opts_remained.size()==0 or vm.count("help")){
      std::cerr<<std::endl<<desc<<std::endl;
      exit(0);
    }
    po::notify(vm);

    string seqfilesfile = vm["seqfiles"].as<string>();
    vector<string> seqfiles;
    ifstream fin;
    string line;
    fin.open(seqfilesfile, ios::in);
    while(getline(fin, line)){
      seqfiles.push_back(line);
    }
    fin.close();

    Params_find20merNGG params;
    params.seqfiles = seqfiles;
    params.output = vm["output"].as<string>();
    params.skiploc = vm.count("skiploc");
    return params;
  }
};

class Params_analyze : public Params_generic{
public:
  vector<string> seqfiles;
  string cmapfile;
  string kmerfile;//20merNGG file.
  int gapsize_min;

  Params_analyze(){}
  Params_analyze(vector<string> seqs): seqfiles(seqs){}
  static Params_analyze parse_opts(const string& program, const string& cmd, vector<string> opts_remained){
    namespace po = boost::program_options;
    po::options_description desc(program+" "+cmd+" <options>\nOptions");
    desc.add_options()
        ("help,h", "print help messages.")
        ("seqfiles,f", po::value<string>()->required(), "file containing sequence files with one file per line.")
        ("cmapfile,c", po::value<string>()->required(), "cmap file showing existing labels.")
        ("kmerfile,m", po::value<string>()->required(), "20merNGG file.")
        ("gapsize_min,g", po::value<int>()->default_value(20000), "minimum size of gap to fill by singleton gRNAs.");

    po::variables_map vm;
    po::store(po::command_line_parser(opts_remained).options(desc).run(), vm);
    if(opts_remained.size()==0 or vm.count("help")){
      std::cerr<<std::endl<<desc<<std::endl;
      exit(0);
    }
    po::notify(vm);

    string seqfilesfile = vm["seqfiles"].as<string>();
    vector<string> seqfiles;
    ifstream fin;
    string line;
    fin.open(seqfilesfile, ios::in);
    while(getline(fin, line)){
      seqfiles.push_back(line);
    }
    fin.close();

    Params_analyze params;
    params.seqfiles = seqfiles;
    params.cmapfile = vm["cmapfile"].as<string>();
    params.kmerfile = vm["kmerfile"].as<string>();
    params.gapsize_min = vm["gapsize_min"].as<int>();
    return params;
  }
};

typedef std::variant<Params_insilicoDigest, Params_find20merNGG, Params_analyze> Params;

// The subcommand idea credit to https://gist.github.com/randomphrase/10801888
Params parse_opts(int argc, const char *argv[], string& cmd){
  namespace po = boost::program_options;

  string program = string(argv[0]);
  po::options_description global(program+" <command> <options>\nUsage");
  global.add_options()
      ("help,h", "print help messages.")
      ("command", po::value<string>()->required(), "pick one command to execute:\n  (1)insilicoDigest\n  (2)find20merNGG")
      ("options", po::value<vector<string>>(), "options for command");

  string helpinfo = "\n" +
    program+" <command> <options>\n"+
    "<command> available: \n"+
    "  insilicoDigest        In silico digestion of DNA sequences.\n"+
    "  find20merNGG          Find all 20merNGG in DNA.\n"+
    "  analyze               For analysis only.\n"
    "#Specify <command> to see <options>.\n";

  po::positional_options_description pos;
  pos.add("command", 1).
      add("options", -1);

  po::variables_map vm;
  po::parsed_options parsed = po::command_line_parser(argc, argv).
        options(global).
        positional(pos).
        allow_unregistered().
        run();
  po::store(parsed, vm);
  if(not vm.count("command")){
    std::cerr<<helpinfo<<endl;
    exit(0);
  }

  cmd = vm["command"].as<std::string>();
  // Collect all the unrecognized options from the first pass. This will include the
  // (positional) command name, so we need to erase that.
  std::vector<std::string> opts = po::collect_unrecognized(parsed.options, po::include_positional);
  opts.erase(opts.begin());

  if(cmd == "insilicoDigest"){
    return Params_insilicoDigest::parse_opts(program, cmd, opts);
  }else if(cmd == "find20merNGG"){
    return Params_find20merNGG::parse_opts(program, cmd, opts);
  }else if(cmd == "analyze"){
    return Params_analyze::parse_opts(program, cmd, opts);
  }else{
    std::cerr<<"[Error] Command "<<cmd<<" unrecognized!"<<endl;
    // std::cerr<<std::endl<<global<<std::endl;
    std::cerr<<helpinfo<<endl;
    exit(0);
  }
}

#endif
