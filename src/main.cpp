#include <string>
#include <iostream>
#include <random>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <getopt.h>
 #include <thread>
#include <time.h>

#include <cstdio>

#include <zlib.h>

#include "common.h"
#include "ProcessReads.h"
#include "KmerIndex.h"
#include "Kmer.hpp"
#include "MinCollector.h"
#include "weights.h"
#include "Inspect.h"
#include "H5Writer.h"


//#define ERROR_STR "\033[1mError:\033[0m"
#define ERROR_STR "Error:"

using namespace std;


void ParseOptionsIndex(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int make_unique_flag = 0;
  const char *opt_string = "i:k:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"make-unique", no_argument, &make_unique_flag, 1},
    // short args
    {"index", required_argument, 0, 'i'},
    {"kmer-size", required_argument, 0, 'k'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 'k': {
      stringstream(optarg) >> opt.k;
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }
  if (make_unique_flag) {
    opt.make_unique = true;
  }

  for (int i = optind; i < argc; i++) {
    opt.transfasta.push_back(argv[i]);
  }
}

void ParseOptionsInspect(int argc, char **argv, ProgramOptions& opt) {


  const char *opt_string = "";
  static struct option long_options[] = {
    // long args
    {"gfa", required_argument, 0, 'g'},
    {0,0,0,0}
  };

  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 'g': {
      opt.gfa = optarg;
      break;
    }
    default: break;
    }
  }
  opt.index = argv[optind];
}



void ParseOptionsEM(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int plaintext_flag = 0;
  int write_index_flag = 0;
  int single_flag = 0;
  int strand_FR_flag = 0;
  int strand_RF_flag = 0;
  int bias_flag = 0;
  int pbam_flag = 0;
  int fusion_flag = 0;

  const char *opt_string = "t:i:l:s:o:n:m:d:b:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"plaintext", no_argument, &plaintext_flag, 1},
    {"write-index", no_argument, &write_index_flag, 1},
    {"single", no_argument, &single_flag, 1},
    {"fr-stranded", no_argument, &strand_FR_flag, 1},
    {"rf-stranded", no_argument, &strand_RF_flag, 1},
    {"bias", no_argument, &bias_flag, 1},
    {"pseudobam", no_argument, &pbam_flag, 1},
    {"fusion", no_argument, &fusion_flag, 1},
    {"seed", required_argument, 0, 'd'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"index", required_argument, 0, 'i'},
    {"fragment-length", required_argument, 0, 'l'},
    {"sd", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {"iterations", required_argument, 0, 'n'},
    {"min-range", required_argument, 0, 'm'},
    {"bootstrap-samples", required_argument, 0, 'b'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 'l': {
      stringstream(optarg) >> opt.fld;
      break;
    }
    case 's': {
      stringstream(optarg) >> opt.sd;
      break;
    }
    case 'o': {
      opt.output = optarg;
      break;
    }
    case 'n': {
      stringstream(optarg) >> opt.iterations;
      break;
    }
    case 'm': {
      stringstream(optarg) >> opt.min_range;
    }
    case 'b': {
      stringstream(optarg) >> opt.bootstrap;
      break;
    }
    case 'd': {
      stringstream(optarg) >> opt.seed;
      break;
    }
    default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }

  if (verbose_flag) {
    opt.verbose = true;
  }

  if (plaintext_flag) {
    opt.plaintext = true;
  }

  if (write_index_flag) {
    opt.write_index = true;
  }

  if (single_flag) {
    opt.single_end = true;
  }

  if (strand_FR_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::FR;
  }
  
  if (strand_RF_flag) {
    opt.strand_specific = true;
    opt.strand = ProgramOptions::StrandType::RF;
  }

  if (bias_flag) {
    opt.bias = true;
  }

  if (pbam_flag) {
    opt.pseudobam = true;
  }

  if (fusion_flag) {
    opt.fusion = true;
  }
}

void ParseOptionsEMOnly(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int plaintext_flag = 0;

  const char *opt_string = "t:s:l:s:o:n:m:d:b:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"plaintext", no_argument, &plaintext_flag, 1},
    {"seed", required_argument, 0, 'd'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"fragment-length", required_argument, 0, 'l'},
    {"sd", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {"iterations", required_argument, 0, 'n'},
    {"min-range", required_argument, 0, 'm'},
    {"bootstrap-samples", required_argument, 0, 'b'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'l': {
      stringstream(optarg) >> opt.fld;
      break;
    }
    case 's': {
      stringstream(optarg) >> opt.sd;
      break;
    }
    case 'o': {
      opt.output = optarg;
      break;
    }
    case 'n': {
      stringstream(optarg) >> opt.iterations;
      break;
    }
    case 'm': {
      stringstream(optarg) >> opt.min_range;
    }
    case 'b': {
      stringstream(optarg) >> opt.bootstrap;
      break;
    }
    case 'd': {
      stringstream(optarg) >> opt.seed;
      break;
    }
    default: break;
    }
  }

  if (verbose_flag) {
    opt.verbose = true;
  }

  if (plaintext_flag) {
    opt.plaintext = true;
  }
}

void ParseOptionsPseudo(int argc, char **argv, ProgramOptions& opt) {
  int verbose_flag = 0;
  int single_flag = 0;
  int strand_flag = 0;
  int pbam_flag = 0;
  int umi_flag = 0;

  const char *opt_string = "t:i:l:s:o:b:";
  static struct option long_options[] = {
    // long args
    {"verbose", no_argument, &verbose_flag, 1},
    {"single", no_argument, &single_flag, 1},
    //{"strand-specific", no_argument, &strand_flag, 1},
    {"pseudobam", no_argument, &pbam_flag, 1},
    {"umi", no_argument, &umi_flag, 'u'},
    {"batch", required_argument, 0, 'b'},
    // short args
    {"threads", required_argument, 0, 't'},
    {"index", required_argument, 0, 'i'},
    {"fragment-length", required_argument, 0, 'l'},
    {"sd", required_argument, 0, 's'},
    {"output-dir", required_argument, 0, 'o'},
    {0,0,0,0}
  };
  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 't': {
      stringstream(optarg) >> opt.threads;
      break;
    }
    case 'i': {
      opt.index = optarg;
      break;
    }
    case 'l': {
      stringstream(optarg) >> opt.fld;
      break;
    }
    case 's': {
      stringstream(optarg) >> opt.sd;
      break;
    }
    case 'o': {
      opt.output = optarg;
      break;
    }
    case 'b': {
      opt.batch_mode = true;
      opt.batch_file_name = optarg;
      break;
    }
    default: break;
    }
  }
  
  if (umi_flag) {
    opt.umi = true;
    opt.single_end = true; // UMI implies single end reads
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
 
  if (verbose_flag) {
    opt.verbose = true;
  }

  if (single_flag) {
    opt.single_end = true;
  }

  if (strand_flag) {
    opt.strand_specific = true;
  }

  if (pbam_flag) {
    opt.pseudobam = true;
  }
  
  
}


void ParseOptionsH5Dump(int argc, char **argv, ProgramOptions& opt) {
  int peek_flag = 0;
  const char *opt_string = "o:";
  static struct option long_options[] = {
    // long args
    {"peek", no_argument, &peek_flag, 1},
    // short args
    {"output-dir", required_argument, 0, 'o'},
    {0,0,0,0}
  };

  int c;
  int option_index = 0;
  while (true) {
    c = getopt_long(argc, argv, opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0:
      break;
    case 'o': {
      opt.output = optarg;
      break;
    }
    default: break;
    }
  }

  // there should only be one thing here...
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }

  if (peek_flag) {
    opt.peek = true;
  }
}

bool CheckOptionsIndex(ProgramOptions& opt) {

  bool ret = true;

  if (opt.k <= 1 || opt.k >= Kmer::MAX_K) {
    cerr << "Error: invalid k-mer length " << opt.k << ", minimum is 3 and  maximum is " << (Kmer::MAX_K -1) << endl;
    ret = false;
  }

  if (opt.k % 2 == 0) {
    cerr << "Error: k needs to be an odd number" << endl;
    ret = false;
  }

  if (opt.transfasta.empty()) {
    cerr << "Error: no FASTA files specified" << endl;
    ret = false;
  } else {

    for (auto& fasta : opt.transfasta) {
      // we want to generate the index, check k, index and transfasta
      struct stat stFileInfo;
      auto intStat = stat(fasta.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: FASTA file not found " << fasta << endl;
        ret = false;
      }
    }
  }

  if (opt.index.empty()) {
    cerr << "Error: need to specify kallisto index name" << endl;
    ret = false;
  }

  return ret;
}

bool CheckOptionsEM(ProgramOptions& opt, bool emonly = false) {

  bool ret = true;

  cerr << endl;
  // check for index
  if (!emonly) {
    if (opt.index.empty()) {
      cerr << ERROR_STR << " kallisto index file missing" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      auto intStat = stat(opt.index.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << ERROR_STR << " kallisto index file not found " << opt.index << endl;
        ret = false;
      }
    }
  }

  // check for read files
  if (!emonly) {
    if (opt.files.size() == 0) {
      cerr << ERROR_STR << " Missing read files" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      for (auto& fn : opt.files) {
        auto intStat = stat(fn.c_str(), &stFileInfo);
        if (intStat != 0) {
          cerr << ERROR_STR << " file not found " << fn << endl;
          ret = false;
        }
      }
    }

    /*
    if (opt.strand_specific && !opt.single_end) {
      cerr << "Error: strand-specific mode requires single end mode" << endl;
      ret = false;
    }*/

    if (!opt.single_end) {
      if (opt.files.size() % 2 != 0) {
        cerr << "Error: paired-end mode requires an even number of input files" << endl
             << "       (use --single for processing single-end reads)" << endl;
        ret = false;
      }
    }
  }

  if ((opt.fld != 0.0 && opt.sd == 0.0) || (opt.sd != 0.0 && opt.fld == 0.0)) {
    cerr << "Error: cannot supply mean/sd without supplying both -l and -s" << endl;
    ret = false;
  }

  if (opt.single_end && (opt.fld == 0.0 || opt.sd == 0.0)) {
    cerr << "Error: fragment length mean and sd must be supplied for single-end reads using -l and -s" << endl;
    ret = false;
  } else if (opt.fld == 0.0 && ret) {
    // In the future, if we have single-end data we should require this
    // argument
    cerr << "[quant] fragment length distribution will be estimated from the data" << endl;
  } else if (ret && opt.fld > 0.0 && opt.sd > 0.0) {
    cerr << "[quant] fragment length distribution is truncated gaussian with mean = " <<
      opt.fld << ", sd = " << opt.sd << endl;
  }

  if (!opt.single_end && (opt.fld > 0.0 && opt.sd > 0.0)) {
    cerr << "[~warn] you specified using a gaussian but have paired end data" << endl;
    cerr << "[~warn] we suggest omitting these parameters and let us estimate the distribution from data" << endl;
  }

  if (opt.fld < 0.0) {
    cerr << "Error: invalid value for mean fragment length " << opt.fld << endl;
    ret = false;
  }

  if (opt.sd < 0.0) {
    cerr << "Error: invalid value for fragment length standard deviation " << opt.sd << endl;
    ret = false;
  }

  if (opt.iterations <= 0) {
    cerr << "Error: invalid number of iterations " << opt.iterations << endl;
    ret = false;
  }

  if (opt.min_range <= 0) {
    cerr << "Error: invalid value for minimum range " << opt.min_range << endl;
    ret = false;
  }
  

  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      } else if (emonly) {
        // check for directory/counts.txt
        struct stat stCountInfo;
        auto intcountstat = stat((opt.output + "/counts.txt" ).c_str(), &stCountInfo);
        if (intcountstat != 0) {
          cerr << "Error: could not find file " << opt.output << "/counts.txt" << endl;
          ret = false;
        }

        // check for directory/index.saved
        struct stat stIndexInfo;
        auto intindexstat = stat((opt.output + "/index.saved").c_str(), &stIndexInfo);
        if (intindexstat != 0) {
          cerr << "Error: could not find index " << opt.output << "/index.saved" << endl;
          ret = false;
        }
        opt.index = (opt.output + "/index.saved");
      }
    } else {
      if (emonly) {
        cerr << "Error: output directory needs to exist, run quant first" << endl;
        ret = false;
      } else {
        // create directory
        if (mkdir(opt.output.c_str(), 0777) == -1) {
          cerr << "Error: could not create directory " << opt.output << endl;
          ret = false;
        }
      }
    }
  }

  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "Warning: you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }
    if (opt.threads > 1 && opt.pseudobam) {
      cerr << "Error: pseudobam is not compatible with running on many threads."<< endl;
      ret = false;
    }
  }

  if (opt.bootstrap < 0) {
    cerr << "Error: number of bootstrap samples must be a non-negative integer." << endl;
    ret = false;
  }

  return ret;
}



bool CheckOptionsPseudo(ProgramOptions& opt) {

  bool ret = true;

  cerr << endl;
  // check for index
  if (opt.index.empty()) {
    cerr << ERROR_STR << " kallisto index file missing" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << ERROR_STR << " kallisto index file not found " << opt.index << endl;
      ret = false;
    }
  }

  // check for read files
  if (!opt.batch_mode) {
    if (opt.umi) {
      cerr << ERROR_STR << " UMI must be run in batch mode, use --batch option" << endl;
      ret = false;      
    }
    
    if (opt.files.size() == 0) {
      cerr << ERROR_STR << " Missing read files" << endl;
      ret = false;
    } else {
      struct stat stFileInfo;      
      for (auto& fn : opt.files) {        
        auto intStat = stat(fn.c_str(), &stFileInfo);
        if (intStat != 0) {
          cerr << ERROR_STR << " file not found " << fn << endl;
          ret = false;
        }
      }
    }
  } else {
    if (opt.files.size() != 0) {
      cerr << ERROR_STR << " cannot specify batch mode and supply read files" << endl;
      ret = false;
    } else {
      // check for batch files
      if (opt.batch_mode) {
        struct stat stFileInfo;
        auto intstat = stat(opt.batch_file_name.c_str(), &stFileInfo);
        if (intstat != 0) {
          cerr << ERROR_STR << " file not found " << opt.batch_file_name << endl;
          ret = false;
        }
        // open the file, parse and fill the batch_files values
        std::ifstream bfile(opt.batch_file_name);
        std::string line;
        std::string id,f1,f2;
        while (std::getline(bfile,line)) {
          if (line.size() == 0) {
            continue;
          }
          std::stringstream ss(line);
          ss >> id;
          if (id[0] == '#') {
            continue;
          }
          opt.batch_ids.push_back(id);
          if (opt.single_end && !opt.umi) {
            ss >> f1;
            opt.batch_files.push_back({f1});
            intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f1 << endl;
              ret = false;
            }
          } else {
            ss >> f1 >> f2;
            if (!opt.umi) {
              opt.batch_files.push_back({f1,f2});
            } else {
              opt.umi_files.push_back(f1);
              opt.batch_files.push_back({f2});
            }
            intstat = stat(f1.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f1 << endl;
              ret = false;
            }
            intstat = stat(f2.c_str(), &stFileInfo);
            if (intstat != 0) {
              cerr << ERROR_STR << " file not found " << f2 << endl;
              ret = false;
            }
          }
        }
      }
    }
  }


  /*
  if (opt.strand_specific && !opt.single_end) {
    cerr << "Error: strand-specific mode requires single end mode" << endl;
    ret = false;
  }*/

  if (!opt.single_end) {
    if (opt.files.size() % 2 != 0) {
      cerr << "Error: paired-end mode requires an even number of input files" << endl
           << "       (use --single for processing single-end reads)" << endl;
      ret = false;
    }
  }
  
  if (opt.umi) {
    opt.single_end = true;
    if (opt.fld != 0.0 || opt.sd != 0.0) {
      cerr << "[~warn] you supplied fragment length information for UMI data which will be ignored" << endl;
    }
  } else {
    if ((opt.fld != 0.0 && opt.sd == 0.0) || (opt.sd != 0.0 && opt.fld == 0.0)) {
      cerr << "Error: cannot supply mean/sd without supplying both -l and -s" << endl;
      ret = false;
    }

    if (opt.single_end && (opt.fld == 0.0 || opt.sd == 0.0)) {
      cerr << "Error: fragment length mean and sd must be supplied for single-end reads using -l and -s" << endl;
      ret = false;
    } else if (opt.fld == 0.0 && ret) {
      // In the future, if we have single-end data we should require this
      // argument
      cerr << "[quant] fragment length distribution will be estimated from the data" << endl;
    } else if (ret && opt.fld > 0.0 && opt.sd > 0.0) {
      cerr << "[quant] fragment length distribution is truncated gaussian with mean = " <<
        opt.fld << ", sd = " << opt.sd << endl;
    }

    if (!opt.single_end && (opt.fld > 0.0 && opt.sd > 0.0)) {
      cerr << "[~warn] you specified using a gaussian but have paired end data" << endl;
      cerr << "[~warn] we suggest omitting these parameters and let us estimate the distribution from data" << endl;
    }
  }

  if (opt.fld < 0.0) {
    cerr << "Error: invalid value for mean fragment length " << opt.fld << endl;
    ret = false;
  }

  if (opt.sd < 0.0) {
    cerr << "Error: invalid value for fragment length standard deviation " << opt.sd << endl;
    ret = false;
  }

  if (opt.output.empty()) {
    cerr << "Error: need to specify output directory " << opt.output << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.output.c_str(), &stFileInfo);
    if (intStat == 0) {
      // file/dir exits
      if (!S_ISDIR(stFileInfo.st_mode)) {
        cerr << "Error: file " << opt.output << " exists and is not a directory" << endl;
        ret = false;
      }
    } else {
      // create directory
      if (mkdir(opt.output.c_str(), 0777) == -1) {
        cerr << "Error: could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  }

  if (opt.threads <= 0) {
    cerr << "Error: invalid number of threads " << opt.threads << endl;
    ret = false;
  } else {
    unsigned int n = std::thread::hardware_concurrency();
    if (n != 0 && n < opt.threads) {
      cerr << "[~warn]  you asked for " << opt.threads
           << ", but only " << n << " cores on the machine" << endl;
    }
    if (opt.threads > 1 && opt.pseudobam) {
      cerr << "Error: pseudobam is not compatible with running on many threads."<< endl;
      ret = false;
    }
  }

  return ret;
}


bool CheckOptionsInspect(ProgramOptions& opt) {

  bool ret = true;
  // check for index
  if (opt.index.empty()) {
    cerr << "Error: kallisto index file missing" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    auto intStat = stat(opt.index.c_str(), &stFileInfo);
    if (intStat != 0) {
      cerr << "Error: kallisto index file not found " << opt.index << endl;
      ret = false;
    }
  }

  return ret;
}

bool CheckOptionsH5Dump(ProgramOptions& opt) {
  bool ret = true;
  if (!opt.peek) {
    if ( opt.output.size() == 0) {
      cerr << "Error: You must specify an output directory." << endl;
      ret = false;
    } else {
      struct stat stFileInfo;
      auto intStat = stat(opt.output.c_str(), &stFileInfo);
      if (intStat == 0) {
        // file/dir exits
        if (!S_ISDIR(stFileInfo.st_mode)) {
          cerr << "Error: tried to open " << opt.output << " but another file already exists there" << endl;
          ret = false;
        }
      } else if (mkdir(opt.output.c_str(), 0777) == -1) {
        cerr << "Error: could not create directory " << opt.output << endl;
        ret = false;
      }
    }
  } else {
    if (opt.output.size() > 0) {
      cerr << "Error: Cannot specify output directory and '--peek'. Please specify only one." << endl;
      ret = false;
    }
  }

  if (opt.files.size() == 0) {
    cerr << "Error: Missing H5 files" << endl;
    ret = false;
  } else if (opt.files.size() > 1) {
    cerr << "Error: Please specify only one H5 file" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    for (auto& fn : opt.files) {
      auto intStat = stat(fn.c_str(), &stFileInfo);
      if (intStat != 0) {
        cerr << "Error: H5 file not found " << fn << endl;
        ret = false;
      }
    }
  }

  return ret;
}

void PrintCite() {
  cout << "When using this program in your research, please cite" << endl << endl
       << "  Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L." << endl
       << "  Near-optimal probabilistic RNA-seq quantification, "<< endl
       << "  Nature Biotechnology 34, 525-527(2016), doi:10.1038/nbt.3519" << endl
       << endl;
}

void PrintVersion() {
  cout << "kallisto, version " << 	KALLISTO_VERSION << endl;
}

void usage() {
  cout << "kallisto " << KALLISTO_VERSION << endl << endl
       << "Usage: kallisto <CMD> [arguments] .." << endl << endl
       << "Where <CMD> can be one of:" << endl << endl
       << "    index         Builds a kallisto index "<< endl
       << "    quant         Runs the quantification algorithm " << endl
       << "    pseudo        Runs the pseudoalignment step " << endl
       << "    h5dump        Converts HDF5-formatted results to plaintext" << endl
       << "    version       Prints version information"<< endl
       << "    cite          Prints citation information" << endl << endl
       << "Running kallisto <CMD> without arguments prints usage information for <CMD>"<< endl << endl;
}


void usageIndex() {
  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Builds a kallisto index" << endl << endl
       << "Usage: kallisto index [arguments] FASTA-files" << endl << endl
       << "Required argument:" << endl
       << "-i, --index=STRING          Filename for the kallisto index to be constructed " << endl << endl
       << "Optional argument:" << endl
       << "-k, --kmer-size=INT         k-mer (odd) length (default: 31, max value: " << (Kmer::MAX_K-1) << ")" << endl
       << "    --make-unique           Replace repeated target names with unique names" << endl
       << endl;

}

void usageh5dump() {
  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Converts HDF5-formatted results to plaintext" << endl << endl
       << "Usage:  kallisto h5dump [arguments] abundance.h5" << endl << endl
       << "Required argument:" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl;
}

void usageInspect() {
  cout << "kallisto " << KALLISTO_VERSION << endl << endl
       << "Usage: kallisto inspect INDEX-file" << endl << endl
       << "Optional arguments:" << endl
       << "    --gfa=STRING              Filename for GFA output of T-DBG" << endl << endl;
}

void usageEM(bool valid_input = true) {
  if (valid_input) {

  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Computes equivalence classes for reads and quantifies abundances" << endl << endl;
  }
  //      "----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|"
  cout << "Usage: kallisto quant [arguments] FASTQ-files" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << endl
       << "                              quantification" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "    --bias                    Perform sequence based bias correction" << endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)" << endl
       << "    --seed=INT                Seed for the bootstrap sampling (default: 42)" << endl
       << "    --plaintext               Output plaintext instead of HDF5" << endl
       << "    --fusion                  Search for fusions for Pizzly" << endl
       << "    --single                  Quantify single-end reads" << endl
       << "    --fr-stranded             Strand specific reads, first read forward" << endl
       << "    --rf-stranded             Strand specific reads, first read reverse" << endl
       << "-l, --fragment-length=DOUBLE  Estimated average fragment length" << endl
       << "-s, --sd=DOUBLE               Estimated standard deviation of fragment length" << endl
       << "                              (default: -l, -s values are estimated from paired" << endl
       << "                               end data, but are required when using --single)" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl
       << "    --pseudobam               Output pseudoalignments in SAM format to stdout" << endl;

}

void usagePseudo(bool valid_input = true) {
  if (valid_input) {
    cout << "kallisto " << KALLISTO_VERSION << endl
         << "Computes equivalence classes for reads and quantifies abundances" << endl << endl;
  }

  cout << "Usage: kallisto pseudo [arguments] FASTQ-files" << endl << endl
       << "Required arguments:" << endl
       << "-i, --index=STRING            Filename for the kallisto index to be used for" << endl
       << "                              pseudoalignment" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "-u  --umi                     First file in pair is a UMI file" << endl
       << "-b  --batch=FILE              Process files listed in FILE" << endl
       << "    --single                  Quantify single-end reads" << endl
       << "-l, --fragment-length=DOUBLE  Estimated average fragment length" << endl
       << "-s, --sd=DOUBLE               Estimated standard deviation of fragment length" << endl
       << "                              (default: -l, -s values are estimated from paired" << endl
       << "                               end data, but are required when using --single)" << endl
       << "-t, --threads=INT             Number of threads to use (default: 1)" << endl
       << "    --pseudobam               Output pseudoalignments in SAM format to stdout" << endl;

}

void usageEMOnly() {
  cout << "kallisto " << KALLISTO_VERSION << endl
       << "Computes equivalence classes for reads and quantifies abundance" << endl << endl
       << "Usage: kallisto quant-only [arguments]" << endl << endl
       << "Required argument:" << endl
       << "-o, --output-dir=STRING       Directory to write output to" << endl << endl
       << "Optional arguments:" << endl
       << "-l, --fragment-length=DOUBLE  Estimated fragment length (default: value is estimated from the input data)" << endl
       << "-b, --bootstrap-samples=INT   Number of bootstrap samples (default: 0)" << endl
       << "    --seed=INT                Seed for the bootstrap sampling (default: 42)" << endl
       << "    --plaintext               Output plaintext instead of HDF5" << endl << endl;
}

std::string argv_to_string(int argc, char *argv[]) {
  std::string res;
  for (int i = 0; i < argc; ++i) {
    res += argv[i];
    if (i + 1 < argc) {
      res += " ";
    }
  }

  return res;
}

std::string get_local_time() {
  time_t rawtime;
  struct tm * timeinfo;

  time( &rawtime );
  timeinfo = localtime( &rawtime );
  std::string ret(asctime(timeinfo));

  // chomp off the newline
  return ret.substr(0, ret.size() - 1);
}

int main(int argc, char *argv[])
{
    ProgramOptions opt;
    string cmd(argv[1]);
   if (cmd == "quant") {
      if (argc==2) {
        usageEM();
        return 0;
      }
      ParseOptionsEM(argc-1,argv+1,opt);
      if (!CheckOptionsEM(opt)) {
        cerr << endl;
        usageEM(false);
        exit(1);
      } else {
        KmerIndex index(opt);
          
          extern void KMInit();
          KMInit();
          
        MinCollector collection(index, opt);
        int num_processed = 0;
        num_processed = ProcessReads(index, opt, collection);
          
          extern void KMAll();
          extern void KMResults();
          
          KMAll();
          //KMResults();
          
          num_processed = num_processed;
      }
    }
   
  return 0;
}
