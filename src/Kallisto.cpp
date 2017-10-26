#include <map>
#include <string>
#include <iostream>
#include <getopt.h>
#include <sys/stat.h>
#include "ProcessReads.h"

using namespace std;

void kallisto(const std::string &i, const std::string &p1, const std::string &p2)
{
    ProgramOptions opt;
    
    opt.index = i;
    opt.files.push_back(p1);
    opt.files.push_back(p2);

    KmerIndex index(opt);
    
    extern void KMInit();
    
    // Initalize for reference k-mers
    KMInit();
    
    MinCollector collection(index, opt);
    ProcessReads(index, opt, collection);
    
    extern void KMAll();
    extern void KMResults(unsigned &, unsigned &, std::map<std::string, unsigned> &);
    
    // Extract all k-mers (only for debugging)
    KMAll();

    unsigned nGen, nSeq;
    std::map<std::string, unsigned> span;
    
    KMResults(nGen, nSeq, span);
    assert((nGen + nSeq) > 0);
    
    std::cout << nGen << std::endl;
    std::cout << nSeq << std::endl;
    std::cout << (float)nSeq / (nSeq + nGen) << std::endl;
}

int main(int argc, char *argv[])
{
    kallisto(argv[1], argv[2], argv[3]);
    return 0;
}
