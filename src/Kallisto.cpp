#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "KmerIterator.hpp"

// Reference kmers counting for sequins
static std::map<std::string, unsigned> __seqs__;

// All k-mers (for debugging)
static std::map<std::string, unsigned> __all__;

template <typename Out> void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

void KMInit()
{
    std::ifstream r("CancerKM.txt");
    
    if (!r.good())
    {
        throw std::runtime_error("Reference file for k-mers is missing");
    }
    
    std::string line;
    while (std::getline(r, line))
    {
        std::vector<std::string> toks;
        split(line, '\t', std::back_inserter(toks));
        assert(!toks.empty());
        
        if (toks[0] == "Name")
        {
            continue;
        }
        
        __seqs__[toks[1]] = 0; // Normal
        __seqs__[toks[2]] = 0; // Reverse complement
    }
    
    r.close();
}

static void KMCount(const char *s1)
{
    Kmer::k = 31;
    KmerIterator iter(s1), end;
    
    for (int i = 0; iter != end; ++i, ++iter)
    {
        const auto s = iter->first.rep().toString();
        
        if (__seqs__.count(s))
        {
            __seqs__[s]++;
        }
        
#ifdef DEBUG
        __all__[s]++;
#endif
    }
}

void KMCount(const char *s1, const char *s2)
{
    KMCount(s1);
    KMCount(s2);
}

void KMAll()
{
#ifdef DEBUG
    std::ofstream w("KMAll.txt");
    
    for (const auto &i : __all__)
    {
        w << i.first << "\t" << i.second << std::endl;
    }
    
    w.close();
#endif
}

void KMResults()
{
    std::ofstream w("KallistoCount.txt");
    
    for (const auto &i : __seqs__)
    {
        w << i.first << "\t" << i.second << std::endl;
    }
    
    w.close();
}
