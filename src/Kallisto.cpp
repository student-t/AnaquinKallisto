#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "KmerIterator.hpp"

// Kmers counting
static std::map<std::string, unsigned> __kcounts__;

template <typename Out> void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

/*
 * Build list of reference k-mers for sequins and human genome
 */

void Anaquin_buildRef()
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
        
        __kcounts__[toks[1]] = 0; // Normal
        __kcounts__[toks[2]] = 0; // Reverse complement
    }
    
    r.close();
}

static void Anaquin_countKM(const char *s1)
{
    KmerIterator iter(s1), end;
    
    for (int i = 0;  iter != end; ++i, ++iter)
    {
        const auto s = iter->first.rep().toString();
        
        if (__kcounts__.count(s))
        {
            __kcounts__[s]++;
        }
    }
}

void Anaquin_countKM(const char *s1, const char *s2)
{
    Anaquin_countKM(s1);
    Anaquin_countKM(s2);
}

void Anaquin_printResults()
{
    std::ofstream w("KallistoCount.txt");
    
    for (const auto &i : __kcounts__)
    {
        w << i.first << "\t" << i.second << std::endl;
    }
    
    w.close();
}


