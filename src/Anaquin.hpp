#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "KmerIndex.h"
#include "KmerIterator.hpp"

// All reference k-mers
static std::map<std::string, unsigned> __all__;

// Kallisto index for all sequin k-mers
static std::shared_ptr<KmerIndex> __seqsIndex__;

// Reference k-mers for spanning variants
static std::map<std::string, unsigned> __span__;

// All k-mers (for debugging)
static std::map<std::string, unsigned> __debug__;

/*
 * Statistics
 */

// Number of reads estimated to be sequins
static unsigned __nSeq__ = 0;

// Number of reads estimated to be genome (not sequins)
static unsigned __nGen__ = 0;

template <typename Out> void split(const std::string &s, char delim, Out result)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

/*
 * Initlaize k-mers spanning variants, they are useful for estimating allele frequency
 */

static void KMCancerSpan()
{
    std::ifstream r("CancerKMSpan.txt");
    
    if (!r.good())
    {
        throw std::runtime_error("Invalid CancerKMSpan.txt");
    }
    
    std::string line;
    while (std::getline(r, line))
    {
        std::vector<std::string> toks;
        split(line, '\t', std::back_inserter(toks));
        assert(!toks.empty());
        
        if (toks[0] == "Name") { continue; }

        __span__[toks[1]] = 0; // Normal
        __span__[toks[2]] = 0; // Reverse complement
    }
    
    r.close();
    assert(!__span__.empty());
}

/*
 * Initalize all k-mers, useful for estimating dilution
 */

static void KMSequins()
{
    std::ifstream r("SequinsKMAll.txt");
    
    if (!r.good())
    {
        throw std::runtime_error("Invalid SequinsKMAll.txt");
    }
    
    std::string line;
    while (std::getline(r, line))
    {
        std::vector<std::string> toks;
        split(line, '\t', std::back_inserter(toks));
        assert(!toks.empty());
        
        if (toks[0] == "Name") { continue; }
        
        __all__[toks[1]] = 0; // Normal
        __all__[toks[2]] = 0; // Reverse complement
    }
    
    r.close();
    assert(!__all__.empty());
}

void KMInit()
{
    ProgramOptions o;
    o.k = 31;
    o.index = "sequins.fa.index";

    __seqsIndex__ = std::shared_ptr<KmerIndex>(new KmerIndex(o));
    __seqsIndex__->load(o);
    
    KMSequins();
    KMCancerSpan();
 }

static void KMCount(const char *s)
{
    Kmer::k = 31;
    KmerIterator iter(s), end;
    
    // Number of k-mers that are sequins
    unsigned isSeq = 0;
    
    // Number of k-mers that are genome
    unsigned isGen = 0;
    
    for (int i = 0; iter != end; ++i, ++iter)
    {
        auto k = iter->first.rep().toString();
        
        /*
         * Does the k-mer span sequin variants?
         */
        
        if (__span__.count(k))
        {
            __span__[k]++;
        }

        /*
         * Is this one of the reference k-mers?
         */
        
        if (__all__.count(k))
        {
            std::vector<std::pair<KmerEntry, int> > v;
            __seqsIndex__->match(k.c_str(), 31, v);
            
            if (v.empty())
            {
                assert(false);
            }
            else
            {
               // std::cout << "GOOD" << std::endl;
            }
            
            __all__[k]++;
            isSeq++;
        }
        else
        {
            isGen++;
        }

#ifdef DEBUG
        __debug__[k]++;
#endif
    }
    
    if (isSeq > isGen)
    {
        __nSeq__++;
    }
    else
    {
        __nGen__++;
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
    
    for (const auto &i : __debug__)
    {
        w << i.first << "\t" << i.second << std::endl;
    }
    
    w.close();
#endif
}

void KMResults(unsigned &nGen, unsigned &nSeq, std::map<std::string, unsigned> &span)
{
    nGen = __nGen__;
    nSeq = __nSeq__;
    span = __span__;

    //std::map<std::string, unsigned>
//    std::ofstream w("KallistoCount.txt");
//
//    for (const auto &i : __span__)
//    {
//        w << i.first << "\t" << i.second << std::endl;
//    }
//    w.close();
}
