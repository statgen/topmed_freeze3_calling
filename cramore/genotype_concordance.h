/* The MIT License

   Copyright (c) 2015 Hyun Min Kang <hmkang@umich.edu>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef GENOTYPE_CONCORDANCE_H
#define GENOTYPE_CONCORDANCE_H

#include <cstdlib>
#include <stdint.h>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <list>
#include <map>
#include <queue>

#include "hts_utils.h"
#include "Error.h"

#define GENO_MASK_MIS 0x01
#define GENO_MASK_REF 0x02
#define GENO_MASK_HET 0x04
#define GENO_MASK_ALT 0x08

#define GENO2MASK(a) (0x01 << ( (a) & 0x11))

class FamilyConcordance {
 public:
  std::vector<int32_t> counts;
  int32_t nKids;
  int64_t nCells;
  
  FamilyConcordance(int32_t numKids = 0);

  // assume that genotypes are unphased, 0(missing), 1, 2, 3
  void addGenotype(int32_t gDad, int32_t gMom, int32_t gKid);  
  void addGenotype(int32_t gDad, int32_t gMom, std::vector<int32_t>& gKids);
  void addGenotype(std::vector<int32_t>& gFams);  
  //int32_t getCount(int32_t indexKid, int32_t maskDad, int32_t maskMom, int32_t maskKid);
  int32_t fillTrioCount(int32_t indexKid, std::vector<int32_t>& c64);
};

class DupConcordance {
 public:
  std::vector<int32_t> counts;
  int32_t nDups;
  int64_t nCells;

  DupConcordance(int32_t numDups = 0);
  void addGenotype(int32_t g1, int32_t g2);  
  void addGenotype(std::vector<int32_t>& gDups);
  //int32_t getCount(int32_t index1, int32_t index2, int32_t mask1, int32_t mask2);
  int32_t fillDupCount(int32_t index1, int32_t index2, std::vector<int32_t>& c16);
};

void printTrioDupCount(htsFile* fp, std::string& hdr, std::vector<int32_t>& cnts);

#endif
