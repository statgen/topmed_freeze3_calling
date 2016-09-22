/* The MIT License

   Copyright (c) 2014 Adrian Tan <atks@umich.edu> and Hyun Min Kang <hmkang@umich.edu>

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

#ifndef SEX_PLOIDY_MAP_H
#define SEX_PLOIDY_MAP_H

#include "htslib/kseq.h"
#include "htslib/vcf.h"
#include "hts_utils.h"

#define PLOIDY_TYPE_AUTOSOME 0
#define PLOIDY_TYPE_X 1
#define PLOIDY_TYPE_Y 2
#define PLOIDY_TYPE_MT 3

class sex_ploidy_map {
 public:
  std::string x_label;
  std::string y_label;
  std::string mt_label;
  int32_t x_rid;
  int32_t y_rid;
  int32_t mt_rid;
  int32_t x_start;
  int32_t x_stop;

  std::map<std::string,int> mSex;  // sample ID to sex info
  std::vector<int32_t> vSex;       // vector of sex info
  int32_t n_males;

  // setters
 sex_ploidy_map(std::string xLabel, std::string yLabel, std::string mtLabel, int32_t xStart, int32_t xStop) :
  x_label(xLabel), y_label(yLabel), mt_label(mtLabel), x_start(xStart), x_stop(xStop), n_males(0) {}
 
  int32_t set_rids_from_hdr(bcf_hdr_t* hdr);
  int32_t load_sex_map_file(const char* filename, bcf_hdr_t* hdr);

  // getters
  int32_t get_ploidy_type(bcf1_t* v);
};

#endif
