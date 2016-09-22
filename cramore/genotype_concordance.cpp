#include "genotype_concordance.h"

FamilyConcordance::FamilyConcordance(int32_t numKids) {
  nKids = numKids;
  nCells = ( 0x01 << (numKids*2+4) ); // 4^(nKids+2)
  counts.resize(nCells,0);
}

// assume that genotypes are unphased, 0(missing), 1, 2, 3
void FamilyConcordance::addGenotype(int32_t gDad, int32_t gMom, int32_t gKid) {
  std::vector<int32_t> gs(3);
  gs[0] = gDad;
  gs[1] = gMom;
  gs[2] = gKid;
  addGenotype(gs);
}

// assume that genotypes are unphased, 0(missing), 1, 2, 3
void FamilyConcordance::addGenotype(int32_t gDad, int32_t gMom, std::vector<int32_t>& gKids) {
  int64_t cell = 0;
  for(int32_t i=0; i < (int32_t)gKids.size()+2; ++i) {
    int32_t g = (i == 0) ? gDad : ( (i == 1) ? gMom : gKids[i-2] );
    if ( ( g < 0 ) || ( g > 3 ) )
      error("Cannot recognize genotype %d - i=%d.. from FamilyConcordance::addGenotype(%d,%d,%d,..)",g,i,gDad,gMom,gKids[i-2]);
    
    cell = ( ( cell << 2 ) | g );
  }
  ++counts[cell];
}

// assume that genotypes are unphased, 0(missing), 1, 2, 3
void FamilyConcordance::addGenotype(std::vector<int32_t>& gFams) {
  int64_t cell = 0;
  for(int32_t i=0; i < (int)gFams.size(); ++i) {
    int32_t g = gFams[i];
    if ( ( g < 0 ) || ( g > 3 ) )
      error("Cannot recognize genotype %d - i=%d.. from FamilyConcordance::addGenotype(%d,%d,%d,...)",g,i,gFams[0],gFams[1],gFams[i]);      
    
    cell = ( ( cell << 2 ) | g );
  }
  ++counts[cell];
}

/*
int32_t FamilyConcordance::getCount(int32_t indexKid, int32_t maskDad, int32_t maskMom, int32_t maskKid) {
  int32_t sum = 0;
  for(int64_t i=0; i < nCells; ++i) {
    if ( counts[i] > 0 ) {
      if ( ( GENO2MASK(i >> (nKids*2+2)) & maskDad ) &&
	   ( GENO2MASK(i >> (nKids*2))   & maskMom ) &&
	   ( GENO2MASK(i >> (nKids-indexKid+1)*2) & maskKid ) ) {
	sum += counts[i];
      }
    }
  }
  return sum;
}
*/

int32_t FamilyConcordance::fillTrioCount(int32_t indexKid, std::vector<int32_t>& c64) {
  int32_t sum = 0;
  c64.resize(64,0);
  std::fill(c64.begin(), c64.end(), 0);
  for(int64_t i=0; i < nCells; ++i) {
    if ( counts[i] > 0 ) {
      int64_t cell = ( ( ( i >> (nKids*2-2) ) & 0x30 ) |
		       ( ( i >> (nKids*2-2) ) & 0x0c ) |
		       ( ( i >> ((nKids-indexKid-1)*2) ) & 0x03 ) );     // convert (dad,mom,kid) to 0-63 cell
      c64[cell] += counts[i];
      sum += counts[i];
    }
  }
  return sum;
}

DupConcordance::DupConcordance(int32_t numDups) {
  nDups = numDups;
  nCells = ( 0x01 << (2*numDups) );  // 4^numDups
  counts.resize(nCells,0);
}

void DupConcordance::addGenotype(std::vector<int32_t>& gDups) {
  int64_t cell = 0;
  for(int32_t i=0; i < (int)gDups.size(); ++i) {
    int32_t g = gDups[i];
    if ( ( g < 0 ) || ( g > 3 ) )
      error("Cannot recognize genotype %d - i=%d.. from DupConcordance::addGenotype(%d..%d)",g,i,gDups[0],gDups[i]);
    
    cell = ( ( cell << 2 ) | g );
  }
  ++counts[cell];
}

void DupConcordance::addGenotype(int32_t g1, int32_t g2) {
  std::vector<int32_t> gs(1,g1);
  gs.push_back(g2);
  addGenotype(gs);
}

/*
int32_t DupConcordance::getCount(int32_t index1, int32_t index2, int32_t mask1, int32_t mask2) {
  int32_t sum = 0;
  for(int64_t i=0; i < nCells; ++i) {
    if ( counts[i] > 0 ) {
      if ( ( GENO2MASK(i >> ((nDups-index1-1)*2)) & mask1 ) &&
	   ( GENO2MASK(i >> ((nDups-index2-1)*2)) & mask2 ) ) {
	sum += counts[i];
      }
    }
  }
  return sum;
}
*/

int32_t DupConcordance::fillDupCount(int32_t index1, int32_t index2, std::vector<int32_t>& c16) {
  int32_t sum = 0;
  c16.resize(16,0);
  std::fill(c16.begin(), c16.end(), 0);  
  for(int64_t i=0; i < nCells; ++i) {
    if ( counts[i] > 0 ) {
      int64_t cell = ( ( ( ( i >> ((nDups-index1-1)*2) ) & 0x03 ) << 2 ) | ( ( i >> ((nDups-index2-1)*2) ) & 0x03 ) );
      c16[cell] += counts[i];
      sum += counts[i];
    }
  }
  return sum;
}

void printTrioDupCount(htsFile* fp, std::string& hdr, std::vector<int32_t>& cnts) {
  int32_t sum = 0;
  for(int32_t i=0; i < (int32_t)cnts.size(); ++i) {
    sum += cnts[i];
  }    
  hprintf(fp, "%s\t%d", hdr.c_str(), sum);
  for(int32_t i=0; i < (int32_t)cnts.size(); ++i) {
    hprintf(fp, "\t%d", cnts[i]);
  }
  hprintf(fp,"\n");
}
