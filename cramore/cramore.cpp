//#include <iostream>
//#include <string>
//#include <map>
//#include <vector>
//#include <ctime>
//#include <cstdlib>
#include <getopt.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <set>

#ifdef _OPENMP
  #include <omp.h>
#else
  #define omp_get_thread_num() 0
#endif  

#include "params.h"
#include "Error.h"
#include "PhredHelper.h"
#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/khash.h"
#include "genome_interval.h"
#include "hts_utils.h"
#include "bcf_ordered_reader.h"
#include "bcf_ordered_writer.h"
#include "estimator.h"
#include "dropseq.h"
#include "nuclear_pedigree.h"
#include "genotype_concordance.h"
#include "filter.h"
#include "utils.h"

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

int32_t runSimulContam(int32_t argc, char** argv) {
  std::string inSam;
  std::string inVcf;
  std::string inContamMap;
  std::string outPrefix;
  std::string tagGroup;
  std::string tagUMI;
  std::vector<std::string> smContams;
  int32_t capBQ = 40;
  int32_t seed = 0;
  int32_t verbose = 10000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("tag-group",&tagGroup, "Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups")
    LONG_STRING_PARAM("tag-UMI",&tagUMI, "Tag representing UMIs")

    LONG_PARAM_GROUP("Options for input VCF/BCF", NULL)
    LONG_STRING_PARAM("vcf",&inVcf, "Input VCF/BCF file")

    LONG_PARAM_GROUP("Options for specifying distribution of contamination", NULL)
    LONG_STRING_PARAM("contam-map",&inContamMap, "Each line contains [BARCODE] [ID1,FRAC_READS_1] [ID2,FRAC_READ_2] ...")
    LONG_MULTI_STRING_PARAM("sm-contam",&smContams, "List of string representing [SAMPLE_ID,FRAC_READS_1]")
    LONG_INT_PARAM("seed",&seed,"Randomization seed")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Out prefix")

    LONG_PARAM_GROUP("Analysis Options", NULL)
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  // load VCF files. This VCF should only contain hard genotypes in GT field
  std::vector<GenomeInterval> intervals;    
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  samFile* in = NULL;
  bam_hdr_t *header = NULL;
  samFile* out = NULL;

  if ( ( in = sam_open(inSam.c_str(), "r") ) == 0 ) {
    error("Cannot open file %s for reading\n",inSam.c_str());    
  }

  if ( ( out = sam_open(outPrefix.c_str(), "w") ) == 0 ) {
    error("Cannot open file %s for writing\n",outPrefix.c_str());    
  }  

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("Cannot open header from %s\n",inSam.c_str());
  }

  if ( sam_hdr_write(out, header) < 0 ) {
    error("Cannot write header to %s\n",outPrefix.c_str());    
  }

  if ( seed == 0 )
    srand(time(NULL));
  else
    srand(seed);

  bam1_t *b = bam_init1();
  //int32_t r;    

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  std::vector< std::vector<uint8_t> > v_gts;
  std::vector<double> v_afs;
  std::vector<int32_t> v_rids;
  std::vector<int32_t> v_poss;
  std::vector<char> v_refs;
  std::vector<char> v_alts;
  
  // identify samples to focus on
  std::map<std::string, int32_t> sm2idx;
  for(int32_t i=0; i < nsamples; ++i) {
    sm2idx[odr.hdr->samples[i]] = i;
  }

  // load contaminaton map
  std::map<std::string, std::vector< std::pair<int32_t, double> > > mapContam;
  
  if ( inContamMap.empty() ) {
    if ( smContams.empty() ) {
      error("Either --contan-map or --sm-contam arguments are required");
    }
    std::vector< std::pair<int32_t, double> >& v = mapContam["."]; // default barcode
    for(int32_t i=0; i < (int32_t)smContams.size(); ++i) {
      uint32_t icomma = smContams[i].find(',');
      if ( icomma == std::string::npos ) {
	error("Cannot recognize --sm-contam %s. Must contain comma",smContams[i].c_str());
      }
      std::string id = smContams[i].substr(0, icomma);
      double alpha = atof(smContams[i].substr(icomma+1).c_str());
      if ( sm2idx.find(id) == sm2idx.end() )
	error("Cannot find sample ID %s from the VCF",id.c_str());
      v.push_back( std::pair<int32_t, double>(sm2idx[id], alpha) );
    }
  }
  else {
    htsFile* hp = hts_open(inContamMap.c_str(), "r");
    if ( hp == NULL )
      error("Cannot open file %s for reading", inContamMap.c_str());

    kstring_t str = {0,0,0};
    int32_t lstr;
    int32_t nfields = 0;
    int32_t* fields = NULL;
    while( ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0 ) {
      fields = ksplit(&str, 0, &nfields);
      if ( nfields < 2 )
	error("Cannot parse line %s to extract contamination fraction",str.s);

      std::vector< std::pair<int32_t, double> >& v = mapContam[&str.s[fields[0]]];
      
      for(int32_t i=1; i < nfields; ++i) {
	char* pcomma = strchr(&str.s[fields[i]], ',');
	if ( pcomma == NULL ) {
	  error("Cannot recognize %s. Must contain comma",&str.s[fields[i]]);  
	}
	std::string id = std::string(&str.s[fields[i]], pcomma-&str.s[fields[i]]);
	double alpha = atof(pcomma+1);
	if ( sm2idx.find(id) == sm2idx.end() )
	  error("Cannot find sample ID %s from the VCF",id.c_str());
	v.push_back( std::pair<int32_t, double>(sm2idx[id], alpha) );
      }
    }
  }

  // normalize the contamination fraction
  for(std::map<std::string, std::vector< std::pair<int32_t, double> > >::iterator it = mapContam.begin(); it != mapContam.end(); ++it) {
    std::vector< std::pair<int32_t, double > >& v = it->second;
    double sumAlpha = 0;
    for(int32_t i=0; i < (int32_t)v.size(); ++i) {
      sumAlpha += v[i].second;
    }
    if ( fabs(sumAlpha - 1) > 1e-8 )
      notice("Sum of alphas for barcode %s is %lg, normalizing to be 1..", it->first.c_str(), sumAlpha);

    for(int32_t i=0; i < (int32_t)v.size(); ++i) {
      v[i].second /= sumAlpha;
    }    
  }

  notice("Started reading from the VCF file %s", inVcf.c_str());
  
  // read VCF and store genotypes
  while( odr.read(iv) ) { // read marker
    //bool skip = false;
    bcf_unpack(iv, BCF_UN_ALL);
    
    if ( iv->n_allele > 2 ) continue; // skip multi-allelics
    if ( !bcf_is_snp(iv) ) continue;  // focus only on SNPs
    
    // chrom is iv->rid
    // position is iv->pos
    // ref is iv->d.allele[0]
    // read genotypes
    int32_t rid = iv->rid;
    int32_t pos = iv->pos;
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];

    // store marker information
    v_rids.push_back(rid);
    v_poss.push_back(pos);
    v_refs.push_back(ref);
    v_alts.push_back(alt);
    
    uint32_t* p_gts = (uint32_t*)calloc(nsamples * 2, sizeof(uint32_t));
    int32_t n_gts = 0;
    
    // extract genotypes fpr selected individuals
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gts, &n_gts) < 0 ) {
      error("Cannot extract genotypes at %s:%d %c/%c", bcf_hdr_id2name(odr.hdr, rid), pos, ref, alt);
    }
    if ( n_gts != nsamples * 2 ) {
      error("Cannot extract %d genotypes at %s:%d %c/%c. Extracted only %d", nsamples * 2, bcf_hdr_id2name(odr.hdr, rid), pos, ref, alt, n_gts); 
    }

    v_gts.resize( v_gts.size() + 1 );
    //v_gts.push_back( std::vector<uint8_t>(nv, 0) );
    std::vector<uint8_t>& v = v_gts.back();
    v.resize(nsamples);

    int32_t ac = 0, an = 0;
    for(int32_t i=0; i < nsamples; ++i) {   // bi-allelic encoding of variant
      int32_t g1 = p_gts[2*i];
      int32_t g2 = p_gts[2*i+1];
      uint8_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	geno = 0;
      }
      else {
	geno = (uint8_t)(((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0) + 1);
	ac += (geno-1);
	an += 2;
      }
      v[i] = geno;
    }

    if ( ( verbose > 0 ) && ( v_rids.size() % verbose == 0 ) )
      notice("Reading %d markers from the VCF file", (int32_t)v_rids.size());

    v_afs.push_back((double)(ac+5e-7)/(double)(an+1e-6));

    free(p_gts);
  }

  notice("Finished reading %d markers from the VCF file", (int32_t)v_rids.size());

  std::map<int32_t, int32_t> rid_vcf2sam;
  vdict_t *d = (vdict_t*)odr.hdr->dict[BCF_DT_CTG];
  for(khint_t k=kh_begin(d); k < kh_end(d); ++k) {
    if ( !kh_exist(d,k) ) continue;
    int32_t vtid = kh_val(d,k).id;
    int32_t btid = bam_name2id(header, kh_key(d,k));
    if ( btid >= 0 ) 
      rid_vcf2sam[vtid] = btid;
  }

  // ensure that the chromosome orders are increasing
  int32_t prev = -1;
  for(std::map<int32_t,int32_t>::const_iterator it = rid_vcf2sam.begin(); it != rid_vcf2sam.end(); ++it) {
    if ( prev >= it->second )
      error("The references sequences are not ordered consistently between BAM and VCF file");
    prev = it->second;
  }
  
  // reading from BAM files
  int32_t ibeg = -1;
  int32_t iend = -1;
  int32_t ichr = -1;

  char gtag[2] = {0,0};
  char utag[2] = {0,0};
  if ( tagGroup.empty() ) { // do nothing
  }
  else if ( tagGroup.size() == 2 ) {
    gtag[0] = tagGroup.at(0);
    gtag[1] = tagGroup.at(1);    
  }
  else {
    error("Cannot recognize group tag %s. It is suppose to be a length 2 string",tagGroup.c_str());
  }

  if ( tagUMI.empty() ) { // do nothing
  }
  else if ( tagUMI.size() == 2 ) {
    utag[0] = tagUMI.at(0);
    utag[1] = tagUMI.at(1);    
  }
  else {
    error("Cannot recognize UMI tag %s. It is suppose to be a length 2 string",tagUMI.c_str());
  }

  int32_t nReadsKept = 0, nReadsChanged = 0;
  int32_t ret = 0;
  kstring_t readseq = {0,0,0};
  kstring_t readqual = {0,0,0};
  char base, qual;
  int32_t rpos;

  while( ( ret = sam_read1(in, header, b) ) >= 0 ) {
    // check whether any variant overlaps with the reads
    int32_t tid = b->core.tid;
    int32_t beg1 = bam_get_pos1(b);
    int32_t end1 = bam_get_end_pos1(b);

    // advance indices
    while( ( ichr < tid ) || ( ( ichr == tid ) && ( v_poss[ibeg] < beg1 ) ) ) {
      ++ibeg;
      ichr = (rid_vcf2sam.find(v_rids[ibeg]) == rid_vcf2sam.end() ? -1 : rid_vcf2sam[v_rids[ibeg]]);
    }

    if ( ichr == tid ) {
      if ( v_poss[ibeg] < beg1 ) {
	++nReadsKept;
	// do nothing, just print it
	sam_write1(out, header, b);
      }
      else {
	iend = ibeg;
	while( ( v_rids[iend] == v_rids[ibeg] ) && ( v_poss[iend] <= end1 ) ) {
	  ++iend;
	}

	// determine the originating sample first
	uint8_t *bcd = (*gtag) ? (uint8_t*) bam_aux_get(b, gtag) : NULL;
	std::vector< std::pair<int32_t,double> >& v = mapContam["."];
	if ( bcd != NULL ) {
	  if ( *bcd == 'Z' ) {
	    const char* sbcd = bam_aux2Z(bcd);
	    v = mapContam[sbcd];
	    if ( v.size() == 0 ) {
	      notice("Cannot find barcode %s. Skippingg...",sbcd);
	      continue;
	    }
	  }
	}

	if ( v.size() == 0 ) {
	  notice("Cannot find barcode . Skippingg...");
	  continue;	  
	}

	// randomly sample the originting sample
	double r = (rand()+0.5)/(RAND_MAX+1.);
	double ir = v.size()-1;
	for(int32_t i=0; i < (int32_t)v.size()-1; ++i) {
	  if ( r < v[i].second ) {
	    ir = i;
	    break;
	  }
	  r -= v[i].second;
	}


	// modify corresponding bases
	for(int32_t i=ibeg; i < iend; ++i) {
	  uint16_t iref = seq_nt16_table[(int32_t)v_refs[i]];
	  uint16_t ialt = seq_nt16_table[(int32_t)v_alts[i]];	
	  
	  // modify corresponding bases	  
	  bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)v_poss[i], base, qual, rpos, &readseq, &readqual);

	  if ( rpos >= 0 ) {
	    r = (rand()+0.5)/(RAND_MAX+1.);		
	    
	    uint8_t gt = v_gts[i][v[ir].first];
	    bool ref = false;
	    if ( gt == 0 ) {
	      ref = (r < v_afs[i]) ? false : true;
	    }
	    else if ( gt == 1 ) {
	      ref = true;
	    }
	    else if ( gt == 2 ) {
	      ref = ( r < 0.5 ) ? false : true;
	    }
	    else if ( gt == 3 ) {
	      ref = false;
	    }

	    r = (rand()+0.5)/(RAND_MAX+1.);
	    if ( r < phredConv.phred2Err[readqual.s[rpos]-33] ) {
	      ref = ( ref ? false : true );
	    }
	    bam_get_seq(b)[rpos] = (ref ? iref : ialt);

	    notice("Modified: %d %d %d %d\n", ichr, v_poss[i], rpos, bam_get_seq(b)[rpos]);
	  }
	}
	sam_write1(out, header, b);	  
      }
    }
  }

  notice("Finished writing BAM file. %d reads unchanged, %d reads changed", nReadsKept, nReadsChanged);

  odr.close();
  hts_close(in);
  hts_close(out);

  return 0;
}

int32_t runVerifyPairID(int32_t argc, char** argv) {
  std::string inSam; // SAM, BAM, or CRAM
  std::string inVcf; // VCF or VCF with GT fields
  std::string outPrefix;
  std::string tagGroup;
  std::string tagUMI;
  int32_t capBQ = 40;
  int32_t minBQ = 13;
  int32_t minMQ = 20;
  int32_t minTD = 0;
  int32_t qcExclFlag = 0x0f04;
  std::vector<std::string> smIDs;
  std::vector<double> gridAlpha;
  std::vector<double> gridASE;
  int32_t verbose = 10000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&inSam, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("tag-group",&tagGroup, "Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups")
    LONG_STRING_PARAM("tag-UMI",&tagUMI, "Tag representing UMIs")

    LONG_PARAM_GROUP("Options for input VCF/BCF", NULL)
    LONG_STRING_PARAM("vcf",&inVcf, "Input VCF/BCF file")
    LONG_MULTI_STRING_PARAM("sm",&smIDs, "List of sample IDs to compare to (default: use all)")    

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Out prefix")

    LONG_PARAM_GROUP("Analysis Options", NULL)
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
    LONG_INT_PARAM("min-BQ", &minBQ, "Minimum base quality to consider (lower BQ will be skipped)")
    LONG_INT_PARAM("min-MQ", &minMQ, "Minimum mapping quality to consider (lower MQ will be ignored)")
    LONG_INT_PARAM("min-TD", &minTD, "Minimum distance to the tail (lower will be ignored)")
    LONG_INT_PARAM("excl-flag", &qcExclFlag, "SAM/BAM FLAGs to be excluded")    
    LONG_MULTI_INT_PARAM("alpha",&gridAlpha, "Grid of alpha to search for (default is 0, 0.1, 0.2, 0.3, 0.4  0.5)")
    LONG_MULTI_INT_PARAM("ase",&gridASE, "Grid of allele-specific expression to search for (default is 0, 0.2, 0.3, 0.4, 0.5) -- Not implemented")
    LONG_DOUBLE_PARAM("verbose",&verbose, "Verbose message frequency")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( gridAlpha.empty() ) {
    gridAlpha.push_back(0);    
    gridAlpha.push_back(0.02);
    gridAlpha.push_back(0.05);
    gridAlpha.push_back(0.1);
    gridAlpha.push_back(0.15);
    gridAlpha.push_back(0.2);
    gridAlpha.push_back(0.25);        
    gridAlpha.push_back(0.3);
    gridAlpha.push_back(0.35);    
    gridAlpha.push_back(0.4);
    gridAlpha.push_back(0.45);    
    gridAlpha.push_back(0.5);    
  }

  if ( gridASE.empty() ) {
    gridASE.push_back(0.1);
    gridASE.push_back(0.2);
    gridASE.push_back(0.3);
    gridASE.push_back(0.4);
    gridASE.push_back(0.5);    
  }

  // load VCF files. This VCF should only contain hard genotypes in GT field
  std::vector<GenomeInterval> intervals;    
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  samFile* in = NULL;
  bam_hdr_t *header = NULL;

  if ( ( in = sam_open(inSam.c_str(), "r") ) == 0 ) {
    error("Cannot open file %s\n",inSam.c_str());    
  }

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("Cannot open header from %s\n",inSam.c_str());
  }

  if ( outPrefix.empty() )
    error("--out parameter is missing");

  bam1_t *b = bam_init1();
  //int32_t r;    

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  std::vector< std::vector<uint8_t> > v_gts;
  std::vector<double> v_afs;
  std::vector<int32_t> v_rids;
  std::vector<int32_t> v_poss;
  std::vector<char> v_refs;
  std::vector<char> v_alts;
  
  // identify samples to focus on
  std::vector<int> vids;
  std::vector<std::string> vSMs;
  if ( smIDs.empty() ) {
    for(int32_t i=0; i < nsamples; ++i) {
      vids.push_back(i);
      vSMs.push_back(odr.hdr->samples[i]);
    }
  }
  else {
    std::set<std::string> smMap;
    for(int32_t i=0; i < (int32_t)smIDs.size(); ++i) {
      smMap.insert(smIDs[i]);
    }
    for(int32_t i=0; i < nsamples; ++i) {
      if ( smMap.find(odr.hdr->samples[i]) != smMap.end() ) {
	vids.push_back(i);
	vSMs.push_back(odr.hdr->samples[i]);
      }
    }
  }
  int32_t nv = (int32_t)vids.size();

  notice("Started reading from the VCF file %s", inVcf.c_str());
  
  // read VCF and store genotypes
  while( odr.read(iv) ) { // read marker
    //bool skip = false;
    bcf_unpack(iv, BCF_UN_ALL);
    
    if ( iv->n_allele > 2 ) continue; // skip multi-allelics
    if ( !bcf_is_snp(iv) ) continue;  // focus only on SNPs
    
    // chrom is iv->rid
    // position is iv->pos
    // ref is iv->d.allele[0]
    // read genotypes
    int32_t rid = iv->rid;
    int32_t pos = iv->pos;
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];

    // store marker information
    v_rids.push_back(rid);
    v_poss.push_back(pos);
    v_refs.push_back(ref);
    v_alts.push_back(alt);
    
    uint32_t* p_gts = (uint32_t*)calloc(nsamples * 2, sizeof(uint32_t));
    int32_t n_gts = 0;
    
    // extract genotypes fpr selected individuals
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gts, &n_gts) < 0 ) {
      error("Cannot extract genotypes at %s:%d %c/%c", bcf_hdr_id2name(odr.hdr, rid), pos, ref, alt);
    }
    if ( n_gts != nsamples * 2 ) {
      error("Cannot extract %d genotypes at %s:%d %c/%c. Extracted only %d", nsamples * 2, bcf_hdr_id2name(odr.hdr, rid), pos, ref, alt, n_gts); 
    }

    v_gts.resize( v_gts.size() + 1 );
    //v_gts.push_back( std::vector<uint8_t>(nv, 0) );
    std::vector<uint8_t>& v = v_gts.back();
    v.resize(nv);

    int32_t ac = 0, an = 0;
    for(int32_t i=0; i < nv; ++i) {   // bi-allelic encoding of variant
      int32_t g1 = p_gts[2*vids[i]];
      int32_t g2 = p_gts[2*vids[i]+1];
      uint8_t geno;
      if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	geno = 0;
      }
      else {
	geno = (uint8_t)(((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0) + 1);
	ac += (geno-1);
	an += 2;
      }
      v[i] = geno;
    }

    if ( ( verbose > 0 ) && ( v_rids.size() % verbose == 0 ) )
      notice("Reading %d markers from the VCF file", (int32_t)v_rids.size());

    v_afs.push_back((double)(ac+5e-7)/(double)(an+1e-6));

    free(p_gts);
  }

  notice("Finished reading %d markers from the VCF file", (int32_t)v_rids.size());  

  // calculate genotype likelihoods from BAM/CRAMs
  // we expect
  // Pr(B|AA,AA) = (e/3) ~ e
  // Pr(B|AA,AB) = a(1/2-e/6) ~ a/2
  // Pr(B|AA,BB) = a(1-e) ~ a
  // Pr(B|AB,AA) = (1-a)(1/2-e/6) ~ (1-a)/2
  // Pr(B|AB,AB) = (1/2-e/6) ~ 1/2
  // Pr(B|AB,BB) = (1-a)(1/2-e/6)+(1-e)a ~ (1-a)/2+a = (1+a)/2
  // Pr(B|BB,AA) = (1-a)(1-e) ~ (1-a)
  // Pr(B|BB,AB) = (1-a)(1-e)+a(1/2-e/t) ~ 1-a/2
  // Pr(B|BB,BB) = (1-e) ~ 1
  // llk for [n * n * m * a * bc] ??  0.1 0.3 0.5 0.7 0.9
  
  int32_t nAlpha = (int32_t)gridAlpha.size();

  kstring_t readseq = {0,0,0};
  kstring_t readqual = {0,0,0};	
  
  hts_idx_t *idx = sam_index_load(in, inSam.c_str());
  if ( idx == NULL )
    error("Cannot load index file for %s",inSam.c_str());

  std::vector< std::vector<uint8_t> > bqs;
  std::vector< std::vector<uint8_t> > alleles; // 0 - REF, 1 - ALT, 2 - OTHER
  std::vector< std::vector<uint32_t> > ibcds;
  std::map< std::string, uint32_t > bcMap;
  std::vector< int32_t > bcReads(1,0);
  std::vector< int32_t > bcVars(1,0);  
  bcMap["."] = 0;
  
  char reg[255];
  char gtag[2] = {0,0};
  char utag[2] = {0,0};  

  int32_t nReadsAll = 0, nReadsPass = 0, nReadsRedundant = 0, nReadsN = 0, nReadsLQ = 0;

  if ( tagGroup.empty() ) { // do nothing
  }
  else if ( tagGroup.size() == 2 ) {
    gtag[0] = tagGroup.at(0);
    gtag[1] = tagGroup.at(1);    
  }
  else {
    error("Cannot recognize group tag %s. It is suppose to be a length 2 string",tagGroup.c_str());
  }

  if ( tagUMI.empty() ) { // do nothing
  }
  else if ( tagUMI.size() == 2 ) {
    utag[0] = tagUMI.at(0);
    utag[1] = tagUMI.at(1);    
  }
  else {
    error("Cannot recognize UMI tag %s. It is suppose to be a length 2 string",tagUMI.c_str());
  }  
  
  for(int32_t j=0; j < (int32_t)v_poss.size(); ++j) {  // process each marker separately
    if ( (j+1) % 10000 == 0 )
      notice("Finished Processing %d reads across %d variants across %d barcodes, filtering %d (%.2lf%%) reads, including %d (%.2lf%%) gapped reads, %d (%.2lf%%) low quality reads, and %d (%.2lf%%) redundant/qcfail reads from the BAM file %s", nReadsPass, j+1, (int32_t)bcMap.size(), nReadsAll - nReadsPass, 100.0 * (nReadsAll - nReadsPass) / nReadsAll, nReadsN, 100.0 * nReadsN / nReadsAll, nReadsLQ, 100.0 * nReadsLQ / nReadsAll, nReadsRedundant, 100.0 * nReadsRedundant / nReadsAll,  inSam.c_str());            
    
    bqs.resize( bqs.size() + 1 );
    alleles.resize( alleles.size() + 1 );
    ibcds.resize( ibcds.size() + 1 );

    sprintf(reg, "%s:%d-%d", bcf_hdr_id2name(odr.hdr, v_rids[j]), v_poss[j], v_poss[j]);

    std::vector<uint8_t>& v_bq = bqs.back();
    std::vector<uint8_t>& v_al = alleles.back();
    std::vector<uint32_t>& v_ibcd = ibcds.back();
    std::set<std::string> sUMI;    
    
    char ref = v_refs[j];
    char alt = v_alts[j];
    char base, qual;
    int32_t rpos;

    std::set<int32_t> s_ibcds;

    //std::string sal;
    //sal += ref;
    //sal += alt;
    //sal += "/";
    
    hts_itr_t* itr = bam_itr_querys(idx, header, reg);
    while( sam_itr_next(in, itr, b) >= 0 ) {
      ++nReadsAll;

      //fprintf(stderr,"\n%d ",b->core.qual);

      if ( b->core.flag & qcExclFlag ) {
	++nReadsRedundant;	
	continue;
      }
      
      if ( b->core.qual < minMQ ) {
	++nReadsN;
	continue;
      }
           
 
      bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)v_poss[j], base, qual, rpos, &readseq, &readqual);
      
      if ( base == 'N' ) {
	++nReadsN;
	continue;
      }

      if ( qual-33 < minBQ ) { ++nReadsLQ; continue; }
      if ( rpos < minTD-1 ) { ++nReadsLQ; continue; }
      if ( rpos + minTD > b->core.l_qseq ) { ++nReadsLQ; continue; }

      uint8_t *bcd = (*gtag) ? (uint8_t*) bam_aux_get(b, gtag) : NULL;
      uint32_t ibcd = 0;
      char* sbcd = ".";
      if ( bcd != NULL ) {
	if ( *bcd == 'Z' ) {
	  sbcd = bam_aux2Z(bcd);
	  if ( bcMap.find(sbcd) == bcMap.end() ) {
	    ibcd = bcMap.size();
	    bcMap[sbcd] = ibcd;
	    bcReads.resize(bcMap.size(),0);
	    bcVars.resize(bcMap.size(),0);	    
	  }
	  else {
	    ibcd = bcMap[sbcd];
	  }
	}
      }

      uint8_t *umi = (*utag) ? (uint8_t*) bam_aux_get(b, utag) : NULL;
      if ( umi != NULL ) {
	if ( *umi == 'Z' ) {
	  char* sumi = bam_aux2Z(umi);
	  std::string bcumi = std::string(sbcd)+sumi;
	  if ( sUMI.find(sumi) == sUMI.end() ) {
	    sUMI.insert(sumi);
	  }
	  else {
	    ++nReadsRedundant;
	    continue; // skip UMI already seen at the same position
	  }
	}
      }

      v_al.push_back( ( base == ref ) ? 0 : ( ( base == alt ) ? 1 : 2 ) );
      v_bq.push_back( qual-33 > capBQ ? capBQ : qual-33 );
      v_ibcd.push_back( ibcd );
      ++bcReads[ibcd];
      s_ibcds.insert(ibcd);

      //sal += (( base == ref ) ? "0" : ( ( base == alt ) ? "1" : "2" ) );
    //sal += qual;
    //sal += (char)(33+rpos);
    //sal += ",";
      
      ++nReadsPass;
    }
    sam_itr_destroy(itr);

    //if ( !v_al.empty() ) 
    //notice("%s",sal.c_str());

    for(std::set<int32_t>::const_iterator it = s_ibcds.begin(); it != s_ibcds.end(); ++it) {
      ++bcVars[*it];
    }
  }

  odr.close();

  notice("Finished processing %d reads across %d variants across %d barcodes, filtering %d (%.2lf%%) reads, including %d (%.2lf%%) gapped reads, %d (%.2lf%%) low quality reads, and %d (%.2lf%%) redundant/qcfail reads from the BAM file %s", nReadsPass, (int32_t)v_poss.size(), (int32_t)bcMap.size(), nReadsAll - nReadsPass, 100.0 * (nReadsAll - nReadsPass) / nReadsAll, nReadsN, 100.0 * nReadsN / nReadsAll, nReadsLQ, 100.0 * nReadsLQ / nReadsAll, nReadsRedundant, 100.0 * nReadsRedundant / nReadsAll,  inSam.c_str());      

  // start evaluating genotype concordances
  // calculate for (nBcd) x (nInds) to find the best matching genotypes first
  notice("Starting to identify best matching individual IDs");

  htsFile* wsingle = hts_open((outPrefix+".single").c_str(),"w");
  htsFile* wpair = hts_open((outPrefix+".pair").c_str(),"w");
  htsFile* wbest = hts_open((outPrefix+".best").c_str(),"w");

  hprintf(wsingle, "BARCODE\tSM_ID\tN_READS\tN_SNPS\tLLK1\tLLK0\n");
  hprintf(wpair,   "BARCODE\tSM1_ID\tSM2_ID\tALPHA\tN_READS\tN_SNPS\tLLK12\tLLK1\tLLK0\tLLK10\tLLK00\n");
  hprintf(wbest,   "BARCODE\tSM1_ID\tSM2_ID\tALPHA\tN_READS\tN_SNPS\tLLK12\tLLK1\tLLK0\tLLK10\tLLK00\n");    

  if ( ( wbest == NULL ) || ( wpair == NULL ) || ( wsingle == NULL ) )
    error("Cannot create %s.single, %s.pair, and %s.best files",outPrefix.c_str(), outPrefix.c_str());
  
  int32_t nbcd = (int32_t)bcMap.size();
  std::vector<double> llks(nbcd * nv, 0);
  std::vector<double> llk0s(nbcd, 0);  
  int32_t nsnp = (int32_t)bqs.size();
  for(int32_t i=0; i < nsnp; ++i) {
    if ( ( verbose > 0 ) && ( (i+1) % verbose == 0 ) )
      notice("Processing %d markers...",i+1);

    if ( bqs[i].size() == 0 ) continue;
    
    std::vector<double> GLs(nbcd * 4, 0);

    std::vector<uint8_t>& vgt = v_gts[i];
    
    for(int32_t j=0; j < (int32_t)bqs[i].size(); ++j) {
      uint8_t bq = bqs[i][j];
      uint8_t al = alleles[i][j];
      uint32_t ibcd = ibcds[i][j];

      if ( al == 2 ) continue;

      GLs[ibcd * 4 + 1] += ( (al == 0) ? phredConv.phred2LogMat3[bq] : -0.1*bq );
      GLs[ibcd * 4 + 2] += phredConv.phred2HalfLogMat3[bq];
      GLs[ibcd * 4 + 3] += ( (al == 1) ? phredConv.phred2LogMat3[bq] : -0.1*bq );      
    }

    double p0 = (1.-v_afs[i])*(1.-v_afs[i]);
    double p1 = 2. * v_afs[i] * (1.-v_afs[i]);
    double p2 = v_afs[i] * v_afs[i];

    // calculate genotype likelihoods for each batcodes
    for(int32_t j=0; j < nbcd; ++j) {
      GLs[j * 4] = ( GLs[j * 4 + 1] > GLs[j * 4 + 2] ) ? ( GLs[j * 4 + 1] > GLs[j * 4 + 3] ? GLs[j * 4 + 1] : GLs[j * 4 + 3] ) : ( GLs[j * 4 + 2] > GLs[j * 4 + 3] ? GLs[j * 4 + 2] : GLs[j * 4 + 3] );
      GLs[j*4+1] -= GLs[j*4];
      GLs[j*4+2] -= GLs[j*4];
      GLs[j*4+3] -= GLs[j*4];
      GLs[j*4] = log10(p0*pow(10.,GLs[j*4+1]) + p1*pow(10.,GLs[j*4+2]) + p2*pow(10.,GLs[j*4+3]));

      for(int32_t k=0; k < nv; ++k) {
	llks[j * nv + k] += GLs[j*4 + vgt[k]];
      }
      llk0s[j] += GLs[j*4];
    }
  }

  std::vector<std::string> sbcd(bcMap.size());
  for(std::map<std::string,uint32_t>::const_iterator it = bcMap.begin(); it != bcMap.end(); ++it) {
    sbcd[it->second] = it->first;
  }

  // find the best matching individual
  std::vector<int32_t> iBest(nbcd,0);
  std::vector<double> llkBest(nbcd,0);
  notice("Identifying best-matching individual..");
  for(int32_t i=0; i < nbcd; ++i) {
    double imax = -1;
    double maxLLK = -1e300;
    double inext = -1;
    double nextLLK = -1e300;
    for(int32_t j=0; j < nv; ++j) {
      hprintf(wsingle,"%s\t%s\t%d\t%d\t%.5lf\t%.5lf\n", sbcd[i].c_str(), vSMs[j].c_str(), bcReads[i], bcVars[i], llks[i*nv+j], llk0s[i]);
      
      if ( llks[i*nv + j] > maxLLK ) {
	inext = imax;
	nextLLK = maxLLK;
	imax = j;
	maxLLK = llks[i*nv+j];
      }
      else if ( llks[i*nv + j] > nextLLK ) {
	inext = j;
	nextLLK = llks[i*nv+j];
      }
    }

    //notice("Barcode:%s\tnReads:%d\tnVariants:%d\tBest:%s\tNext:%s\tmaxLLK=%.5lf\tnextLLK=%.5lf\t%LLK0=%.5lf", sbcd[i].c_str(), bcReads[i], bcVars[i], vSMs[imax].c_str(), vSMs[inext].c_str(), llks[i*nv+imax], nextLLK, llk0s[i]);
    iBest[i] = imax;
    llkBest[i] = maxLLK;
  }

  hts_close(wsingle);

  // start finding the next-best matching individual
  std::vector<double> pllks(nAlpha * nbcd * nv, 0);
  std::vector<double> pllk0s(nAlpha * nbcd , 0);
  std::vector<double> pllk00s(nAlpha * nbcd , 0);  
  
  for(int32_t i=0; i < nsnp; ++i) {
    if ( ( verbose > 0 ) && ( (i+1) % verbose == 0 ) )
      notice("Processing %d markers...",i+1);

    if ( bqs[i].size() == 0 ) continue;
    
    std::vector<double> pGs(nbcd * nAlpha * 16, 1.);
    std::vector<uint8_t>& vgt = v_gts[i];
    std::vector<int32_t> cbcds(nbcd,0);
    for(int32_t j=0; j < (int32_t)bqs[i].size(); ++j) {
      uint8_t bq = bqs[i][j];
      uint8_t al = alleles[i][j];
      uint32_t ibcd = ibcds[i][j];

      if ( al == 2 ) continue;

      double pR = (al == 0) ? phredConv.phred2Mat3[bq] : phredConv.phred2Err[bq];
      double pA = (al == 1) ? phredConv.phred2Mat3[bq] : phredConv.phred2Err[bq];

      double maxpG = 0;      
      for(int32_t k=0; k < nAlpha; ++k) {
	for(int32_t l=0; l < 3; ++l) {  // 1-Alpha
	  for(int32_t m=0; m < 3; ++m) { // Alpha
	    double p = 0.5*l + (m-l)*0.5*gridAlpha[k]; // %A
	    double* pG = &pGs[ibcd*16*nAlpha + k*16 + (l+1)*4 + m+1];
	    *pG *= (pR * (1-p) + pA * p);
	    if ( maxpG < *pG )
	      maxpG = *pG;
	  }
	}
      }

      for(int32_t k=0; k < nAlpha; ++k) {
	// normalize
	for(int32_t l=0; l < 3; ++l) {  // 1-Alpha
	  for(int32_t m=0; m < 3; ++m) { // Alpha
	    pGs[ibcd*16*nAlpha + k*16 + (l+1)*4 + m+1] /= maxpG;
	  }
	}
      }
      ++cbcds[ibcd];
    }

    double p0 = (1.-v_afs[i])*(1.-v_afs[i]);
    double p1 = 2.*v_afs[i]*(1.-v_afs[i]);
    double p2 = v_afs[i] * v_afs[i];
    for(int32_t j=0; j < nbcd; ++j) {
      if ( cbcds[j] == 0 ) continue;
      for(int32_t k=0; k < nAlpha; ++k) {
	for(int32_t l=0; l < 3; ++l) {
	  pGs[j*16*nAlpha + k*16 + (l+1)*4] = p0*pGs[j*16*nAlpha + k*16 + (l+1)*4 +1] + p1*pGs[j*16*nAlpha + k*16 + (l+1)*4 +2] + p2*pGs[j*16*nAlpha + k*16 + (l+1)*4 +2];
	  pGs[j*16*nAlpha + k*16 + (l+1)] = p0*pGs[j*16*nAlpha + k*16 + l + 5] + p1*pGs[j*16*nAlpha + k*16 + l + 9] + p2*pGs[j*16*nAlpha + k*16 + l + 13];
	}
	pGs[j*16*nAlpha + k*16] = p0*p0*pGs[j*16*nAlpha + k*16 + 5] + p0*p1*pGs[j*16*nAlpha + k*16 + 6] + p0*p2*pGs[j*16*nAlpha + k*16 + 7] + p1*p0*pGs[j*16*nAlpha + k*16 + 9] + p1*p1*pGs[j*16*nAlpha + k*16 + 10] + p1*p2*pGs[j*16*nAlpha + k*16 + 11] + p2*p0*pGs[j*16*nAlpha + k*16 + 13] + p2*p1*pGs[j*16*nAlpha + k*16 + 14] + p2*p2*pGs[j*16*nAlpha + k*16 + 15];
	
	for(int32_t l=0; l < 16; ++l) {
	  pGs[j*16*nAlpha + k*16 + l] = log10(pGs[j*16*nAlpha + k*16 + l]);
	}

	for(int32_t l=0; l < nv; ++l) {
	  pllks[j*nAlpha*nv + k*nv + l] += pGs[j*16*nAlpha + k*16 + vgt[iBest[j]]*4 + vgt[l]];
	}
	pllk0s[j*nAlpha + k] += pGs[j*16*nAlpha + k*16 + vgt[iBest[j]]*4];
	pllk00s[j*nAlpha + k] += pGs[j*16*nAlpha + k*16];	
      }
    }
  }

  for(int32_t i=0; i < nbcd; ++i) {
    double jBest = -1;
    double kBest = -1;
    double maxLLK = -1e300;
    for(int32_t j=0; j < nAlpha; ++j) {
      for(int32_t k=0; k < nv; ++k) {
	hprintf(wpair,"%s\t%s\t%s\t%.3lf\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n", sbcd[i].c_str(), vSMs[iBest[i]].c_str(), vSMs[k].c_str(), gridAlpha[j], bcReads[i], bcVars[i], pllks[i*nAlpha*nv + j*nv +k], pllks[i*nAlpha*nv], pllk00s[i*nAlpha], pllk0s[i*nAlpha + j], pllk00s[i*nAlpha + j]);
	
	if ( pllks[i*nAlpha*nv + j*nv + k] > maxLLK ) {
	  jBest = j;
	  kBest = k;
	  maxLLK = pllks[i*nAlpha*nv + j*nv + k];
	}
      }
    }

    //notice("Barcode: %s\tBest: (%s, %s, %.1lf%%)\tmaxLLK=%.5lf", sbcd[i].c_str(), vSMs[iBest[i]].c_str(), vSMs[kBest].c_str(), 100*gridAlpha[jBest], maxLLK);
    hprintf(wbest,"%s\t%s\t%s\t%.3lf\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n", sbcd[i].c_str(), vSMs[iBest[i]].c_str(), vSMs[kBest].c_str(), gridAlpha[jBest], bcReads[i], bcVars[i], pllks[i*nAlpha*nv + jBest*nv +kBest], pllks[i*nAlpha*nv], pllk00s[i*nAlpha], pllk0s[i*nAlpha + jBest], pllk00s[i*nAlpha + jBest]);    
  }

  notice("Finished writing output files");

  hts_close(wpair);
  hts_close(wbest);
  
  return 0;
}

struct drop_read {
  int32_t tid;
  std::string tname;
  std::string umi;
  std::string barcode;
  std::string seq;
  std::string qual;
  int32_t nh;
  int32_t nm;

  drop_read() : tid(0), nh(0), nm(0) {}
};

typedef struct drop_read drop_read_t;

struct drop_count {
  int32_t umis;
  int32_t reads;
  double fumis;
  int32_t etc;
  double etc2;

  drop_count() : umis(0), reads(0), fumis(0), etc(0), etc2(0) {}
};

typedef struct drop_count drop_count_t;


drop_read_t* bam_drop_read(bam_hdr_t* h, bam1_t* b, const char* NH, const char* NM) {
  uint8_t* s = bam_aux_get(b, NH);
  if ( !s ) {
    error("Cannot find %c%c tag in record\n", NH[0], NH[1]);
    return NULL;
  }
  int32_t vNH = bam_aux2i(s);
    
  s = bam_aux_get(b, NM);
  if ( !s ) {
    error("Cannot find %c%c tag in record\n", NM[0], NM[1]);
    return NULL;
  }
  int32_t vNM = bam_aux2i(s);
  
  drop_read_t* read = new drop_read_t;
  read->nh = vNH;
  read->nm = vNM;
  
  // extract barcode
  char *prn = bam_get_qname(b);
  char *pbc = NULL;
  char *pumi = prn;
  char *ptmp = NULL;
  while( ( ptmp = strchr(pumi, ':') ) != NULL ) {
    pbc = pumi+1;
    pumi = ptmp+1;
  }

  read->barcode.assign(pbc,pumi-pbc-1);
  std::transform(read->barcode.begin(), read->barcode.end(), read->barcode.begin(), ::toupper);
  
  read->umi.assign(pumi);
  read->tid = b->core.tid;
  read->tname = h->target_name[read->tid];

  return read;
}

int32_t runClusterMultinomEM(int32_t argc, char** argv) {
  std::string inMatrix;
  std::string outPrefix;
  double alpha = 1.0;       // pseudo-count per cell
  double thresDiff = 1e-10; // threshold to stop EM iteration
  int32_t maxIter = 100;    // maximum number of EM iteration
  int32_t nClust = 0;       // Number of clusters required
  int32_t nRestarts = 1;    // Number of restarts to pick the best model
  int32_t seed = 0;         // random seed
  int32_t nCollapseGenes = 0; // collapse genes into a specific number

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Required Options", NULL)
    LONG_STRING_PARAM("in",&inMatrix, "Input matrix int the format of R-compatible text matrix (can be gzipped)")
    LONG_STRING_PARAM("out",&outPrefix, "Output file prefix")
    LONG_INT_PARAM("k",&nClust, "Number of clusters")

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_DOUBLE_PARAM("alpha",&alpha, "Pseudo-count per cell")    
    LONG_DOUBLE_PARAM("thres",&thresDiff, "Threshold of LLK difference to terminate the EM iteration")
    LONG_INT_PARAM("restarts",&nRestarts, "Number of restarts to pick the best model")
    LONG_INT_PARAM("max-iter",&maxIter, "Number of maximum E-M iterations")
    LONG_INT_PARAM("collapse-genes",&nCollapseGenes,"Number of genes to be collapsed into to reduce parameter space")
    LONG_INT_PARAM("seed",&seed, "Seed for random number generator (default uses clock)")        
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( nClust == 0 ) {
    error("--k is a required parameter");
  }

  htsFile* wf = hts_open((outPrefix+".pis").c_str(),"w");
  if ( wf == NULL )
    error("Cannot open file %s for writing", (outPrefix+".pis").c_str());

  htsFile* hp = hts_open(inMatrix.c_str(), "r");
  if ( hp == NULL )
    error("Cannot open file %s for reading",inMatrix.c_str());

  kstring_t str = {0,0,0};  
  int32_t lstr = 0;

  // read and parse header columns
  lstr = hts_getline(hp, KS_SEP_LINE, &str);
  if ( lstr < 0 )
    error("Cannot find header line from %s",inMatrix.c_str());

  
  int32_t nfields = 0;
  int32_t* fields = NULL;  
  fields = ksplit(&str, 0, &nfields);

  std::vector<std::string> hdrs;
  for(int32_t i=0; i < nfields; ++i) {
    hdrs.push_back(std::string(&str.s[fields[i]]));
  }

  std::vector<int32_t*> R;
  std::vector<std::string> genes;
  std::vector<int64_t> rowSums;
  int64_t* colSums = NULL;

  notice("%d columns found in the header",(int32_t)hdrs.size());

  // read and parse the matrix elements
  int64_t nZero = 0;
  int64_t nSum = 0;
  int32_t nEmptyRows = 0;
  while( ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0 ) {
    fields = ksplit(&str, 0, &nfields);

    if ( R.empty() ) {
      notice("%d columns found in the second line",nfields);      
      if ( nfields == (int32_t)hdrs.size() ) { // remove one from the header
	hdrs.erase(hdrs.begin());
	notice("Ignoring the first column in the header");	
      }
      colSums = (int64_t*)calloc((int32_t)hdrs.size(),sizeof(int64_t));
    }
    else {
      if ( ( nfields != (int32_t)hdrs.size() + 1 ) && ( nfields != (int32_t)hdrs.size() + 0 ) )
	error("Inconsistent number of headers. Expected %d but observed %d",(int32_t)hdrs.size()+1, nfields);
    }

    int32_t* cnts = (int32_t*)malloc(sizeof(int32_t)*(nfields-1));
    int64_t rowSum = 0;
    for(int32_t i=1; i < nfields; ++i) {
      cnts[i-1] = atoi(&str.s[fields[i]]);
      if ( cnts[i-1] == 0 ) ++nZero;
      else {
	rowSum += cnts[i-1];
	colSums[i-1] += cnts[i-1];
      }
    }
    
    if ( rowSum == 0 ) {
      free(cnts);
      ++nEmptyRows;
    }
    else {
      genes.push_back(std::string(&str.s[fields[0]]));      
      R.push_back(cnts);
      rowSums.push_back(rowSum);
      nSum += rowSum;
    }
  }
  hts_close(hp);


  int32_t nRow = (int32_t)genes.size();
  int32_t nCol = (int32_t)hdrs.size();
  int64_t nCell = (int64_t)nRow * (int64_t)nCol;

  notice("Loaded a matrix with %d rows and %d columns after ignoring %d empty rows. Sparsity is %.5lg. Average of non-empty cells is %.5lg", nRow, nCol, nEmptyRows, (double)nZero/(double)(nCell+nEmptyRows*nCol), (double)nSum/(double)(nCell+nEmptyRows*nCol-nZero));

  if ( nCollapseGenes > 0 ) {
    std::vector< std::vector<int32_t> > group2Gene( nCollapseGenes );
    std::vector< int32_t > gene2Group( nRow, 0 );

    for(int32_t i=0; i < nRow; ++i) {
      int32_t g = (rand() % nCollapseGenes);
      group2Gene[g].push_back(i);
    }

    nEmptyRows = 0;
		
    for(int32_t i=nCollapseGenes-1; i >= 0; --i) {
      if ( group2Gene[i].empty() ) {
	++nEmptyRows;
	group2Gene.erase(group2Gene.begin() + i);
      }
      else {
	for(int32_t j=0; j < (int32_t)group2Gene[i].size(); ++j) {
	  gene2Group[group2Gene[i][j]] = i;
	}
      }
    }

    std::vector<std::string> newGenes(nCollapseGenes-nEmptyRows);
    std::vector<int64_t> newRowSums(nCollapseGenes-nEmptyRows, 0);
    std::vector<int32_t*> newR(nCollapseGenes-nEmptyRows, NULL);

    for(int32_t i=0; i < nRow; ++i) {
      int32_t g = gene2Group[i];
      if ( newGenes[g].empty() ) {
	newGenes[g] = genes[i];
	newR[g] = (int32_t*) calloc(sizeof(int32_t), nCol);
      }
      else {
	newGenes[g] += ",";
	newGenes[g] += genes[i];
      }
      newRowSums[g] += rowSums[i];
      for(int32_t j=0; j < nCol; ++j) {
	newR[g][j] += R[i][j];
      }
      free(R[i]);
    }

    genes = newGenes;
    rowSums = newRowSums;
    R = newR;
    nRow = nCollapseGenes-nEmptyRows;
    nCell = (int64_t)nRow * (int64_t)nCol;    

    notice("Collapsed the matrix with %d rows and %d columns after ignoring %d additional empty rows created during the random collpaing procedure", nRow, nCol, nEmptyRows);
  }

  // calculate the global proportion matrix
  double* p0 = new double[nRow];
  for(int32_t i=0; i < nRow; ++i) {
    p0[i] = (double)rowSums[i]/(double)nSum;
  }

  // create multiple copies of parameters for simultaneous EM
  int32_t nCxR = nClust * nRestarts;
  double* pis = (double*)calloc(nCxR,sizeof(double));
  double* Ps = (double*)calloc(nCxR * nRow,sizeof(double));
  double* Zs = (double*)calloc(nCxR * nCol,sizeof(double));
  double* llks = new double[nRestarts];
  double* llk0s = new double[nRestarts];

  if ( seed == 0 )
    srand(time(NULL));
  else
    srand(seed);

  // randomize class assignments
  for(int32_t c=0; c < nCol; ++c) {
    double* z = &Zs[c * nCxR];
    for(int32_t r=0; r < nRestarts; ++r) {
      z[(int32_t)(floor((rand()+0.5)/(RAND_MAX+1.)*nClust)) * nRestarts + r] = 1.;
    }
  }

  //notice("foo");

  // run EM iteration
  for(int32_t iter=0; iter < maxIter; ++iter) {
    // M-step for pi
    for(int32_t c=0; c < nCol; ++c) {
      double* z = &Zs[c * nCxR];    
      for(int32_t k=0; k < nCxR; ++k) {
	pis[k] += z[k];  // pi_k = \sum_c Pr(z_c = k)
      }
    }

    //for(int32_t k=0; k < nCxR; ++k) {
    //  notice("k=%d\tu_pi=%lg",k,pis[k]);
    //}    

    // normalize pi_k, and take a log
    for(int32_t r=0; r < nRestarts; ++r) {
      double sum = 0;
      for(int32_t k=0; k < nClust; ++k) {
	sum += pis[k*nRestarts+r];
      }
      for(int32_t k=0; k < nClust; ++k) {
	pis[k*nRestarts+r] = log( pis[k*nRestarts+r]/sum );
      }	
    }

    //for(int32_t k=0; k < nCxR; ++k) {
    //  notice("k=%d\tpi=%lg",k,exp(pis[k]));
    //}

    //notice("bar");    
    
    // M-step for P (without normalization)
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &Ps[g * nCxR];
      
      for(int32_t k=0; k < nCxR; ++k)
	p[k] = 0;
      
      for(int32_t c=0; c < nCol; ++c) {
	double* z = &Zs[c * nCxR];
	double r = R[g][c] + p0[g]*alpha;
	//double r = R[g][c] + alpha/nRow;
	for(int32_t k=0; k < nCxR; ++k) {
	  //double t = z[k] * r;
	  p[k] += (z[k] * r); // not normalized   \Pr(x_g|z_c=k) \propt \sum_c R_gc Pr(z_c=k)
	}
      }
    }

    //notice("goo");

    // normalize P
    for(int32_t k=0; k < nCxR; ++k) {
      double sumP = 0;
      for(int32_t g=0; g < nRow; ++g) {
	sumP += Ps[g*nCxR + k];
      }
      
      for(int32_t g=0; g < nRow; ++g) {
	Ps[g*nCxR + k] /= sumP;
      }
    }
    

    // transform p into logp
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &(Ps[g * nCxR]);
      for(int32_t k=0; k < nCxR; ++k) {
	p[k] = log(p[k]); // + pis[k];  // pi*P in log-scale
      }
    }


    for(int32_t c=0; c < nCol; ++c) {
      double* z = &(Zs[c*nCxR]);      
      for(int32_t k=0; k < nCxR; ++k) {
	z[k] = pis[k];
      }      
    }
    
    
    // E-step for Z : t(R) %*% logP
    for(int32_t g=0; g < nRow; ++g) {
      double* p = &(Ps[g * nCxR]);      
      for(int32_t c=0; c < nCol; ++c) {
	int32_t r = R[g][c];
	double* z = &(Zs[c*nCxR]);
	for(int32_t k=0; k < nCxR; ++k) {
	  z[k] += (r*p[k]);  // \log Pr(z_c = k | x) \propt \sum_g [ \log \Pr(R_gc|z_c=k) ] + \pi_k
	}
      }
    }

    for(int32_t r=0; r < nRestarts; ++r) {
      if ( iter == 0 )
	llk0s[r] = -1e300;
      else
	llk0s[r] = llks[r];
      
      llks[r] = 0;      
    }

    for(int32_t c=0; c < nCol; ++c) {
      double* z = &(Zs[c*nCxR]);      
      for(int32_t r=0; r < nRestarts; ++r) {
	double maxZ = z[r];
	for(int32_t k=1; k < nClust; ++k) {
	  if ( maxZ < z[k*nRestarts + r])
	    maxZ = z[k*nRestarts + r];
	}

	double sumZ = 0;
	for(int32_t k=0; k < nClust; ++k) {
	  double zdiff = z[k*nRestarts+r] - maxZ;	  
	  //sumZ += (z[k*nRestarts+r] = ((zdiff < -100) ? 3.7e-44 : exp(zdiff)));
	  sumZ += (z[k*nRestarts+r] = exp(zdiff));
	}

	//for(int32_t k=0; k < nClust; ++k) {
	//  z[k*nRestarts+r] /= sumZ;
	//}
	
	llks[r] += (maxZ + log(sumZ));
	//notice("%d\t%lg\t%lg\t%lg\t%lg\t%lg",c,llks[r],maxZ,sumZ,z[0],z[1]);
	if ( isnan(llks[r]) )
	  abort();
      }
    }
    double maxDiff = -1e300;
    for(int32_t r=0; r < nRestarts; ++r) {
      notice("Iter:%d\tThread %d:\tLLK=%.5lf\tDiff=%.5lg", iter, r, llks[r], llks[r]-llk0s[r]); //, exp(pis[0]), exp(pis[nRestarts]));
      if ( maxDiff < llks[r]-llk0s[r] ) {
	maxDiff = llks[r]-llk0s[r];
      }
    }

    if ( maxDiff < thresDiff ) {
      notice("All LLK differences are less than %.5lg < %.5lg", maxDiff, thresDiff);
      break;
    }

    if ( iter + 1 == maxIter )
      notice("Reached maximum iteration %d",maxIter);      
  }

  int32_t iMin = 0;
  for(int32_t r=1; r < nRestarts; ++r) {
    if ( llks[iMin] < llks[r] )
      iMin = r;
  }

  // transform P to linear scale
  for(int32_t k=0; k < nClust; ++k) {
    int32_t kr = k*nRestarts+iMin;
    double maxP = Ps[kr];
    for(int32_t g=0; g < nRow; ++g) {
      if ( maxP < Ps[g*nCxR + kr] )
	maxP = Ps[g*nCxR + kr];
    }

    double sumP = 0;
    for(int32_t g=0; g < nRow; ++g) {
      sumP += (Ps[g*nCxR + kr] = exp(Ps[g*nCxR + kr] - maxP));
    }

    for(int32_t g=0; g < nRow; ++g) {
      Ps[g*nCxR + kr] /= sumP;
    }
  }

  for(int32_t k=0; k < nClust; ++k) 
    hprintf(wf, "%g\n",exp(pis[k*nRestarts+iMin]));
  hts_close(wf);

  wf = hts_open((outPrefix+".Ps").c_str(),"w");
  for(int32_t g=0; g < nRow; ++g) {
    hprintf(wf, "%s",genes[g].c_str());    
    for(int32_t k=0; k < nClust; ++k) {
      hprintf(wf, "\t%.5lg",Ps[g*nCxR + k*nRestarts + iMin]); 
    }
    hprintf(wf, "\n");
  }
  hts_close(wf);

  wf = hts_open((outPrefix+".Zs").c_str(),"w");
  for(int32_t c=0; c < nCol; ++c) {
    double sumZ = 0;
    for(int32_t k=0; k < nClust; ++k)
      sumZ += Zs[c*nCxR + k*nRestarts + iMin];

    int32_t iBest = 0;
    for(int32_t k=1; k < nClust; ++k) {
      if ( Zs[c*nCxR + iBest*nRestarts + iMin] < Zs[c*nCxR + k*nRestarts + iMin] )
	iBest = k;
    }

    hprintf(wf, "%s\t%d",hdrs[c].c_str(), iBest+1);
    for(int32_t k=0; k < nClust; ++k)
      hprintf(wf, "\t%.5lg",Zs[c*nCxR + k*nRestarts + iMin]/sumZ);
    hprintf(wf, "\n");
  }
  hts_close(wf);    
  
  // free up the memories
  for(int32_t i=0; i < nRow; ++i) {
    free(R[i]);
  }
  delete[] llks;
  delete[] p0;
  free(pis);
  free(Zs);
  free(Ps);
  free(colSums);
  
  return 0;
}

// map a 12-letter STAMPs accounting for errors and trimming
// when length is 12, it is using 24-bits in 2bit space
// here are the information collected
// 1. 24bit maps counting the number of reads (randomly resolve Ns)
// 2. Sort the maps based on the counts and generate raw count distribution

int32_t runMapSTAMPs(int32_t argc, char** argv) {
  std::string inFastQ;
  std::string outPrefix;  
  int32_t bcLen = 12;
  int32_t bcMinLen = 9;
  int32_t umiLen = 8;
  int32_t trailLen = 3;
  
  paramList pl;
  
  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input FASTQ Options", NULL)
    LONG_STRING_PARAM("fq",&inFastQ, "Input FASTQ file (plain-text or gzipped)")
    
    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("bc-len",&bcLen, "Barcode length")
    LONG_INT_PARAM("bc-min-len",&bcMinLen, "Minimum Barcode length retained during the trimming")
    LONG_INT_PARAM("umi-len",&umiLen, "UMI length")
    LONG_INT_PARAM("trail-len",&trailLen, "Length of trailing bases")
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outPrefix, "Output prefix")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  htsFile* hp = hts_open(inFastQ.c_str(), "r");
  if ( hp == NULL )
    error("Cannot open file %s for reading",inFastQ.c_str());
  
  // read FASTQ files
  notice("Scanning the FASTQ file first time to construct the barcode map");
  kstring_t lines[4] = { {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0} };
  int32_t lens[4];

  std::vector<uint32_t> nBarcodes(bcLen+1,0);
  std::vector<uint32_t*> cntBarcodes(bcLen+1,NULL);
  std::vector<uint32_t> maxCnts(bcLen+1,0);  

  for(int32_t i=bcMinLen; i <= bcLen; ++i) {
    nBarcodes[i] =  (1 << (i*2));
    cntBarcodes[i] = (uint32_t*) calloc(nBarcodes[i], sizeof(uint32_t) );
  }

  uint32_t n2i[256] = {0};
  //memset(n2i, 0, sizeof(uint32_t)*256);
  n2i['A'] = 0; n2i['a'] = 0;
  n2i['C'] = 1; n2i['c'] = 1;
  n2i['G'] = 2; n2i['g'] = 2;
  n2i['T'] = 3; n2i['t'] = 3;    

  int32_t nRead = 0;
  
  while( ( lens[0] = hts_getline(hp, KS_SEP_LINE, &lines[0]) ) > 0 ) {
    if ( nRead % 10000000 == 0 )
      notice("Reading %d FASTQ Records", nRead);
           
    for(int32_t i=1; i < 4; ++i) {
      if ( ( lens[i] = hts_getline(hp, KS_SEP_LINE, &lines[i]) ) <= 0 ) {
	error("FASTQ number of lines are %d mod 4", i);
      }
    }

    // sequence is stored in lines[1]
    uint32_t bc = 0;
    char* s = lines[1].s;
    if ( lines[1].l < (uint32_t)bcLen )
      error("lines[1].l = %d < bcLen = %d",lines[1].l,bcLen);
    
    for(int32_t j=0; j < bcLen; ++j) {
      bc = ( (bc << 2) + ( (s[j] == 'N') ? (rand() % 4) : n2i[(int32_t)s[j]] ) );
      if ( j+1 >= bcMinLen ) {
	++(cntBarcodes[j+1][bc]);
	if ( maxCnts[j+1] < cntBarcodes[j+1][bc] )
	  maxCnts[j+1] = cntBarcodes[j+1][bc];
      }
    }

    ++nRead;
  }

  notice("Loaded total of %d reads", nRead);

  // calculate the summary statistics of barcode counts
  std::vector< std::vector<uint32_t> > bcHist(bcLen+1);
  for(int32_t i=bcMinLen; i <= bcLen; ++i) {
    bcHist[i].resize(maxCnts[i]+1);
    for(uint32_t j = 0; j < nBarcodes[i]; ++j)
      ++(bcHist[i][cntBarcodes[i][j]]);

    notice("Writing raw barcode histogram for length %d - maxCnt is %d",i, maxCnts[i]);

    char buf[255];
    sprintf(buf,"%d",i);

    htsFile* wf = hts_open((outPrefix+".raw."+buf+".hist").c_str(),"w");
    uint32_t sum = 0;
    for(int32_t j=maxCnts[i]; j >= 0; --j) {
      if ( bcHist[i][j] > 0 ) {
	sum += bcHist[i][j];
	hprintf(wf, "%d\t%u\t%u\n", j, bcHist[i][j], sum);
	//fprintf(stderr, "%u\t%u\n", j, bcHist[i][j]);	
      }
    }
    hts_close(wf);
  }

  for(int32_t i=bcMinLen; i <= bcLen; ++i) {
    free(cntBarcodes[i]);
  }

  return 0;
}

/*
typedef struct {
  int min_baseQ, tid, max_bases;
  uint16_t * bases;
  samFile* fp;
  bam_hdr_t* h;
  char* ref;
  int len;
  faidx_t *fai;
  errmod_t *em;
} ct_t;

static int read_aln(void* data, bam1_t* b) {
  int ret;
  ct_t *g = (ct_t*) data;
  while(1) {
    ret = sam_read(g->fp, g->h, b);
    if ( ret < 0 ) break;
    if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
    if ( g->fai && b->core.tid >= 0 ) {
      if (b->core.tid != g->tid) { // then load the sequence
	free(g->ref);
	g->ref = fai_fetch(g->fai, g->h->target_name[b->core.tid], &g->len);
	g->tid = b->core.tid;
      }
      bam_prob_realn_core(b, g->ref, g->len, 1<<1|1);
    }
    break;
  }
  return ret;
}
*/

std::string bam_hdr_get_sample_name(bam_hdr_t* hdr) {
  if ( !hdr )
    error("Failed to read the BAM header");

  const char *p = hdr->text;
  const char *q, *r;
  int32_t n = 0;
  std::string sm;
  while( ( q = strstr(p, "@RG" ) ) != 0 ) {
    p = q + 3;
    r = q = 0;
    if ( ( q = strstr(p, "\tID:" ) ) != 0 ) q += 4;
    if ( ( r = strstr(p, "\tSM:" ) ) != 0 ) r += 4;
    if ( r && q ) {
      char *u, *v;
      for (u = (char*)q; *u && *u != '\t' && *u != '\n'; ++u);
      for (v = (char*)r; *v && *v != '\t' && *v != '\n'; ++v);
      *u = *v = '\0';
      if ( sm.empty() )
	sm = r;
      else if ( sm.compare(r) != 0 )
	error("Multiple sample IDs are included in one BAM file - %s, %s", sm.c_str(), r);
    }
    else break;
    p = q > r ? q : r;
    ++n;
  }
  if ( sm.empty() )
    error("Sample ID information cannot be found");
  return sm;
}

int32_t getPersonGenoDepth( int32_t* gts, int32_t* dps, NuclearFamilyPerson* pPerson, std::vector<int>& genos, std::vector<int>& depths) {
  genos.clear();
  depths.clear();
  if ( pPerson == NULL ) return 0;
  else {
    int32_t nSamples = 0;
    for(int32_t i=0; i < (int32_t)pPerson->samples.size(); ++i) {
      int32_t idx = pPerson->samples[i]->index;
      if ( idx >= 0 ) {
        int32_t g1 = gts[2*idx];
        int32_t g2 = gts[2*idx+1];
        int32_t geno;
        if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
          geno = 0;
        }
        else {
          geno = bcf_alleles2gt(bcf_gt_allele(g1),bcf_gt_allele(g2))+1;
        }
        //fprintf(stderr,"%s %d %d %d %d\n",pPerson->samples[i]->sampleID.c_str(), idx, g1, g2, geno);
        genos.push_back(geno);
        depths.push_back(dps == NULL ? 0 : dps[idx]);
        ++nSamples;
      }
    }
    return nSamples;
  }
}

int32_t runMendelDupConc(int32_t argc, char** argv) {
  std::string inVcf;
  std::string inPed;
  std::string region;
  std::string outf;
  int32_t minDP = 0;
  int32_t verbose = 1000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Files", NULL)
    LONG_STRING_PARAM("vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("ped",&inPed, "Input PED file")    
    LONG_STRING_PARAM("region",&region,"Genomic region to focus on")
    LONG_INT_PARAM("minDP",&minDP,"Minimum genotype depth")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (randomly 1/n)")    

    LONG_PARAM_GROUP("Output Files", NULL)
    LONG_STRING_PARAM("out",&outf, "Output file prefix")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  notice("Analysis Started");    
  
  // sanity check of input arguments
  if ( inPed.empty() || outf.empty() || inVcf.empty() ) {
    error("--vcf, --out, --ped are required parameters");
  }

  notice("Loading pedigree file %s",inPed.c_str());
  NuclearPedigree* ped = new NuclearPedigree(inPed.c_str());
    
  std::vector<GenomeInterval> intervals;
  if ( !region.empty() ) {
    parse_intervals(intervals, "", region);
  }
  BCFOrderedReader odr(inVcf, intervals);

  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);
  for(int32_t i=0; i < nsamples; ++i)
    ped->setSampleIndex(odr.hdr->samples[i], i);

  notice("In the Original Pedigree: %d families, %d samples, %d unique individuals", (int32_t)ped->famIDmap.size(), (int32_t)ped->smIDmap.size(), ped->numPeople());

  int32_t nremoved = ped->removeSamplesWithoutIndex();
  notice("Removed %d samples not in the VCF file from the pedigree",nremoved);

  notice("Overlapping Samples Only : %d families, %d individuals, %d sequenced samples", (int32_t)ped->famIDmap.size(), ped->numPeople(), ped->numSamplesWithIndex());

  bcf1_t* iv = bcf_init();
  int32_t nread = 0, nskip = 0;

  //int32_t ns = ped->numSamplesWithIndex();
  int32_t* p_gt = NULL;
  int32_t* p_dp = NULL;
  int32_t* p_od = NULL;  
  int32_t np_gt, np_dp, np_od;

  std::vector<int32_t> dadGTs;
  std::vector<int32_t> momGTs;
  std::vector<int32_t> nKids;
  std::vector< std::vector<int32_t> > kidGTs;
  std::vector<int32_t> dadDPs;
  std::vector<int32_t> momDPs;
  std::vector< std::vector<int32_t> > kidDPs;
  std::vector<int32_t> famGTs;
  std::vector<int32_t> trioGTs(3);

  // iterate over each family
  // for each family member, we collect the following metrics
  // (4 x 4) genotype concordance matrix for trios
  // 4^{# dups} matrix for dups
  std::map<std::string, NuclearFamily*>::iterator itF;  
  std::map<NuclearFamily*, FamilyConcordance> famConc;
  std::map<NuclearFamilyPerson*, DupConcordance>    dupConc;
  // print out variant level summary
  std::vector<int32_t> c64, c16;  
  
  int32_t i, j, k;

  htsFile* wf_vf = hts_open((outf+".var.fam.conc").c_str(),"w");
  htsFile* wf_vd = hts_open((outf+".var.dup.conc").c_str(),"w");
  htsFile* wf_if = hts_open((outf+".ind.fam.conc").c_str(),"w");
  htsFile* wf_id = hts_open((outf+".ind.dup.conc").c_str(),"w");

  hprintf(wf_vf,"CHROM\tPOS\tREF\tALT\tTOTAL");
  hprintf(wf_if,"DAD\tMOM\tKID\tTOTAL");
  for(i=0; i < 64; ++i) {
    hprintf(wf_vf,"\tN%d%d%d",(int32_t)(i/16), (int32_t)((i/4) % 4), (int32_t)(i % 4));
    hprintf(wf_if,"\tN%d%d%d",(int32_t)(i/16), (int32_t)((i/4) % 4), (int32_t)(i % 4));    
  }
  hprintf(wf_vd,"CHROM\tPOS\tREF\tALT\tTOTAL");
  hprintf(wf_id,"ID1\tID2\tTOTAL");
  for(i=0; i < 16; ++i) {
    hprintf(wf_vd,"\tN%d%d",(int32_t)(i/4), (int32_t)(i%4));
    hprintf(wf_id,"\tN%d%d",(int32_t)(i/4), (int32_t)(i%4));
  }
  hprintf(wf_vf,"\n");
  hprintf(wf_if,"\n");
  hprintf(wf_vd,"\n");
  hprintf(wf_id,"\n");  
  
  for(itF = ped->famIDmap.begin(); itF != ped->famIDmap.end(); ++itF, ++i) {
    NuclearFamily* pFam = itF->second;

    if ( ( pFam->pDad ) && ( pFam->pDad->samples.size() > 1 ) ) 
      dupConc.insert( std::make_pair(pFam->pDad, DupConcordance((int32_t)pFam->pDad->samples.size())) );
    if ( ( pFam->pMom ) && ( pFam->pMom->samples.size() > 1 ) ) 
      dupConc.insert( std::make_pair(pFam->pMom, DupConcordance((int32_t)pFam->pMom->samples.size())) );
    
    if ( !pFam->pKids.empty() ) {
      famConc.insert( std::make_pair(pFam, FamilyConcordance((int32_t)pFam->pKids.size())) );

      for(i=0; i < (int32_t)pFam->pKids.size(); ++i) {
	if ( pFam->pKids[i]->samples.size() > 1 )
	  dupConc.insert( std::make_pair(pFam->pKids[i], DupConcordance((int32_t)pFam->pKids[i]->samples.size())) );
      }
    }
  }
  
  for(nread=0; odr.read(iv); ++nread) {
    bool skip = false;
    bcf_unpack(iv, BCF_UN_ALL);

    if ( iv->n_allele > 2 )
      skip = true;
    else if ( ( !intervals.empty() ) && ( ( iv->pos + 1 < intervals[0].start1 ) || ( iv->pos + 1 > intervals[0].end1 ) ) )
      skip = true;
    else {
      bool is_vntr = false;
      for(i=0; i < iv->n_allele; ++i) {
	if ( strcmp(iv->d.allele[i],"<VNTR>") == 0 )
	  is_vntr = true;
      }
      if ( is_vntr ) skip = true;      
    }
    if ( skip ) {
      ++nskip;
      continue;
    }

    if ( iv->pos % verbose == 0 ) {
      notice("Reporting whenever the position is a multiple of %d - currently processing [%s %d %s %s]",verbose, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, iv->d.allele[0], iv->d.allele[1]);    
    }

    // extract GT and AD/DP field
    if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &np_gt) < 0 )
      error("Cannot parse GT field");

    if ( minDP > 0 ) {
      if ( bcf_get_format_int32(odr.hdr, iv, "DP", &p_dp, &np_dp) < 0 ) {
	if ( bcf_get_format_int32(odr.hdr, iv, "AD", &p_dp, &np_dp) < 0 ) {
	  error("Cannot parse AD or DP field");	
	}
	else if ( bcf_get_format_int32(odr.hdr, iv, "OD", &p_od, &np_od) < 0 ) {
	  error("Cannot parse AD or DP field");		  
	}
	
	// if AD field is available, use their sum as depth (assuming biallelics);
	for(i=0; i < nsamples; ++i) {
	  p_dp[i] = p_dp[2*i] + p_dp[2*i+1] + p_od[i];
	}
      }
    }

    FamilyConcordance vFam(1);
    DupConcordance vDup(2);
    
    for(itF = ped->famIDmap.begin(); itF != ped->famIDmap.end(); ++itF, ++i) {
      NuclearFamily* pFam = itF->second;

      // calculate genotype concordance
      // first get genotypes
      famGTs.clear();
      
      int32_t nDad = getPersonGenoDepth( p_gt, p_dp, pFam->pDad, dadGTs, dadDPs);
      for(i=0; i < nDad; ++i) {
	if ( dadDPs[i] < minDP )
	  dadGTs[i] = 0;
      }
      famGTs.push_back(nDad > 0 ? dadGTs[0] : 0);
      trioGTs[0] = (nDad > 0 ? dadGTs[0] : 0);

      int32_t nMom = getPersonGenoDepth( p_gt, p_dp, pFam->pMom, momGTs, momDPs);
      for(i=0; i < nMom; ++i) {
	if ( momDPs[i] < minDP )
	  momGTs[i] = 0;
      }
      famGTs.push_back(nMom > 0 ? momGTs[0] : 0);
      trioGTs[1] = (nMom > 0 ? momGTs[0] : 0);      

      if ( nDad > 1 ) {
	dupConc[pFam->pDad].addGenotype(dadGTs);
	for( i=1; i < nDad; ++i )
	  for( k=0; k < i; ++k) 
	    vDup.addGenotype(dadGTs[k],dadGTs[i]);
      }
      if ( nMom > 1 ) {
	dupConc[pFam->pMom].addGenotype(momGTs);
	for( i=1; i < nMom; ++i )
	  for( k=0; k < i; ++k) 
	    vDup.addGenotype(momGTs[k],momGTs[i]);	
      }

      if ( !pFam->pKids.empty() ) {
	nKids.resize(pFam->pKids.size());
	kidGTs.resize(pFam->pKids.size());
	kidDPs.resize(pFam->pKids.size());
	for(j=0; j < (int32_t)pFam->pKids.size(); ++j) {
	  nKids[j] = getPersonGenoDepth( p_gt, p_dp, pFam->pKids[j], kidGTs[j], kidDPs[j]);
	  for(i=0; i < nKids[j]; ++i) {
	    if ( kidDPs[j][i] < minDP )
	      kidGTs[j][i] = 0;
	    trioGTs[2] = kidGTs[j][i];
	    vFam.addGenotype(trioGTs);
	  }
	  famGTs.push_back( (nKids[j] > 0) ? kidGTs[j][0] : 0);
	}
	
	for(j=0; j < (int32_t)pFam->pKids.size(); ++j) {
	  if ( nKids[j] > 1 ) {
	    dupConc[pFam->pKids[j]].addGenotype(kidGTs[j]);
	    for( i=1; i < nKids[j]; ++i )
	      for( k=0; k < i; ++k) 
		vDup.addGenotype(kidGTs[j][k],kidGTs[j][i]);
	  }
	}
	
	// get the duplicate concordance and trio concordance
	famConc[pFam].addGenotype(famGTs);
      }
    }

    std::string hdr;
    int32_t total = vFam.fillTrioCount(0,c64);
    catprintf(hdr, "%s\t%d\t%s\t%s",bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, iv->d.allele[0], iv->d.allele[1]);    
    if ( total > 0 ) {
      printTrioDupCount(wf_vf, hdr, c64);
    }

    total = vDup.fillDupCount(0,1,c16);
    if ( total > 0 ) {    
      printTrioDupCount(wf_vd, hdr, c16);
    }
  }

  for(itF = ped->famIDmap.begin(); itF != ped->famIDmap.end(); ++itF, ++i) {
    NuclearFamily* pFam = itF->second;
    FamilyConcordance& iFam = famConc[pFam];
    int32_t total;

    std::string dadID = pFam->pDad ? pFam->pDad->samples[0]->sampleID : ".";
    std::string momID = pFam->pMom ? pFam->pMom->samples[0]->sampleID : ".";

    std::string hdr;     
    // print duplicates for dad
    if ( ( pFam->pDad ) && ( pFam->pDad->samples.size() > 1 ) ) {
      DupConcordance& iDup = dupConc[pFam->pDad];
      for(i=1; i < (int32_t)pFam->pDad->samples.size(); ++i) {
	for(k=0; k < i; ++k) {
	  total = iDup.fillDupCount(k,i,c16);
	  if ( total > 0 ) {
	    hdr.clear();
	    catprintf(hdr,"%s\t%s",pFam->pDad->samples[k]->sampleID.c_str(),pFam->pDad->samples[i]->sampleID.c_str());
	    printTrioDupCount(wf_id, hdr, c16);
	  }
	}
      }
    }

    // print duplicates for mom
    if ( ( pFam->pMom ) && ( pFam->pMom->samples.size() > 1 ) ) {
      DupConcordance& iDup = dupConc[pFam->pMom];
      for(i=1; i < (int32_t)pFam->pMom->samples.size(); ++i) {
	for(k=0; k < i; ++k) {
	  total = iDup.fillDupCount(k,i,c16);
	  if ( total > 0 ) {
	    hdr.clear();	    
	    catprintf(hdr,"%s\t%s",pFam->pMom->samples[k]->sampleID.c_str(),pFam->pMom->samples[i]->sampleID.c_str());
	    printTrioDupCount(wf_id, hdr, c16);
	  }
	}
      }
    }    

    if ( !pFam->pKids.empty() ) {
      for(j=0; j < (int32_t)pFam->pKids.size(); ++j) {
	std::string kidID = pFam->pKids[j]->samples[0]->sampleID;
	total = iFam.fillTrioCount(j,c64);
	if ( total > 0 ) {
	  hdr.clear();	  
	  catprintf(hdr,"%s\t%s\t%s",dadID.c_str(), momID.c_str(), kidID.c_str());
	  printTrioDupCount(wf_if, hdr, c64);
	}
	if ( pFam->pKids[j]->samples.size() > 1 ) {
	  DupConcordance& iDup = dupConc[pFam->pKids[j]];
	  for(i=1; i < (int32_t)pFam->pKids[j]->samples.size(); ++i) {
	    for(k=0; k < i; ++k) {
	      total = iDup.fillDupCount(k,i,c16);
	      if ( total > 0 ) {
		hdr.clear();		
		catprintf(hdr,"%s\t%s",pFam->pKids[j]->samples[k]->sampleID.c_str(),pFam->pKids[j]->samples[i]->sampleID.c_str());
		printTrioDupCount(wf_id, hdr, c16);
	      }
	    }
	  }	  
	}
      }
    }
  }

  hts_close(wf_vf);
  hts_close(wf_vd);
  hts_close(wf_if);
  hts_close(wf_id);

  bcf_destroy(iv);

  notice("Analysis Finished");

  return 0;
}

struct _bcf_vfilter_arg {
  std::vector<std::string> required_filters; // require at least one of the
  std::string include_expr;
  std::string exclude_expr;
};

struct _bcf_gfilter_arg {
  int32_t minDP;
  int32_t minGQ;
  //int32_t minAD;

  _bcf_gfilter_arg() {
    minDP = 0;
    minGQ = 0;
  }
};

typedef struct _bcf_vfilter_arg bcf_vfilter_arg;
typedef struct _bcf_gfilter_arg bcf_gfilter_arg;

#define MASK_GT_MISS   0x01
#define MASK_GT_HOMREF 0x02
#define MASK_GT_HET    0x04
#define MASK_GT_HOMALT 0x08
#define MASK_GT_NONREF (MASK_GT_HET|MASK_GT_HOMALT)
#define MASK_GT_NOMISS (MASK_GT_HOMREF|MASK_GT_HET|MASK_GT_HOMALT)
#define MASK_GT_ALL    (MASK_GT_MISS|MASK_GT_HOMREF|MASK_GT_HET|MASK_GT_HOMALT)

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2

int32_t runVcfSampleSummary(int32_t argc, char** argv) {
  std::string inVcf;
  std::string out;
  std::string reg;
  std::vector<std::string> sumFields;
  int32_t minDistBp = 0;
  int32_t verbose = 1000;
  bool countVariants = false;  // count variants by variant type, allele type, and genotypes

  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;
  
  std::vector<int32_t> acThres;
  std::vector<double> afThres;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Analysis Options", NULL)
    LONG_MULTI_STRING_PARAM("sum-field",&sumFields, "Field values to calculate the sums")
    LONG_PARAM("count-variants",&countVariants, "Flag to turn on counting variants")    

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_INT_PARAM("min-dist-bp",&minDistBp, "Minimum distance from the previous variant in base-position")    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")    

    LONG_PARAM_GROUP("Genotype Filtering Options", NULL)    
    LONG_MULTI_INT_PARAM("ac",&acThres,"Allele count threshold to count rare/common variants")
    LONG_MULTI_DOUBLE_PARAM("af",&afThres,"Allele frequency threshold to count rare/common variants")
    LONG_INT_PARAM("minDP",&gfilt.minDP,"Minimum depth threshold for counting genotypes")
    LONG_INT_PARAM("minGQ",&gfilt.minGQ,"Minimum depth threshold for counting genotypes") 

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() ) {
    error("--in-vcf, --out are required parameters");
  }



  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;      
    }
    else {
      error("Cannot use both --include-expr and --exclude-expr options");
    }    
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr.hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr.hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }  

  notice("Started Reading site information from VCF file");

  std::map< std::string, std::vector<int64_t> > mapFieldSums;
  std::map< std::string, std::vector<int64_t> > mapFieldVars;  
  int32_t nVariant = 0;
  int32_t nsamples = bcf_hdr_nsamples(odr.hdr);

  for(int32_t i=0; i < (int32_t)sumFields.size(); ++i)
    mapFieldSums[sumFields[i]].resize(nsamples, (int64_t)0);

  std::vector<std::string> varFields;

  varFields.push_back("ALL.SNP");
  varFields.push_back("ALL.OTH");
  varFields.push_back("NREF.SNP");
  varFields.push_back("NREF.OTH");  
  varFields.push_back("REF.SNP");
  varFields.push_back("REF.OTH");  
  varFields.push_back("HET.SNP");
  varFields.push_back("HET.OTH");
  varFields.push_back("ALT.SNP");
  varFields.push_back("ALT.OTH");
  varFields.push_back("MISS.SNP");
  varFields.push_back("MISS.OTH");    

  std::sort(acThres.begin(), acThres.end());
  for(int32_t i=0; i < (int32_t)acThres.size(); ++i) {
    char buf[255];
    sprintf(buf, "AC_%d_%d.SNP", i == 0 ? 1 : acThres[i-1]+1, acThres[i]);
    varFields.push_back(buf);
    sprintf(buf, "AC_%d_%d.OTH", i == 0 ? 1 : acThres[i-1]+1, acThres[i]);
    varFields.push_back(buf);    
  }

  std::sort(afThres.begin(), afThres.end());
  for(int32_t i=0; i < (int32_t)afThres.size(); ++i) {
    char buf[255];
    sprintf(buf, "AF_%f_%f.SNP", i == 0 ? 0 : afThres[i-1], afThres[i]);
    varFields.push_back(buf);
    sprintf(buf, "AF_%f_%f.OTH", i == 0 ? 0 : afThres[i-1], afThres[i]);
    varFields.push_back(buf);    
  }
  
  for(int32_t i=0; i < (int32_t)varFields.size(); ++i) {
    mapFieldVars[varFields[i]].resize(nsamples, (int64_t)0);    
  }

  std::vector<int32_t> varMasks(varFields.size());


  int32_t* p_gt = NULL;
  int32_t n_gt = 0;
  
  int32_t* p_fld = NULL;
  int32_t n_fld = 0;
  int32_t prev_rid = -1, prev_pos = -1;
  int32_t nskip = 0;

  int32_t an = 0, ac_alloc = 0, non_ref_ac = 0;
  int32_t* ac = NULL;
    
  for(int32_t k=0; odr.read(iv); ++k) {  // read marker
    if ( k % verbose == 0 )
      notice("Processing %d markers at %s:%d. Skipped %d markers within %d-bp from the previous marker", k, bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1, nskip, minDistBp);

    // if minimum distance is specified, skip the variant
    
    if ( ( prev_rid == iv->rid ) && ( iv->pos - prev_pos < minDistBp ) ) {
      ++nskip;
      continue;
    }
    
    bcf_unpack(iv, BCF_UN_FLT);

    // check --apply-filters
    bool has_filter = req_flt_ids.empty() ? true : false;
    if ( ! has_filter ) {
      //notice("%d %d", iv->d.n_flt, (int32_t)req_flt_ids.size());
      for(int32_t i=0; i < iv->d.n_flt; ++i) {
	for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
	  if ( req_flt_ids[j] == iv->d.flt[i] )
	    has_filter = true;
	}
      }
    }

    //if ( k % 1000 == 999 ) abort();

    if ( ! has_filter ) { ++nskip; continue; }

    // check filter logic
    if ( filt != NULL ) {
      int32_t ret = filter_test(filt, iv, NULL);
      if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
      else if ( ret ) { has_filter = false; }
    }

    if ( ! has_filter ) { ++nskip; continue; }    

    ++nVariant;

    if ( countVariants ) {
      // calculate AC
      if ( acThres.size() + afThres.size() > 0 ) {
	hts_expand(int, iv->n_allele, ac_alloc, ac);
	an = 0;
	non_ref_ac = 0;
	bcf_calc_ac(odr.hdr, iv, ac, BCF_UN_INFO|BCF_UN_FMT); // get original AC and AN values from INFO field if available, otherwise calculate
	for (int32_t i=1; i<iv->n_allele; i++)
	  non_ref_ac += ac[i];
	for (int32_t i=0; i<iv->n_allele; i++)
	  an += ac[i];      
      }
      
      // determine variant type : flags are   MISSING HOMALT HET HOMREF 
      std::fill(varMasks.begin(), varMasks.end(), 0);
      if ( bcf_is_snp(iv) ) {
	varMasks[0] = MASK_GT_ALL;
	varMasks[2] = MASK_GT_NONREF;
	varMasks[4] = MASK_GT_HOMREF;
	varMasks[6] = MASK_GT_HET;
	varMasks[8] = MASK_GT_HOMALT;
	varMasks[10] = MASK_GT_MISS;
      } // for non-ref genotypes
      else {
	varMasks[1] = MASK_GT_ALL;
	varMasks[3] = MASK_GT_NONREF;
	varMasks[5] = MASK_GT_HOMREF;
	varMasks[7] = MASK_GT_HET;
	varMasks[9] = MASK_GT_HOMALT;
	varMasks[11] = MASK_GT_MISS;
      } // for non-ref genotypes
      for(int32_t i=0, j=12; i < (int32_t)acThres.size(); ++i, j += 2) {
	if ( ( non_ref_ac > (i == 0 ? 0 : acThres[i-1]) ) && ( non_ref_ac <= acThres[i] ) ) {
	  if ( bcf_is_snp(iv) ) {
	    varMasks[j] = MASK_GT_NONREF;
	  }
	  else {
	    varMasks[j+1] = MASK_GT_NONREF;	    
	  }
	}
      }
      for(int32_t i=0, j=12 + 2*acThres.size(); i < (int32_t)afThres.size(); ++i, j += 2) {
	if ( ( non_ref_ac > (i == 0 ? 0 : afThres[i-1]*an) ) && ( non_ref_ac <= afThres[i]*an ) ) {
	  if ( bcf_is_snp(iv) ) {
	    varMasks[j] = MASK_GT_NONREF;
	  }
	  else {
	    varMasks[j+1] = MASK_GT_NONREF;	    
	  }	  
	}
      }
      
      // extract genotype and apply genotype level filter
      if ( bcf_get_genotypes(odr.hdr, iv, &p_gt, &n_gt) < 0 ) {
	error("Cannot find the field GT from the VCF file at position %s:%d", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      for(int32_t i=0; i < nsamples; ++i) {
	int32_t g1 = p_gt[2*i];
	int32_t g2 = p_gt[2*i+1];
	int32_t geno;
	if ( bcf_gt_is_missing(g1) || bcf_gt_is_missing(g2) ) {
	  geno = 0;
	}
	else {
	  geno = ((bcf_gt_allele(g1) > 0) ? 1 : 0) + ((bcf_gt_allele(g2) > 0) ? 1 : 0) + 1;
	  //if ( i == 0 )
	  //notice("g1 = %d, g2 = %d, geno = %d", g1,g2,geno);
	}
	p_gt[i] = geno;
      }

      if ( gfilt.minDP > 0 ) {
	if ( bcf_get_format_int32(odr.hdr, iv, "DP", &p_fld, &n_fld) < 0 ) {
	  error("Cannot find the field DP from the VCF file at position %s:%d", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);

	}
	for(int32_t i=0; i < nsamples; ++i) {
	  if ( p_fld[i] < gfilt.minDP ) p_gt[i] = 0;
	}
      }

      if ( gfilt.minGQ > 0 ) {
	if ( bcf_get_format_int32(odr.hdr, iv, "GQ", &p_fld, &n_fld) < 0 ) {
	  error("Cannot find the field GQ from the VCF file at position %s:%d", bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);

	}
	for(int32_t i=0; i < nsamples; ++i) {
	  if ( p_fld[i] < gfilt.minGQ ) p_gt[i] = 0;
	}
      }

      // update the maps
      for(int32_t i=0; i < (int32_t)varFields.size(); ++i) {
	std::vector<int64_t>& v = mapFieldVars[varFields[i]];
	for(int32_t j=0; j < nsamples; ++j) {
	  if ( varMasks[i] & ( 0x01 << p_gt[j] ) ) {
	    //if ( j == 0 ) notice("VarMask[i] = %d, p_gt[j] = %d", varMasks[i], p_gt[j]);
	    ++v[j];
	  }
	}
      }
      //if ( rand() % 100 == 0 ) abort();
             
    }
    
    // perform sumField tasks
    for(int32_t i=0; i < (int32_t)sumFields.size(); ++i) {
      if ( bcf_get_format_int32(odr.hdr, iv, sumFields[i].c_str(), &p_fld, &n_fld) < 0 ) {
	error("Cannot find the field %s from the VCF file at position %s:%d", sumFields[i].c_str(), bcf_hdr_id2name(odr.hdr, iv->rid), iv->pos+1);
      }
      if ( nsamples != n_fld )
	error("Field %s has multiple elements",sumFields[i].c_str());
      std::vector<int64_t>& v = mapFieldSums[sumFields[i]];
      if ( (int32_t)v.size() != nsamples )
	error("mapFieldSums object does not have %s as key",sumFields[i].c_str());

      for(int32_t j=0; j < nsamples; ++j) {
	//if ( p_fld[j] != bcf_int32_missing ) {
	v[j] += p_fld[j];
	//}
      }
    }

    prev_rid = iv->rid;
    prev_pos = iv->pos;
  }

  htsFile* wf = hts_open(out.c_str(), "w");
  
  hprintf(wf, "ID\tN.VAR");
  for(int32_t i=0; i < (int32_t)varFields.size(); ++i) {
    hprintf(wf, "\t%s",varFields[i].c_str());    
  }
  for(int32_t i=0; i < (int32_t)sumFields.size(); ++i) {
    hprintf(wf, "\tSUM.%s",sumFields[i].c_str());
  }
  hprintf(wf,"\n");
  
  for(int32_t i=0; i < nsamples; ++i) {
    hprintf(wf, "%s", odr.hdr->id[BCF_DT_SAMPLE][i]);
    hprintf(wf, "\t%d", nVariant);
    for(int32_t j=0; j < (int32_t)varFields.size(); ++j) {
      hprintf(wf, "\t%lld", mapFieldVars[varFields[j]][i]);
    }    
    for(int32_t j=0; j < (int32_t)sumFields.size(); ++j) {
      hprintf(wf, "\t%lld", mapFieldSums[sumFields[j]][i]);
    }
    hprintf(wf, "\n");
  }
  hts_close(wf);
  odr.close();

  return 0;
}

int32_t runVcfSqueeze(int32_t argc, char** argv) {
  std::vector<std::string> inVcfs;
  std::string out;
  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;
  int32_t verbose = 1000;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input VCFs", NULL)
    LONG_MULTI_STRING_PARAM("in",&inVcfs, "Input VCF/BCF files")

    LONG_PARAM_GROUP("Variant Filtering Options", NULL)    
    LONG_MULTI_STRING_PARAM("apply-filter",&vfilt.required_filters, "Require at least one of the listed FILTER strings")
    LONG_STRING_PARAM("include-expr",&vfilt.include_expr, "Include sites for which expression is true")
    LONG_STRING_PARAM("exclude-expr",&vfilt.exclude_expr, "Exclude sites for which expression is true")    

    LONG_PARAM_GROUP("Genotype Filtering Options", NULL)    
    LONG_INT_PARAM("minDP",&gfilt.minDP,"Minimum depth threshold for counting genotypes")
    LONG_INT_PARAM("minGQ",&gfilt.minGQ,"Minimum depth threshold for counting genotypes") 

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcfs.empty() || out.empty() ) {
    error("--in, --out are required parameters");
  }

  std::vector<GenomeInterval> intervals;
  BCFOrderedReader* odr = new BCFOrderedReader(inVcfs[0], intervals);
  bcf1_t* iv = bcf_init();

  // handle filter string
  std::string filter_str;
  int32_t filter_logic = 0;
  if ( vfilt.include_expr.empty() ) {
    if ( vfilt.exclude_expr.empty() ) {
      // do nothing
    }
    else {
      filter_str = vfilt.exclude_expr;
      filter_logic |= FLT_EXCLUDE;
    }
  }
  else {
    if ( vfilt.exclude_expr.empty() ) {
      filter_str = vfilt.include_expr;
      filter_logic |= FLT_INCLUDE;      
    }
    else {
      error("Cannot use both --include-expr and --exclude-expr options");
    }    
  }

  filter_t* filt = NULL;
  if ( filter_logic != 0 )
    filter_init(odr->hdr, filter_str.c_str());

  // handle --apply-filtrs
  std::vector<int32_t> req_flt_ids;
  if ( !vfilt.required_filters.empty() ) {
    for(int32_t i=0; i < (int32_t)vfilt.required_filters.size(); ++i) {
      req_flt_ids.push_back(bcf_hdr_id2int(odr->hdr, BCF_DT_ID, vfilt.required_filters[i].c_str()));
    }
  }

  BCFOrderedWriter odw(out.c_str(),0);
  odw.set_hdr(odr->hdr);
  odw.write_hdr();  

  int32_t nsamples = bcf_hdr_nsamples(odr->hdr);
  int32_t n_gts = 0, n_flds = 0;
  int32_t* gts = (int32_t*) calloc( nsamples * 2, sizeof(int32_t) );
  int32_t* flds = (int32_t*) calloc( nsamples, sizeof(int32_t) );  

  // read a specific marker position
  int32_t k = 0, nskip = 0;
 
  for(int32_t l=0; l < (int32_t)inVcfs.size(); ++l) {
    notice("PROCESSING %d-th VCF file %s", l+1, inVcfs[l].c_str());
    
    for(; odr->read(iv); ++k) {  // read marker
      if ( k % verbose == 0 )
	notice("Processing %d markers at %s:%d. Skipped %d markers", k, bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1, nskip);
      
      bcf_unpack(iv, BCF_UN_ALL);
      
      // check --apply-filters
      bool has_filter = req_flt_ids.empty() ? true : false;
      if ( ! has_filter ) {
	//notice("%d %d", iv->d.n_flt, (int32_t)req_flt_ids.size());
	for(int32_t i=0; i < iv->d.n_flt; ++i) {
	  for(int32_t j=0; j < (int32_t)req_flt_ids.size(); ++j) {
	    if ( req_flt_ids[j] == iv->d.flt[i] )
	      has_filter = true;
	  }
	}
      }
      
      if ( ! has_filter ) { ++nskip; continue; }
      
      // check filter logic
      if ( filt != NULL ) {
	int32_t ret = filter_test(filt, iv, NULL);
	if ( filter_logic == FLT_INCLUDE ) { if ( !ret)  has_filter = false; }
	else if ( ret ) { has_filter = false; }
      }
      
      if ( ! has_filter ) { ++nskip; continue; }

      // extract genotype and apply genotype level filter
      if ( bcf_get_genotypes(odr->hdr, iv, &gts, &n_gts) < 0 ) {
	error("Cannot find the field GT from the VCF file at position %s:%d", bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1);
      }

      if ( gfilt.minDP > 0 ) {
	if ( bcf_get_format_int32(odr->hdr, iv, "DP", &flds, &n_flds) < 0 ) {
	  error("Cannot find the field DP from the VCF file at position %s:%d", bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1);
	  
	}
	for(int32_t i=0; i < nsamples; ++i) {
	  if ( flds[i] < gfilt.minDP ) {
	    gts[2*i] = bcf_gt_missing;
	    gts[2*i+1] = bcf_gt_missing;	    
	  }
	}
      }
      
      if ( gfilt.minGQ > 0 ) {
	if ( bcf_get_format_int32(odr->hdr, iv, "GQ", &flds, &n_flds) < 0 ) {
	  error("Cannot find the field GQ from the VCF file at position %s:%d", bcf_hdr_id2name(odr->hdr, iv->rid), iv->pos+1);
	  
	}
	for(int32_t i=0; i < nsamples; ++i) {
	  if ( flds[i] < gfilt.minGQ ) {
	    gts[2*i] = bcf_gt_missing;
	    gts[2*i+1] = bcf_gt_missing;
	  }
	}
      }

      bcf1_t* nv = bcf_init();

      nv->n_sample = iv->n_sample;
      nv->rid = iv->rid;
      nv->pos = iv->pos;
      nv->rlen = iv->rlen;
      nv->qual = iv->qual;
      nv->n_allele = iv->n_allele;

      bcf_update_alleles(odw.hdr, nv, (const char**)iv->d.allele, iv->n_allele);
      bcf_update_filter(odw.hdr, nv, iv->d.flt, iv->d.n_flt);
      
      bcf_unpack(nv, BCF_UN_ALL);

      // transfer INFO fields
      for(int32_t i=0; i < iv->n_info; ++i) {
	bcf_info_t& info = iv->d.info[i];
	if ( info.type != BCF_BT_NULL ) {
	  const char* tag = bcf_hdr_int2id(odr->hdr,BCF_DT_ID,info.key);
	  int32_t htype = bcf_hdr_id2type(odr->hdr,BCF_HL_INFO,info.key);
	  int32_t ntmp_arr = 0;
	  void* tmp_arr = NULL;
	  int32_t ret = bcf_get_info_values(odr->hdr, iv, tag, &tmp_arr, &ntmp_arr, htype);
	  if ( ret > 0 ) {
	    if ( bcf_update_info(odw.hdr, nv, tag, tmp_arr, ntmp_arr, htype) < 0 ) {
	      fprintf(stderr,"Cannot write INFO field %s\n",tag);
	      abort();
	    }
	  }
	  else {
	    fprintf(stderr,"Cannot retrieve INFO field %s\n",tag);
	    abort();
	  }
	  free(tmp_arr);
	}
      }

      bcf_update_format_int32(odw.hdr, nv, "GT", gts, nsamples * 2);
      
      odw.write(nv);
      bcf_destroy(nv);
    }
    
    delete odr;
    if ( l+1 < (int32_t)inVcfs.size() ) 
      odr = new BCFOrderedReader(inVcfs[l+1], intervals);
    else
      odr = NULL;
  }
  odw.close();

  notice("Analysis finished");

  return 0;
}

int32_t runSparseGenotype(int32_t argc, char** argv) {
  std::string inVcf;
  std::vector<std::string> inCrams;
  std::string inCramList;
  std::string out;
  std::string reg;
  int32_t threads = 1;
  double gl_adj = 0.01;
  int32_t xStart = 2699520;
  int32_t xStop = 154931044;
  std::string xLabel("X");
  std::string yLabel("Y");
  std::string mtLabel("MT");
  std::string sexMap;

  #ifdef _OPENMP
  threads = 4;
  #endif
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Input Sequence Files", NULL)
    LONG_MULTI_STRING_PARAM("in-cram",&inCrams, "Input CRAM file(s)")
    LONG_STRING_PARAM("in-cram-list",&inCramList, "File containing input CRAM files")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")

    LONG_PARAM_GROUP("Sex Chromosomes",NULL)
    LONG_STRING_PARAM("xLabel", &xLabel, "Contig name for X chromosome")
    LONG_STRING_PARAM("yLabel", &xLabel, "Contig name for Y chromosome")
    LONG_STRING_PARAM("mtLabel", &xLabel, "Contig name for MT chromosome")
    LONG_INT_PARAM("xStart", &xStart, "Start base position of non-PAR region in X chromosome")
    LONG_INT_PARAM("xStop",  &xStop,  "End base position of non-PAR region in X chromosome")    

    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_INT_PARAM("threads",&threads, "Number of threads to parallelize")
    LONG_STRING_PARAM("sex-map",&sexMap, "Sex map file, containing ID and sex (1 for male and 2 for female) for each individual")
    LONG_DOUBLE_PARAM("gl-adj",&gl_adj, "Genotype likelihood adjustment factor at homozygous sites as Pr(1|0/0) or Pr(0|1/1)")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inVcf.empty() || out.empty() || ( inCrams.empty() && inCramList.empty() ) ) {
    error("--in-vcf, --out, --in-cram (or --in-cram-list) are required parameters");
  }

  if ( ( !inCrams.empty() ) && ( ! inCramList.empty() ) ) {
    error("--in-cram-list and --in-cram options cannot be used together");
  }

  if ( ! inCramList.empty() ) {
    htsFile* fp = hts_open(inCramList.c_str(),"r");
    if ( fp == NULL )
      error("Cannot open file %s for reading", inCramList.c_str());


    int32_t lstr = 0;
    int32_t* fields = NULL;
    int32_t n = 0;
    kstring_t str = {0,0,0};      

    while( ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0 ) {
      fields = ksplit(&str, 0, &n);
      if ( n > 1 )
	error("in-cram-list file %s contains whitespace - # fields = %d, (%s, %s)", inCramList.c_str(), n, str.s + fields[0], str.s + fields[1]);
      inCrams.push_back(std::string(str.s));
    }
    hts_close(fp);
  }

  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  
  int32_t nsamples = inCrams.size();
  std::vector<std::string> v_regs;
  std::vector<int32_t> v_rids;
  std::vector<int32_t> v_poss;
  std::vector<char> v_refs;
  std::vector<char> v_alts;
  std::vector<uint8_t*> v_pls;  // stored PLs, ADs, and ODs together
  std::vector<uint8_t*> v_ads;  // stored PLs, ADs, and ODs together
  std::vector<uint8_t*> v_ods;  // stored PLs, ADs, and ODs together

  std::vector<int32_t> vSex;
  std::map<std::string,int> mSex;

  if ( !sexMap.empty() ) {
    htsFile *file = hts_open(sexMap.c_str(),"r");
    if ( file == NULL ) {
      fprintf(stderr,"ERROR: Cannot open %s\n",sexMap.c_str());
      exit(1);
    }
    kstring_t *s = &file->line;
    while( hts_getline(file,'\n',s) >= 0 ) {
      std::string ss = std::string(s->s);
      size_t idx = ss.find_first_of("\t ");
      if ( idx == std::string::npos ) {
	fprintf(stderr,"ERROR: Cannot parse line %s in %s\n",ss.c_str(), sexMap.c_str());
	exit(1);
      }
      std::string id = ss.substr(0, idx);
      int32_t sex = atoi(ss.substr(idx+1).c_str());
      
      if ( mSex.find(id) != mSex.end() ) {
	fprintf(stderr,"ERROR: Duplicate ID %s in %s\n",id.c_str(), sexMap.c_str());
	exit(1);	      
      }
      
      if ( sex == 0 ) {
	fprintf(stderr,"WARNING: Unknown sex for individual %s, assuming females\n",id.c_str());
	sex = 2;
      }
      else if ( sex > 2 ) {
	fprintf(stderr,"ERROR: Invalid sex %d for individual %s\n",sex,id.c_str());
	exit(1);
      }
      mSex[id] = sex;
    }
  }  

  notice("Started Reading site information from VCF file");
  // load all site information first
  char region[65535];
  
  while( odr.read(iv) ) {  // read marker
    bcf_unpack(iv, BCF_UN_STR);
    if ( iv->n_allele > 2 ) continue; // skip multi-allelics
    if ( !bcf_is_snp(iv) ) continue;  // focus only on SNPs

    int32_t rid = iv->rid;
    int32_t pos = iv->pos;
    char ref = iv->d.allele[0][0];
    char alt = iv->d.allele[1][0];

    sprintf(region, "%s:%d-%d", bcf_hdr_id2name(odr.hdr, rid), pos+1, pos+1);

    uint8_t* p_pls = (uint8_t*)calloc(nsamples * 3, sizeof(uint8_t));
    uint8_t* p_ads = (uint8_t*)calloc(nsamples * 2, sizeof(uint8_t));
    uint8_t* p_ods = (uint8_t*)calloc(nsamples * 1, sizeof(uint8_t));

    v_regs.push_back(region);
    v_rids.push_back(rid);
    v_poss.push_back(pos);
    v_refs.push_back(ref);
    v_alts.push_back(alt);

    v_pls.push_back(p_pls);
    v_ads.push_back(p_ads);
    v_ods.push_back(p_ods);
  }
  
  notice("Finished Reading %d site information from VCF file",(int32_t)v_poss.size());

  #ifdef _OPENMP
  omp_set_num_threads(threads);
  #else
  threads = 1;
  #endif

  std::vector<std::string> v_sms(nsamples);
  
  // Read BAM files for each sample
#pragma omp parallel for schedule(dynamic, 1)
  for(int32_t i=0; i < nsamples; ++i) {  // read each CRAM files in parallel
    samFile* in = NULL;
    bam_hdr_t* hdr;
    char base, qual;
    int32_t rpos;
    
    if ( ( in = sam_open(inCrams[i].c_str(), "r") ) == 0 ) 
      error("Cannot open SAM/BAM/CRAM file %s",inCrams[i].c_str());
    
    if ( ( hdr = sam_hdr_read(in) ) == 0 )
      error("Cannot open header from %s\n",inCrams[i].c_str());

    v_sms[i] = bam_hdr_get_sample_name(hdr);

    notice("Processing reg = %s, i=%d, SM=%s, nsamples %d",reg.c_str(), i,v_sms[i].c_str(),nsamples);
    
    bam1_t *b = bam_init1();

    kstring_t readseq = {0,0,0};
    kstring_t readqual = {0,0,0};	
    
    hts_idx_t *idx = sam_index_load(in, inCrams[i].c_str());
    if ( idx == NULL )
      error("Cannot load index file for %s",inCrams[i].c_str());
    int32_t numReads = 0;

    // read sam
    for(int32_t j=0; j < (int32_t)v_poss.size(); ++j) {  // process each marker separately
      //if ( j % 1000 == 0 )
      //notice("i=%d, j=%d, pos = %d", i, j, v_poss[j]);
      double p[3] = {1,1,1};
      int32_t ads[3] = {0,0,0};
      double pm, pe, sump;
      
      hts_itr_t* itr = bam_itr_querys(idx, hdr, v_regs[j].c_str());
      while( sam_itr_next(in, itr, b) >= 0 ) {
	bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)v_poss[j], base, qual, rpos, &readseq, &readqual);
	//free(readseq.s);
	//free(readqual.s);
	if ( qual < 34 ) qual = 34;
	if ( qual > 73 ) qual = 73;

	pm = phredConv.phred2Mat[qual-33];
	pe = phredConv.phred2Err[qual-33];	

	if ( base == v_refs[j] ) {
	  p[0] *= ( pm * (1-gl_adj) + pe * gl_adj / 3 );
	  p[1] *= ( pm / 2 + pe / 6 );
	  p[2] *= ( pm * gl_adj + pe * (1-gl_adj) / 3 );
	  ++ads[0];
	}
	else if ( base == v_alts[j] ) {
	  p[0] *= ( pm * gl_adj + pe * (1-gl_adj) / 3 );	  	  
	  p[1] *= ( pm / 2 + pe / 6 );
	  p[2] *= ( pm * (1-gl_adj) + pe * gl_adj / 3 );	  
	  ++ads[1];
	}
	else {
	  ++ads[2];
	}
	sump = p[0]+p[1]+p[2]+1e-300;
	p[0] /= sump;
	p[1] /= sump;
	p[2] /= sump;
	
	++numReads;
	//notice("%d\t%d\t%d\t%d\t%c\t%c\t%d",i,omp_get_thread_num(),omp_get_num_threads(),v_poss[j],base,qual,rpos);
      }
      sam_itr_destroy(itr);

      sump = p[0];
      if ( p[1] > sump ) sump = p[1];
      if ( p[2] > sump ) sump = p[2];
      p[0] /= sump;
      p[1] /= sump;
      p[2] /= sump;      

      v_pls[j][i * 3 + 0] = phredConv.err2Phred(p[0]);
      v_pls[j][i * 3 + 1] = phredConv.err2Phred(p[1]);
      v_pls[j][i * 3 + 2] = phredConv.err2Phred(p[2]);
      v_ads[j][i * 2 + 0] = (ads[0] > 255 ? 255 : ads[0]);
      v_ads[j][i * 2 + 1] = (ads[1] > 255 ? 255 : ads[1]);
      v_ods[j][i] = (ads[2] > 255 ? 255 : ads[2]);
    }
    hts_idx_destroy(idx);
    sam_close(in);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);    
  }

  notice("Finished reading all BAM files to calculate genotype likelihoods and allele depths");

  // write BCF output files
  BCFOrderedWriter odw(out.c_str(), 0);
  bcf_hdr_transfer_contigs((const bcf_hdr_t*) odr.hdr, odw.hdr );
  //odw->set_hdr(odr.hdr);

  odr.close();  

  bcf_hdr_append(odw.hdr, "##INFO=<ID=AVGDP,Number=1,Type=Float,Description=\"Average Depth per Sample\">\n");	
  bcf_hdr_append(odw.hdr, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total Number Allele Counts\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequency from Best-guess Genotypes\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=GC,Number=G,Type=Integer,Description=\"Genotype Counts\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=GN,Number=1,Type=Integer,Description=\"Total Number of Genotypes\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=HWEAF,Number=A,Type=Float,Description=\"Genotype likelihood based Allele Frequency assuming HWE\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=HWDGF,Number=G,Type=Float,Description=\"Genotype likelihood based Genotype Frequency ignoring HWE\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=IBC,Number=1,Type=Float,Description=\"Inbreeding Coefficients calculated from genotype likelihoods\">\n");	
  bcf_hdr_append(odw.hdr, "##INFO=<ID=HWE_SLRT,Number=1,Type=Float,Description=\"Signed LRT test statistics based Hardy Weinberg ln(Likelihood Ratio)\">\n");
  bcf_hdr_append(odw.hdr, "##INFO=<ID=ABE,Number=1,Type=Float,Description=\"Expected allele Balance towards Reference Allele on Heterozygous Sites\">\n");
  //bcf_hdr_append(odw.hdr, "##INFO=<ID=NS_NREF,Number=1,Type=Integer,Description=\"Number of samples with non-reference reads\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele Depth\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=OD,Number=1,Type=Integer,Description=\"Other Allele Depth\">\n");
  bcf_hdr_append(odw.hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scale Genotype Likelihoods\">\n");

  for(int32_t i=0; i < nsamples; ++i) {
    bcf_hdr_add_sample(odw.hdr, v_sms[i].c_str());
  }

  odw.write_hdr();
  //odw->close();

  // write each variant in parallel
#pragma omp parallel for ordered schedule(static, 1) // OPENMP parallelization
  for(int32_t j=0; j < (int32_t)v_poss.size(); ++j) {  // process each marker separately
    bcf1_t* nv = bcf_init();
    nv->rid = v_rids[j];
    nv->pos = v_poss[j];
    nv->rlen = 1;
    nv->n_sample = nsamples;

    char tmp_allele_str[4] = {0,',',0,0};
    tmp_allele_str[0] = v_refs[j];
    tmp_allele_str[2] = v_alts[j];    

    //notice("Alleles are %s %s",tmp_d_alleles[0], tmp_d_alleles[1]);
    //bcf_update_alleles(odw->hdr, nv, tmp_d_alleles, 2);
    bcf_update_alleles_str(odw.hdr, nv, tmp_allele_str);
    //notice("Successfully updated alleles to %s %s",tmp_d_alleles[0], tmp_d_alleles[1]);    

    bcf_unpack(nv, BCF_UN_ALL);

    bool isX = ( ( xLabel.compare(0, xLabel.size(), v_regs[j]) == 0 ) && ( v_regs[j].at(xLabel.size()) == ':' ) ) ? true : false;
    // calculate the allele frequencies under HWE. When calculating allele frequencies, the sex information will be ignored
    float MLE_HWE_AF[2];
    float MLE_HWE_GF[3];
    int32_t ploidy = 2; // temporarily constant
    int32_t n = 0;

    // calculate the genotypes (diploid only)
    double gp, gp_sum, max_gp;
    int32_t best_gt;
    int32_t best_a1, best_a2;
    int32_t* pls_i;
    int32_t an = 0;
    int32_t acs[2];
    int32_t gcs[3];
    float afs[3];
    int32_t* gqs = (int32_t*) calloc( nsamples, sizeof(int32_t) );
    int32_t max_gq = 0;
    int32_t* gts = (int32_t*) calloc( nsamples * 2, sizeof(int32_t) );
    int32_t* pls = (int32_t*) calloc( nsamples * 3, sizeof(int32_t) );
    int32_t* ads = (int32_t*) calloc( nsamples * 2, sizeof(int32_t) );
    int32_t* ods = (int32_t*) calloc( nsamples, sizeof(int32_t) );        
    int32_t dp_sum = 0;    

    memset(acs, 0, sizeof(int32_t)*2);
    memset(gcs, 0, sizeof(int32_t)*3);

    for(int32_t i=0; i < nsamples; ++i) {
      pls[3*i  ] = (int32_t)v_pls[j][3*i];
      pls[3*i+1] = (int32_t)v_pls[j][3*i+1];
      pls[3*i+2] = (int32_t)v_pls[j][3*i+2];

      ads[2*i  ] = (int32_t)v_ads[j][2*i];
      ads[2*i+1] = (int32_t)v_ads[j][2*i+1];
      ods[i] = (int32_t)v_ods[j][i]; 

      dp_sum += (v_ads[j][2*i] + v_ads[j][2*i+1] + v_ods[j][i]);
    }
      
    Estimator * est = new Estimator();
    est->compute_gl_af_hwe(pls, nsamples, ploidy, 2, MLE_HWE_AF, MLE_HWE_GF,  n, 1e-20);

    int32_t adSumHet[2] = {0,0};
	  
    for(int32_t i=0; i < nsamples; ++i) {
      pls_i = &pls[ i * 3 ];
      
      if ( isX && (vSex[i] == 1) ) { // haploid
	max_gp = gp_sum = gp = ( est->lt->pl2prob(pls_i[0]) * MLE_HWE_AF[0] );
	best_gt = 0; best_a1 = 0; best_a2 = 0;
	for(size_t l=1; l < 2; ++l) {
	  gp = ( est->lt->pl2prob(pls_i[ l*(l+1)/2 + l]) * MLE_HWE_AF[l] );
	  gp_sum += gp;
	  if ( max_gp < gp ) {
	    max_gp = gp;
	    best_gt = l*(l+1)/2 + l; best_a1 = l; best_a2 = l;
	  }		
	}
      }
      else { // diploid
	max_gp = gp_sum = gp = ( est->lt->pl2prob(pls_i[0]) * MLE_HWE_AF[0] * MLE_HWE_AF[0] );
	best_gt = 0; best_a1 = 0; best_a2 = 0;
	for(size_t l=1; l < 2; ++l) {
	  for(size_t m=0; m <= l; ++m) {
	    gp = ( est->lt->pl2prob(pls_i[ l*(l+1)/2 + m]) * MLE_HWE_AF[l] * MLE_HWE_AF[m] * (l == m ? 1 : 2) );
	    gp_sum += gp;
	    if ( max_gp < gp ) {
	      max_gp = gp;
	      best_gt = l*(l+1)/2 + m; best_a1 = m; best_a2 = l;
	    }
	  }
	}
	if ( best_gt == 1 ) {
	  adSumHet[0] += ads[2*i];
	  adSumHet[1] += ads[2*i+1];
	}
      }
      
      double prob = 1.-max_gp/gp_sum;  // to calculate GQ
      if ( prob <= 3.162278e-26 )
	prob = 3.162278e-26;
      if ( prob > 1 )
	prob = 1;
      
      gqs[i] = (int32_t)est->lt->prob2pl(prob);
      
      if ( ( best_gt > 0 ) && ( max_gq < gqs[i] ) )
	max_gq = gqs[i];

      gts[2*i]   = ((best_a1 + 1) << 1);
      gts[2*i+1] = ((best_a2 + 1) << 1);	    
      an += 2;             // still use diploid representation of chrX for now.
      ++acs[best_a1];
      ++acs[best_a2];
      ++gcs[best_gt];
    }
    
    for(size_t i=0; i < 2; ++i) {
      afs[i] = acs[i]/(float)an;
    }
    
    bcf_update_format_int32(odw.hdr, nv, "GT", gts, nsamples * 2);
    bcf_update_format_int32(odw.hdr, nv, "GQ", gqs, nsamples );	  
    bcf_update_format_int32(odw.hdr, nv, "AD", ads, nsamples * 2);
    bcf_update_format_int32(odw.hdr, nv, "OD", ods, nsamples );	  
    bcf_update_format_int32(odw.hdr, nv, "PL", pls, nsamples * 3);
    
    float avgdp = (float)dp_sum / (float)nsamples;
    
    nv->qual = (float) max_gq;
    bcf_update_info_float(odw.hdr, nv, "AVGDP", &avgdp, 1);	  
    bcf_update_info_int32(odw.hdr, nv, "AC", &acs[1], 1);
    bcf_update_info_int32(odw.hdr, nv, "AN", &an, 1);
    bcf_update_info_float(odw.hdr, nv, "AF", &afs[1], 1);
    bcf_update_info_int32(odw.hdr, nv, "GC", gcs, 3);
    bcf_update_info_int32(odw.hdr, nv, "GN", &nsamples, 1);
    
    if (n) {
      float* MLE_HWE_AF_PTR = &MLE_HWE_AF[1];
      bcf_update_info_float(odw.hdr, nv, "HWEAF", MLE_HWE_AF_PTR, 1);
    }

    // calculate the allele frequencies under HWD	  
    float MLE_AF[2];
    float MLE_GF[3];
    n = 0;
    est->compute_gl_af(pls, nsamples, ploidy, 2, MLE_AF, MLE_GF,  n, 1e-20);
    if (n) {
      //float* MLE_AF_PTR = &MLE_AF[1];
      //bcf_update_info_float(odw->hdr, nv, "HWDAF", MLE_AF_PTR, n_alleles-1);
      bcf_update_info_float(odw.hdr, nv, "HWDGF", &MLE_GF, 3);
    }

    if ( isX && !mSex.empty() ) { // copy only female GLs to calculate IBC and HWE_SLP
      int32_t* p_XX_pls = (int32_t*) malloc(nsamples * 3 * sizeof(int32_t));
      int32_t i, k, l;
      for(i=0, k=0; i < nsamples; ++i) {
	if ( vSex[i] == 2 ) {
	  for(l=0; l < 3; ++l)  {
	    p_XX_pls[3 * k + l] = pls[3 * i + l];
	  }
	  ++k;
	}
      }
      
      float MLE_HWE_AF_XX[2];
      float MLE_HWE_GF_XX[3];
      float MLE_AF_XX[2];
      float MLE_GF_XX[3];
      
      // calculate allele frequencies using females
      est->compute_gl_af_hwe(p_XX_pls, j, ploidy, 2, MLE_HWE_AF_XX, MLE_HWE_GF_XX,  n, 1e-20);
      est->compute_gl_af(p_XX_pls, j, ploidy, 2, MLE_AF_XX, MLE_GF_XX,  n, 1e-20);
      
      for(i=0; i < 2; ++i) {
	if ( MLE_HWE_AF_XX[i] < 1e-6 ) MLE_HWE_AF_XX[i] = 1e-6;
	if ( MLE_AF_XX[i] < 1e-6 ) MLE_AF_XX[i] = 1e-6;	      
      }
      
      for(i=0; i < 3; ++i) {
	if ( MLE_HWE_GF_XX[i] < 1e-10 ) MLE_HWE_GF_XX[i] = 1e-10;
	if ( MLE_GF_XX[i] < 1e-10 ) MLE_GF_XX[i] = 1e-10;	      
      }	    
      
      float fic = 0;
      n = 0;
      est->compute_gl_fic(p_XX_pls, j, ploidy, MLE_HWE_AF_XX, 2, MLE_GF_XX, fic, n);
      if ( isnanf(fic) ) fic = 0;	  
      if (n) {
	bcf_update_info_float(odw.hdr, nv, "IBC", &fic, 1);
      }

      float abe = (adSumHet[0] + 0.5)/(adSumHet[0] + adSumHet[1] + 1.0);
      bcf_update_info_float(odw.hdr, nv, "ABE", &abe, 1);      
      
      // calculate the LRT statistics related to HWE
      float lrts;
      //float logp;
      //int32_t df;
      n = 0;
      est->compute_hwe_lrt(p_XX_pls, j, ploidy, 2, MLE_HWE_GF_XX, MLE_GF_XX, n, lrts);
      if (n) {
	if ( lrts < 0 ) lrts = 0;
	if ( fic < 0 ) lrts = 0-lrts;
	bcf_update_info_float(odw.hdr, nv, "HWE_SLRT", &lrts, 1);
      }
      
      free(p_XX_pls);
    }
    else {
      float fic = 0;
      n = 0;
      est->compute_gl_fic(pls, nsamples, ploidy, MLE_HWE_AF, 2, MLE_GF, fic, n);
      if ( isnanf(fic) ) fic = 0;
      if (n) {
	bcf_update_info_float(odw.hdr, nv, "IBC", &fic, 1);
      }
      
      // calculate the LRT statistics related to HWE
      float lrts;
      //float logp;
      //int32_t df;
      n = 0;
      est->compute_hwe_lrt(pls, nsamples, ploidy, 2, MLE_HWE_GF, MLE_GF, n, lrts);
      if (n) {
	if ( lrts < 0 ) lrts = 0;	
	if ( fic < 0 ) lrts = 0-lrts;
	bcf_update_info_float(odw.hdr, nv, "HWE_SLRT", &lrts, 1);
      }

      float abe = (adSumHet[0] + 0.5)/(adSumHet[0] + adSumHet[1] + 1.0);
      bcf_update_info_float(odw.hdr, nv, "ABE", &abe, 1);            
    }

    delete est;

#pragma omp ordered
    //notice("Writing variant j=%d",j);
    odw.write(nv);
    bcf_destroy(nv);
    free(gts);
    free(gqs);
    free(pls);
    free(ads);
    free(ods);
    delete est;
  }

  odw.close();
  
  return 0;
}

int32_t runKallistoCount(int32_t argc, char** argv) {
  std::string inFile;
  std::string outPrefix;  
  std::string tagNH("NH");
  std::string tagNM("NM");  
  int32_t maxNM = 3;
  
  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("SAM/BAM/CRAM Input Options", NULL)
    LONG_STRING_PARAM("sam",&inFile, "Input SAM/BAM/CRAM file. Sorted by readnames")
    
    LONG_PARAM_GROUP("Additional Options", NULL)
    LONG_STRING_PARAM("tag-NH",&tagNH, "Tag indicating the number of multiple mapping")
    LONG_STRING_PARAM("tag-NM",&tagNM, "Tag indicating the number of mismatches")
    LONG_INT_PARAM("max-NM",&maxNM, "Maximum number of mismatches allowed")    
    
    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &outPrefix, "Output prefix")
  END_LONG_PARAMS();
  
  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();
  
  // sanity check of input arguments
  if ( inFile.empty() || outPrefix.empty()  ) {
    error("--sam, --out are required parameters");
  }

  samFile* in = NULL;
  bam_hdr_t *header = NULL;
  
  if ( ( in = sam_open(inFile.c_str(), "r") ) == 0 ) {
    error("Cannot open file %s\n",inFile.c_str());    
  }

  if ( ( header = sam_hdr_read(in) ) == 0 ) {
    error("Cannot open header from %s\n",inFile.c_str());
  }

  bam1_t *b = bam_init1();
  int32_t r;

  // try to read all files
  char NH[2]; NH[0] = tagNH[0]; NH[1] = tagNH[1];
  char NM[2]; NM[0] = tagNM[0]; NM[1] = tagNM[1];

  std::vector<drop_read_t*> v_reads;
  std::map<int32_t, std::string> tid2name;

  DropLibrary dl;

  while( ( r = sam_read1(in, header, b) ) >= 0 ) {
    drop_read_t* rd = bam_drop_read(header, b, NH, NM);
    if ( rd == NULL ) error("Cannot extract a read from a BAM");
    v_reads.push_back(rd);
    int32_t minNM = rd->nm;
    if ( rd->nh > 1 ) {
      for(int32_t i=1; i < rd->nh; ++i) {
	if ( ( r = sam_read1(in, header, b) ) < 0 ) {
	  error("No more read observed despite NH tag");
	}
	drop_read_t* rd2 = bam_drop_read(header, b, NH, NM);
	if ( ( rd->barcode != rd2->barcode ) || ( rd->umi != rd2->umi ) ) {
	  error("Barcode or UMI does not match\n");
	}
	if ( minNM > rd2->nm )
	  minNM = rd2->nm;
	v_reads.push_back(rd2);
      }
    }

    // count only the reads with minimal matches
    if ( minNM <= maxNM ) {
      for(size_t i=0; i < v_reads.size(); ++i) {
	rd = v_reads[i];
	if ( rd->nm == minNM ) {
	  //printf("%s\t%s\t%d\t%s\t%d\t%d\n", rd->barcode.c_str(), rd->umi.c_str(), rd->tid, rd->tname.c_str(), rd->nh, rd->nm);
	  dl.addRead( rd->barcode, rd->tid, rd->umi, rd->nh );

	  if ( tid2name.find(rd->tid) == tid2name.end() ) {
	    tid2name[rd->tid] = rd->tname;
	  }
	}
	delete rd;
      }
    }
    else {
      for(size_t i=0; i < v_reads.size(); ++i) {
	delete v_reads[i];
      }
    }
    v_reads.clear();
  }

  FILE* fpCT = fopen((outPrefix+".cell.tx.cnts").c_str(), "w");
  if ( fpCT == NULL )
    error("Cannot write file");
  fprintf(fpCT,"#Cell\tTx\tUMIs\tReads\tfUMIs\n");
  
  FILE* fpC = fopen((outPrefix+".cell.cnts").c_str(), "w");
  if ( fpC == NULL )
    error("Cannot write file");    
  fprintf(fpC,"#Cell\tUMIs\tReads\tfUMIs\tnTx\tumiENST\t%%ENST\n");

  FILE* fpT = fopen((outPrefix+".tx.cnts").c_str(), "w");
  if ( fpT == NULL )
    error("Cannot write file");    
  fprintf(fpT,"#Tx\tUMIs\tReads\tfUMIs\tnCell\n");

  std::map<int32_t,drop_count_t> t2cnt;

  notice("Writing digital expression matrices");

  for(sc_map_it_t itC = dl.mapCell.begin(); itC != dl.mapCell.end(); ++itC) {
    drop_count_t cCnt;
    
    it_map_t& mT = itC->second->mapTranscript;
    for(it_map_it_t itT = mT.begin(); itT != mT.end(); ++itT) {
      sr_map_t& mR = itT->second->mapRead;
      drop_count_t ctCnt;
      for(sr_map_it_t itR = mR.begin(); itR != mR.end(); ++itR) {
	++ctCnt.umis;
	ctCnt.reads += itR->second.first;
	ctCnt.fumis += (1.0/itR->second.second);
      }

      std::string& tname = tid2name[itT->first];

      cCnt.umis += ctCnt.umis;
      cCnt.reads += ctCnt.reads;
      cCnt.fumis += ctCnt.fumis;
      ++cCnt.etc;
      if ( tname.compare(0, 4, "ENST") == 0 ) {
	cCnt.etc2 += ctCnt.fumis;
      }

      drop_count_t& tCnt = t2cnt[itT->first];

      tCnt.umis += ctCnt.umis;
      tCnt.reads += ctCnt.reads;
      tCnt.fumis += ctCnt.fumis;
      ++tCnt.etc;      
      
      fprintf(fpCT, "%s\t%s\t%d\t%d\t%.3lf\n",itC->first.c_str(), tname.c_str(), ctCnt.umis, ctCnt.reads, ctCnt.fumis);
    }
    
    fprintf(fpC, "%s\t\t%d\t%d\t%.3lf\t%d\t%.3lf\t%.3lf\n",itC->first.c_str(), cCnt.umis, cCnt.reads, cCnt.fumis, cCnt.etc, cCnt.etc2, cCnt.etc2/cCnt.fumis);   
  }

  for(std::map<int32_t,drop_count_t>::iterator it = t2cnt.begin(); it != t2cnt.end(); ++it) {
    drop_count_t& tCnt = it->second;    
    fprintf(fpT, "%s\t\t%d\t%d\t%.3lf\t%d\n", tid2name[it->first].c_str(), tCnt.umis, tCnt.reads, tCnt.fumis, tCnt.etc);              
  }

  bam_destroy1(b);
  
  fclose(fpCT);
  fclose(fpC);
  fclose(fpT);  

  return 0;
}

/*
int32_t runGeneForest(int32_t argc, char** argv) {
  std::string inVcf;
  std::string mapf;
  std::string reg;
  std::string outf;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Input Sites", NULL)
    LONG_STRING_PARAM("in-vcf",&inVcf, "Input VCF/BCF file")
    LONG_STRING_PARAM("map",&mapf, "Map file containing population of each individual")    
    LONG_STRING_PARAM("region",&reg,"Genomic region to focus on")

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out", &out, "Output VCF file name")
    LONG_INT_PARAM("verbose",&verbose,"Frequency of verbose output (1/n)")    
  END_LONG_PARAMS();  

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();


  // sanity check of input arguments
  if ( inVcf.empty() || mapf.empty() || out.empty() ) {
    error("--in-vcf, --map, --out are required parameters");
  }
  
  htsFile* fp = hts_open(mapf.c_str(),"r");
  if ( fp == NULL )
    error("Cannot open file %s for reading", mapf.c_str());

  kstring_t str = {0,0,0};
  int32_t lstr = 0;
  int32_t* flds = NULL;
  int32_t nflds = 0;

  std::map<std::string, std::string> id2pop;
  std::map<std::string, std::string> id2con;  

  notice("Started Reading sample map information");
  
  while( ( lstr = hts_getline(fp, KS_SEP_LINE, &str) ) >= 0 ) {
    flds = ksplit(&str, 0, &nflds); // expects ID, POPULATION, CONTINENT
    if ( nflds != 3 )
      error("Expected 3 columns from %s but observed %d", mapf.c_str(), nflds);

    id2pop[flds[0]] = flds[1];
    id2con[flds[0]] = flds[2];
  }

  notice("Started Reading site information from VCF file");

  std::vector<GenomeInterval> intervals;
  if ( !reg.empty() ) {
    parse_intervals(intervals, "", reg);
  }
  BCFOrderedReader odr(inVcf, intervals);
  bcf1_t* iv = bcf_init();
  
  uint32_t *fls = NULL, *frs = NULL, *rls = NULL, *rrs = NULL, *als = NULL;
  uint32_t nfl = 0, nfr = 0, nrl = 0, nrr = 0, nal = 0;

  int32_t nhaps = bcf_hdr_nsamples(odr.hdr);

  std::vector<int32_t> lois(nhaps,0);
  std::vector<int32_t> rois(nhaps,0);  
    
  while( odr.read(iv) ) {
    bcf_unpack(iv, BCF_UN_ALL);

    if ( bcf_get_format_int32(odr.hdr, iv, "RR", &rrs, &nrr) < 0 )
      error("Cannot parse RR field");

    if ( bcf_get_format_int32(odr.hdr, iv, "RL", &rrs, &nrr) < 0 )
      error("Cannot parse RL field");

    if ( bcf_get_format_int32(odr.hdr, iv, "FR", &rrs, &nrr) < 0 )
      error("Cannot parse FR field");

    if ( bcf_get_format_int32(odr.hdr, iv, "FL", &rrs, &nrr) < 0 )
      error("Cannot parse FL field");

    if ( bcf_get_format_int32(odr.hdr, iv, "AL", &rrs, &nrr) < 0 )
      error("Cannot parse AL field");

    // first, order the individuals based on ranks - should be linear time
    for(int32_t i=0; i < nhaps; ++i) {
      rois[rrs[i]] = i;
      lois[lrs[i]] = i;
    }

    // next, order by matching length - should be nlog(n)
  }
}
*/

int32_t main(int32_t argc, char** argv) {
  if ( argc < 2 ) {
    printf("cramore -- Fast analytic tools for analyzing and manipulating SAM/BAM/CRAM files\n");
    printf("Copyright (c) 2016 Hyun Min Kang and Adrian Tan\n");
    printf("Usage : %s [command] [options]\n",argv[0]);
    printf("\tType one of the following commands below to get detailed usage\n");
    printf("\t%s kallisto-count [options] : Count digital expressions from kallisto-aligned dropseq data\n",argv[0]);
    printf("\t%s sparse-genotypes [options] : Sparse genotyping from BAM files\n",argv[0]);
    printf("\t%s mendel-dup-conc [options] : Mendelian/Duplicate genotype concordance\n",argv[0]);
    printf("\t%s vcf-sample-summary [options] : Sample level summary of VCF file\n",argv[0]);
  }
  else {
    std::string cmd(argv[1]);
    if ( cmd == "kallisto-count" ) {
      return runKallistoCount(argc-1,argv+1);
    }
    else if ( cmd == "map-stamps" ) {
      return runMapSTAMPs(argc-1,argv+1);
    }    
    else if ( cmd == "verify-pair-id" ) {
      return runVerifyPairID(argc-1,argv+1);  
    }
    else if ( cmd == "sparse-genotype" ) {
      return runSparseGenotype(argc-1,argv+1);
    }
    else if ( cmd == "mendel-dup-conc") {
      return runMendelDupConc(argc-1,argv+1);      
    }
    else if ( cmd == "vcf-sample-summary") {
      return runVcfSampleSummary(argc-1,argv+1);      
    }
    else if ( cmd == "vcf-squeeze") {
      return runVcfSqueeze(argc-1,argv+1);      
    }
    else if ( cmd == "multinom-em") {
      return runClusterMultinomEM(argc-1, argv+1);
    }
    else if ( cmd == "simul-contam") {
      return runSimulContam(argc-1, argv+1);
    }    
    else {
      error("Unrecognized command %s\n",argv[1]);
    }
  }
  return 0;
}
