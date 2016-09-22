#include "sex_ploidy_map.h"

int32_t sex_ploidy_map::get_ploidy_type(bcf1_t* v) {
  if ( v ) {
    if ( v->rid == x_rid ) {
      if ( ( v->pos >= x_start ) && ( v->pos <= x_stop ) ) {
	return PLOIDY_TYPE_X;
      }
      else
	return PLOIDY_TYPE_AUTOSOME;
    }
    else if ( v->rid == y_rid )
      return PLOIDY_TYPE_Y;
    else if ( v->rid == mt_rid )
      return PLOIDY_TYPE_MT;
    else
      return PLOIDY_TYPE_AUTOSOME;
  }
  else
    return PLOIDY_TYPE_AUTOSOME;
}

int32_t sex_ploidy_map::set_rids_from_hdr(bcf_hdr_t* hdr) {
  x_rid = bcf_hdr_name2id(hdr, x_label.c_str());
  y_rid = bcf_hdr_name2id(hdr, y_label.c_str());
  mt_rid = bcf_hdr_name2id(hdr, mt_label.c_str());

  return ( ( x_rid < 0 ? 0 : 1 ) + ( y_rid < 0 ? 0 : 1 ) + ( mt_rid < 0 ? 0 : 1 ) );
}

int32_t sex_ploidy_map::load_sex_map_file(const char* filename, bcf_hdr_t* hdr) {
  set_rids_from_hdr(hdr);
  int32_t nsamples = bcf_hdr_nsamples(hdr);

  //fprintf(stderr,"load_sex_map_files() : nsamples=%d\n",nsamples);
  
  if ( ( filename == NULL ) || ( strlen(filename) == 0 ) ) {
    for(int32_t i=0; i < nsamples; ++i) {
      vSex.push_back(2);
    }
  }
  else {
    htsFile *file = hts_open(filename,"r");
    if ( file == NULL ) {
      fprintf(stderr,"ERROR: Cannot open %s\n",filename);
      exit(1);
    }
    
    kstring_t *s = &file->line;
    while( hts_getline(file,'\n',s) >= 0 ) {
      std::string ss = std::string(s->s);
      size_t idx = ss.find_first_of("\t ");
      if ( idx == std::string::npos ) {
	fprintf(stderr,"ERROR: Cannot parse line %s in %s\n",ss.c_str(), filename);
	exit(1);
      }
      std::string id = ss.substr(0, idx);
      int32_t sex = atoi(ss.substr(idx+1).c_str());
      
      if ( mSex.find(id) != mSex.end() ) {
	fprintf(stderr,"ERROR: Duplicate ID %s in %s\n",id.c_str(), filename);
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

    for(int32_t i=0; i < nsamples; ++i) {
      const char* sid = bcf_hdr_int2id(hdr, BCF_DT_SAMPLE, i);
      std::map<std::string,int>::iterator it = mSex.find(sid);
      if ( it == mSex.end() ) { // not found
	fprintf(stderr,"WARNING: No sex information is available for %s, treating as female\n",sid);
	vSex.push_back(2);
      }
      else {
	vSex.push_back(it->second);
	if ( it->second == 1 )
	  ++n_males;
      }
    }
  }
    
  return nsamples;
}

