OPTFLAG = -g -O3 -fopenmp
INCLUDES = -I. -I../htslib -I../htslib/htslib
CFLAGS = -pipe -std=c++0x $(OPTFLAG) $(INCLUDES) -D__STDC_LIMIT_MACROS
CXX = g++

SOURCES = Error PhredHelper params dropseq hts_utils bcf_ordered_reader interval interval_tree utils genome_interval log_tool estimator bcf_ordered_writer nuclear_pedigree genotype_concordance filter

SOURCESONLY = cramore.cpp

TARGET = cramore
TOOLSRC = $(SOURCES:=.cpp) $(SOURCESONLY)
TOOLOBJ = $(TOOLSRC:.cpp=.o)
LIBHTS = ../htslib/libhts.a

all : $(TARGET)

${LIBHTS} :
	cd ../htslib; $(MAKE) || exit 1; 

$(TARGET) : ${LIBHTS} $(TOOLOBJ)
	$(CXX) $(CFLAGS) -o $@ $(TOOLOBJ) $(LIBHTS) -lz -lpthread -fopenmp

$(TOOLOBJ): $(HEADERSONLY)

.cpp.o :
	$(CXX) $(CFLAGS) -o $@ -c $*.cpp

.PHONY: clean cleanvt test

clean :
	cd ../htslib; $(MAKE) clean
	rm -rf $(TARGET) $(TOOLOBJ)

test: vt
	for x in ./test/*/_run.sh; do \
		$${x}; \
	done
