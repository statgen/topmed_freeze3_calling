SUBDIRS = htslib samtools bcftools cramore gotcloud vt

.PHONY: subdirs $(SUBDIRS)

subdirs: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@

SUBCLEAN = $(addsuffix .clean,$(SUBDIRS))

.PHONY: clean $(SUBCLEAN)

clean: $(SUBCLEAN)
	find . -name '*.o' | xargs -I {} rm -f {}

$(SUBCLEAN): %.clean:
	$(MAKE) -C $* clean


