PACKAGES := $(wildcard *.m)
PACKAGEDIR := $(HOME)

PACKAGEDIR_FILES := $(foreach file,$(PACKAGES),$(PACKAGEDIR)/$(file))

all : $(PACKAGEDIR_FILES)

$(PACKAGEDIR)/%.m : %.m
	cp $< $@
