PACKAGES := $(wildcard *.m)
PACKAGEDIR := $(HOME)
BASEVER := c8cbfe9

PACKAGEDIR_FILES := $(foreach file,$(PACKAGES),$(PACKAGEDIR)/$(file))

all : $(PACKAGEDIR_FILES)

$(PACKAGEDIR)/%.m : %.m
	cp $< $@
