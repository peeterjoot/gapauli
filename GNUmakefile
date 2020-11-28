PACKAGES := $(wildcard *.m)

# MacOS:
#PACKAGEDIR := $(HOME)

# $UserBaseDirectory as viewed from WSL2 bash on Windows:
#PACKAGEDIR := /mnt/c/Users/peete/AppData/Roaming/Mathematica
PACKAGEDIR := /mnt/c/ProgramData/Mathematica

PACKAGEDIR_FILES := $(foreach file,$(PACKAGES),$(PACKAGEDIR)/$(file))

all : $(PACKAGEDIR_FILES)

$(PACKAGEDIR)/%.m : %.m
	cp $< $@
