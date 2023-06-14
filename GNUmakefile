PACKAGES := $(wildcard *.m)

# MacOS:
PACKAGEDIR := $(HOME)

# $UserBaseDirectory as viewed from WSL2 bash on Windows:
#PACKAGEDIR := /mnt/c/Users/peete
#PACKAGEDIR := /mnt/c/Users/peete/AppData/Roaming/Mathematica
#PACKAGEDIR := /mnt/c/ProgramData/Mathematica

#"C:\\Users\\peete\\AppData\\Roaming\\Mathematica\\DocumentationIndices"
#"C:\\Program Files\\Wolfram\ Research\\Mathematica\\12.2\\SystemFiles\\Links"
#"C:\\Users\\peete\\AppData\\Roaming\\Mathematica\\Kernel"
#"C:\\Users\\peete\\AppData\\Roaming\\Mathematica\\Autoload"
#"C:\\Users\\peete\\AppData\\Roaming\\Mathematica\\Applications"
#"C:\\ProgramData\\Mathematica\\Kernel"
#"C:\\ProgramData\\Mathematica\\Autoload"
#"C:\\ProgramData\\Mathematica\\Applications"
#"C:\\Users\\peete"
#"C:\\Program Files\\Wolfram Research\\Mathematica\\12.2\\AddOns\\Packages"

PACKAGEDIR_FILES := $(foreach file,$(PACKAGES),$(PACKAGEDIR)/$(file))

all : $(PACKAGEDIR_FILES)

$(PACKAGEDIR)/%.m : %.m
	cp $< $@

clean:
	rm -f $(PACKAGEDIR_FILES)
