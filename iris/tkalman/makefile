	#Rep.
	REP=$(subst objects/ ,, $(subst tests/ , , $(subst include/ ,, $(wildcard ls -d */) ) ) )
	
	#make
	MAKE_ROOT = $(addprefix cd ,$(addsuffix ;make;, $(REP)) )
	MAKE_ALL = $(subst ;make;,;make;cd ..;,$(MAKE_ROOT))
	MAKE_MR_PROPER = $(subst ;make;,;make MrProper;cd ..;,$(MAKE_ROOT))
	MAKE_CLEAN = $(subst ;make;,;make clean;cd ..;,$(MAKE_ROOT))
	MAKE_CP = $(subst ;make;,;make cp;cd ..;,$(MAKE_ROOT))	
	MAKE_ADD= $(subst ;make;,;make add;cd ..;,$(MAKE_ROOT))	
#Rep courant	
	CUR=$(shell echo `pwd`)
all:
	@$(MAKE_ALL)
	
clean:
	@$(MAKE_CLEAN)
	
MrProper:
	@$(MAKE_MR_PROPER)

forced:
	@$(MAKE_CLEAN)
	@$(MAKE_ALL)

add:
	@cd $(D_PATH); git add $(CUR)/makefile
	@$(MAKE_ADD)

git-pull:
	git pull origin master;
	
git-push:
	git push origin master;
	
git-commit:
	git add -A;
	git commit -m"dtc"
