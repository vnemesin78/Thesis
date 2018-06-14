TARGET=iris_default.hpp

APP=~/Applications/convert_cfg_2_hpp.run

CONFIG_DIR=~/Config ~/Config/image

CONFIG=$(wildcard $(addsuffix /*.cfg, $(CONFIG_DIR)))

#Rep courant	
	CUR=$(shell echo `pwd`)

$(TARGET) : $(CONFIG)
	$(APP) $(TARGET) $(CONFIG) ~/Config/fusion/lin.cfg ~/Config/matching/hamming.cfg

add :
	@(cd ..; git add $(CUR)/$(TARGET) $(CUR)/makefile)

clean:
	@rm -f $(TARGET)
	
forced:
	@make clean;
	@make

git-pull:
	git pull origin master;
	
git-push:
	git push origin master;
	
git-commit:
	git add -A;
	git commit -m"dtc"
