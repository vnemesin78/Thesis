#Application
	APP_DIR=App

	
all:
	@rm -f $(LOG_FILE)
	@(cd tkalman; make)
	@(cd image; make)
	@(cd config; make )
	@(cd iris;make )
	@(cd $(APP_DIR);make)
clean:
	@(cd tkalman; make clean)
	@(cd image; make clean)
	@(cd iris;make clean) 
	@(cd $(APP_DIR);make clean)
MrProper:
	@(cd tkalman; make MrProper)
	@(cd image; make MrProper)
	@(cd iris;make MrProper)
	@(cd $(APP_DIR);make MrProper)
forced:
	@make clean
	@make
