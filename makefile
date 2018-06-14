all:
	@(cd API; make)
	@(cd tkalman_c;make)
	@(cd iris;make)
clean:
	@(cd API; make clean)
	@(cd tkalman_c;make clean)
	@(cd iris;make clean)
MrProper:
	@(cd API; make MrProper)
	@(cd tkalman_c;make MrProper)
	@(cd iris;make MrProper)
forced:
	make clean
	make all
