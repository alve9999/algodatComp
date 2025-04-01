
all:
	mkdir -p ../data
	@echo "generating data 10k"
	cargo r -r -- 10000 > 10k.txt
	@echo "generating data 1M"
	cargo r -r -- 1000000 > 1M.txt
	@echo "generating data 10M"
	cargo r -r -- 10000000 > 10M.txt
	cp *.txt ../data
	@echo "Done!"
