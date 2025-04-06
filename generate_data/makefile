
all:
	mkdir -p ../data
	@echo " ::: Generating 1M files :::"
	@for i in A B C ; do \
		name=$${i}_1M.txt ; \
		echo "generating $$name" ; \
		cargo r -r -- 1000000 > $$name ; \
		echo "copying $$name to ../data" ; \
		cp $$name ../data ; \
	done
	@echo "Done!"
