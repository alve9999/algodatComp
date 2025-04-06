
all:
	@echo "TESTING THE AI C FILES"
	@for file in *.c; do \
		out="$${file}.out"; \
		echo " ------- COMPILING $$file to $$out ------------"; \
		gcc ../parse_and_run.c $$file -o $$out -O3 -g -fopenmp; \
		for data in ../data/*.txt; do \
		    echo "running $$data"; \
			./$$out < $$data; \
		done; \
	done
