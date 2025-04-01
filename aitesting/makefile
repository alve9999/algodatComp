
all:
	@echo "testing the ai c files"
	@for file in *.c; do \
		out="$${file}.out"; \
		echo " ------- COMPILING $$file to $$out ------------"; \
		gcc $$file -o $$out -O3 -g -fopenmp; \
		for data in ../data/*; do \
		    echo "running $$data"; \
			time  ./$$out < $$data; \
		done; \
	done
