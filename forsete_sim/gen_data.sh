#!/bin/bash

echo "GENERATING DATA FOR FORSETE SIM."

mkdir -p data

for i in {1..100}
do
	echo "<< processing number $i >> "
	name=${i}_1M.txt
	echo "     creating $name "
	cargo r -r --manifest-path ../generate_data/Cargo.toml 1000000 > data/$name
done
