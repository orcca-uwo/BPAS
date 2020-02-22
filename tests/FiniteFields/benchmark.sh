#!/bin/bash

output_dir="./benchmark_results"

suffix=""
if [[ $# -eq 1 ]]; then
	suffix="_$1"	
fi

echo $output_dir
#######################################
old_ifs="$IFS"
IFS=" "
cpu_model=( $(cat /proc/cpuinfo  | grep 'name'| uniq) )
cpu_model="${cpu_model[4]}""_""${cpu_model[7]}"
cpu_model="${cpu_model/-/_}"
cpu_model="${cpu_model/ /_}"

host_name=$(hostname)
host_name="${host_name/-/_}"
host_name="${host_name/ /_}"

output_name=""
output_name+="$host_name""_"
output_name+="$cpu_model"""
output_name+="$suffix"

#output_name+=".bnch"
echo $output_name
#######################################
node_control="numactl -N 0 "
bin="gfpf_six_step_fft_test.bin"

K_list=(8 16 32 64 128 256);
e_list=(2 3);
output="$output_dir""/gfpf_$output_name.bnch"

cat /dev/null > $output

for K in ${K_list[@]}; do	
	for e in ${e_list[@]}; do	
		if [[ $e -gt 2 ]] && [[ $K -gt 64 ]]; then
			printf "[SKIP K=$K and e=$e] ...\n";
			continue;
		fi;
		output_tmp="$output_dir""/gfpf_$output_name""_K"$K"_e"$e".bnch"
		printf "testing [$K ^ $e]...\n" > $output_tmp
		$node_control ./$bin $K $e >> $output_tmp
		printf "\n===================================\n" >>$output_tmp
		printf "===================================\n" >>$output_tmp
		cat $output_tmp >> $output
		cat $output_tmp
	done;
done;
