
# benchmark_outputs/out_d10_explode1_gvcf1_rep1_qs0.log
# benchmark_outputs/out_d10_explode1_gvcf1_rep2_qs0.log
# benchmark_outputs/out_d10_explode1_gvcf1_rep3_qs0.log
# benchmark_outputs/out_d10_explode1_gvcf1_rep4_qs0.log
# benchmark_outputs/out_d10_explode1_gvcf1_rep5_qs0.log
# benchmark_outputs/out_d1_explode1_t1_rep1_qs0.log
# benchmark_outputs/out_d1_explode1_t1_rep1_qs25.log
# benchmark_outputs/out_d1_explode1_t1_rep2_qs0.log
# benchmark_outputs/out_d1_explode1_t1_rep2_qs25.log
# benchmark_outputs/out_d1_explode1_t1_rep3_qs0.log
# benchmark_outputs/out_d1_explode1_t1_rep3_qs25.log
# benchmark_outputs/out_d1_explode1_t1_rep4_qs0.log
# benchmark_outputs/out_d1_explode1_t1_rep4_qs25.log
# benchmark_outputs/out_d1_explode1_t1_rep5_qs0.log
# benchmark_outputs/out_d1_explode1_t1_rep5_qs25.log
# benchmark_outputs/out_d1_explode1_t8_rep1_qs0.log
# benchmark_outputs/out_d1_explode1_t8_rep1_qs25.log
# benchmark_outputs/out_d1_explode1_t8_rep2_qs0.log
# benchmark_outputs/out_d1_explode1_t8_rep2_qs25.log
# benchmark_outputs/out_d1_explode1_t8_rep3_qs0.log
# benchmark_outputs/out_d1_explode1_t8_rep3_qs25.log
# benchmark_outputs/out_d1_explode1_t8_rep4_qs0.log
# benchmark_outputs/out_d1_explode1_t8_rep4_qs25.log
# benchmark_outputs/out_d1_explode1_t8_rep5_qs0.log
# benchmark_outputs/out_d1_explode1_t8_rep5_qs25.log

# from these names get mean of replicates when the rest of the name is the same 




# ====================

# for i in benchmark_outputs/out_d1_explode1_t1_rep*_qs0.log;do

#     timems=$(cat $i | grep "Elapsed" | cut -d: -f5-)

#     echo $timems
# done | sort |cut -d':' -f1|datamash mean 1

# # -> mean result for replicates:
# format M:S
# # 45 minutes for no multithreading, d1, whole chr22, qs0


# ====================


#format H:M:S
# for i in benchmark_outputs/out_d1_explode1_t1_rep*_qs25.log;do

#     timems=$(cat $i | grep "Elapsed" | cut -d: -f5-)

#     echo $timems
# done | cut -d':' -f1 

# all 1 hr


# for i in benchmark_outputs/out_d1_explode1_t1_rep*_qs25.log;do

#     timems=$(cat $i | grep "Elapsed" | cut -d: -f5-)

#     echo $timems
# done | cut -d':' -f2 | datamash mean 1
# # 36

# -> mean result for replicates:
# 1:36 (1 hr 36 min) for no multithreading, d1, whole chr22, qs25


# ====================

#format: M:S

# for i in benchmark_outputs/out_d1_explode1_t8_rep*_qs0.log;do

#     timems=$(cat $i | grep "Elapsed" | cut -d: -f5-)

#     echo $timems
# done  | cut -d':' -f1 | datamash mean 1

# -> mean result for replicates:
# 11.4 minutes for 8 threads, d1, whole chr22, qs0

# 1 thread: 45 minutes -> 8 threads: 11.4 minutes

# ====================

#format: M:S

# for i in benchmark_outputs/out_d1_explode1_t8_rep*_qs25.log;do

#     timems=$(cat $i | grep "Elapsed" | cut -d: -f5-)

#     echo $timems
# done | cut -d':' -f1 | datamash mean 1

# -> mean result for replicates:
# 53 minutes for 8 threads, d1, whole chr22, qs25

# 1 thread: 1:36 (96 minutes) -> 8 threads: 53 minutes

# ====================


# format: H:M:S
# benchmark_outputs/out_d10_explode1_gvcf1_rep1_qs0.log

for i in benchmark_outputs/out_d10_explode1_gvcf1_rep*_qs0.log;do

    timems=$(cat $i | grep "Elapsed" | cut -d: -f5-)

    echo $timems
done  | cut -d':' -f1,2 | awk -F: '{sum1+=$1; sum2+=$2} END {print sum1/NR, sum2/NR}'


# -> mean result for replicates:
# 1:39 (1 hr 39 min) for no multithreading, d10, whole chr22, qs0, gvcf blocking





# #collect time -v output for runtime:
# for i in benchmark_outputs/out_d1*log; do 
#     # echo $i;

#     DP=$(echo $i|cut -d'/' -f2| cut -d_ -f2 | tr -d 'd')
#     REP=$(echo $i|cut -d'/' -f2| cut -d_ -f5|sed 's/rep//g')
#     # echo $DP
#     # echo $REP
    

#     timems=$(cat $i | grep "Elapsed" | cut -d: -f5-)

# done



# Latex:

The maximum runtime for these simulations

96 minutes with quality score errors

45 minutes without quality score errors 