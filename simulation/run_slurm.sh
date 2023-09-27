
logdir="sim/logs/slurm_logs/logs_"$(date +'%y%m%d_%H%M%S')
mkdir -pv $logdir

snakemake --cluster-config cluster.json --cluster "sbatch --partition=compregular --nodes={cluster.n} --tasks-per-node={cluster.tasks} --mem={cluster.mem_mb} --job-name={cluster.name} --output=$logdir/{cluster.logfile} --error=$logdir/{cluster.errfile} -t 1-03:00:00" --jobs 200 --notemp --nolock --latency-wait 60 --rerun-incomplete --immediate-submit -s ${1}


