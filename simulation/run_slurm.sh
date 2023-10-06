projectID="vcfgl"

configFile=${projectID}_cluster_config.yaml

logdir="sim/slurm_logs/logs_"$(date +'%y%m%d_%H%M%S')

echo "{
	\"__default__\" :
	{
		\"mem_mb\" : 4000,
		\"n\" : 1,
		\"tasks\" : 1,
		\"name\" : \"PROJECT_${projectID}-RULE_{rule}-JOB_{jobid}\",
		\"logfile\": \"${logdir}/smk_{rule}_%j.out\",
		\"errfile\": \"${logdir}/smk_{rule}_%j.err\",
		\"time\": \"13:00:00\"
	}
}" > ${configFile}

# \"time\": \"1-03:00:00\"

mkdir -pv ${logdir}

snakemake --cluster-config ${configFile} \
	--cluster "sbatch --partition=compregular --nodes={cluster.n} --tasks-per-node={cluster.tasks} --mem={cluster.mem_mb} --job-name={cluster.name} --output={cluster.logfile} --error={cluster.errfile} -t {cluster.time} " \
	--jobs 5000 \
	--notemp \
	--nolock \
	--immediate-submit \
	--latency-wait 60 \
	--rerun-incomplete \
	-s "${@}"

#231002/10:48 disabled this for now to see if slurm can work with multiple run alls when this is disabled. enabled: it cant
# --immediate-submit \

