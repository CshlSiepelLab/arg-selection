### WORKFLOW TO RUN RELATE ON DISCOAL SIMULATIONS ###

1) Partition pickle files to enable parallel processing

	*Make sure the `handles.txt` file contains the following information in
two tab-delimited columns, each row corresponds to one pickle file
	<suffix>	<threads>
	where discoal_<suffix>.pkl is the pickle file name and <threads> is
the number of files each pickle file is partitioned into

	*Make sure the for loop in `partition.sh` loops through all lines
in `handles.txt`

	*Make sure the `partition.sh` file has the correct directory to the
pickle files, and also the desired name of the output directory

	Having `handles.txt`, `partition.sh` and
`partition.py` in the same directory
	RUN: `qsub partition.sh` in that directory

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

2) Run relate and extract features

	For each SET of pickle files (created by the partition step above), do
the following:

	*Make sure the number of array tasks (specified by #$ -t in
`sim2arg_multThr.sh`) match the number of partitioned pickle files

	*Make sure the desired version of feature extraction python script is
called in `sim2arg_multThr.sh`

	*Make sure the global var `RELATE_PATH` at the beginning of `sim2arg.py` is valid`

	Having `time.txt`, `sim2arg.py` (or another version of it) in the same
directory
	RUN: `qsub sim2arg_multThr.sh <suffix> <Ne> <identifier>`
	where <suffix> is the same as described above, Ne is the corresponding
pop size and <identifier> can be specified as you wish

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3) IF any thread fails due to memory issues, you will need to re-run the specific thread

	*Remove the temporary folders by running `rm -r */*_temp_*`, make sure no sim2arg jobs is running when doing this
	
	For each failed thread,
	RUN: `qsub sim2arg_patch.sh <suffix> <Ne> <failed_thread_#> <identifier>`

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4) Combine the extracted features

        usage: $./combine.py <handle> <TAG> <total threads>
            Combine discoal_<handle>/discoal_<handle>_<TAG>_inf_fea_*.pickle
                    discoal_<handle>_neutral/discoal_<handle>_neutral_<TAG>_inf_fea_*.pickle
            into inf_fea_<handle>_<TAG>.pkl
