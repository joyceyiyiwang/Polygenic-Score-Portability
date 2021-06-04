import os 
import pathlib

#Script used to manually submit jobs, just modify range values
start = 1
end = 10 + 1
temp_dir = pathlib.Path('temp_fst_path/')
for group_label in range(start,end):
    script_path = temp_dir.joinpath(f'run_{group_label}.sh')
    os.system(f'sbatch {script_path.as_posix()}')
