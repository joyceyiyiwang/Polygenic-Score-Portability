import itertools
import os
import pathlib

import pandas as pd
from random import sample


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return itertools.zip_longest(*args, fillvalue=fillvalue)


def submit_batch_fst_computation(group_label, test_fids, temp_dir):
    """
    Batch Fst computations. We have ~32,000 test non-EUR individuals and would rather
    submit O(100) jobs than 32,000 jobs.
    """
    formatted_fids = ' '.join([str(i) for i in test_fids])
    script_path = temp_dir.joinpath(f'run_{group_label}.sh')
    if group_label % 5 == 0:
        script_string = (
            "#!/bin/bash\n"
            f"#SBATCH -J fst{group_label}\n"
            f"#SBATCH -o fst{group_label}.o%j\n"
            f"#SBATCH -e fst{group_label}.o%j\n"
            "#SBATCH -p normal\n"
            "#SBATCH -N 1\n"
            "#SBATCH -n 1\n"
            "#SBATCH -t 4:00:00\n"
            "#SBATCH -A Harpak-Lab-GWAS\n"
            "#SBATCH --mail-user=joyce.wang@utexas.edu\n"
            "#SBATCH --mail-type=begin\n"
            "#SBATCH --mail-type=end\n\n"
            "plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'\n"
            f"group={group_label}\n"
            f"for test_fid in {formatted_fids}\n"
            "do\n"
            "  $plink "
            "  --bim ../data/ukb_merged/merged.bim "
            "  --bed ../data/ukb_merged/merged.bed "
            "  --fam ../data/ukb_merged/merged_fst.fam "
            "  --family "
            "  --fst "
            "  --keep-cluster-names train ${test_fid} "
            "  --out ../data/fst_full/fst${test_fid}\n\n"
            "  rm ../data/fst_full/fst${test_fid}.fst ../data/fst_full/fst${test_fid}.nosex\n"
            "  echo ${test_fid} >> ../data/fst_full/fst${group}.est\n"
            "  grep 'Fst estimate:' ../data/fst_full/fst${test_fid}.log >> ../data/fst_full/fst${group}.est\n"
            "  rm ../data/fst_full/fst${test_fid}.log\n"
            "done\n"
            # Delete the script itself
            )
        with open(script_path, 'w') as f:
            f.write(script_string)
#        os.system(f'sbatch {script_path.as_posix()}')
    else:
        script_string = (
            "#!/bin/bash\n"
            f"#SBATCH -J fst{group_label}\n"
            f"#SBATCH -o fst{group_label}.o%j\n"
            f"#SBATCH -e fst{group_label}.o%j\n"
            "#SBATCH -p normal\n"
            "#SBATCH -N 1\n"
            "#SBATCH -n 1\n"
            "#SBATCH -t 4:00:00\n"
            "#SBATCH -A Harpak-Lab-GWAS\n\n"
            "plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'\n"
            f"group={group_label}\n"
            f"for test_fid in {formatted_fids}\n"
            "do\n"
            "  $plink "
            "  --bim ../data/ukb_merged/merged.bim "
            "  --bed ../data/ukb_merged/merged.bed "
            "  --fam ../data/ukb_merged/merged_fst.fam "
            "  --family "
            "  --fst "
            "  --keep-cluster-names train ${test_fid} "
            "  --out ../data/fst_full/fst${test_fid}\n\n"
            "  rm ../data/fst_full/fst${test_fid}.fst ../data/fst_full/fst${test_fid}.nosex\n"
            "  echo ${test_fid} >> ../data/fst_full/fst${group}.est\n"
            "  grep 'Fst estimate:' ../data/fst_full/fst${test_fid}.log >> ../data/fst_full/fst${group}.est\n"
            "  rm ../data/fst_full/fst${test_fid}.log\n"
            "done\n"
            # Delete the script itself
            )
        with open(script_path, 'w') as f:
            f.write(script_string)
#        os.system(f'sbatch {script_path.as_posix()}')
    


def main():
    temp_dir = pathlib.Path('temp_fst_path_full/')
    temp_dir.mkdir(exist_ok=True)

    all_test = (
        pd.read_csv('data/ukb_merged/merged_fst.fam', sep='\t')
        .loc[lambda df: df['#FID'] != 'train', '#FID']
        .values
        .tolist()
    )

    processed = (
        pd.read_csv('data/fst/final_fst_full.tsv', sep='\t')
        .[['IID']]
        .values
        .tolist()
    )

    all_test = all_test - processed
    all_test = sample(all_test, 13973-7665)

    group_size = 40
    for i, test_fid_group in enumerate(grouper(all_test, group_size, '')):
        submit_batch_fst_computation(i + 192, test_fid_group, temp_dir)


if __name__ == "__main__":
    main()
