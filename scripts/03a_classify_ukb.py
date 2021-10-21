import pickle

import pandas as pd


def classify_ukbb(df):
    """
    Classifies whether a participant is WB based on UKBB data field 22006
    """
    # Ancestry data
    wb_df = pd.read_csv(
        'data/extracted_phenotypes/ukb.wb.id.txt',
        sep=' ',
        header=0
    )
    wb_df = wb_df.loc[:, ['IID']]

    wb_df = wb_df.assign(WB = ['WB'] * wb_df.shape[0])

    df = df.merge(wb_df, on='IID', how='left')
    df = df.fillna('NWB')
    
    return df


if __name__ == '__main__':

    labels_df = (
        pd.read_csv('data/ukb_merged/merged.psam', sep='\t')
        #.drop_duplicates(subset=['#FID', 'IID'])
    )
    labels_df = labels_df.loc[:, ['#FID', 'IID']]

    ukb_labels_df = classify_ukbb(labels_df)
    ukb_labels_df.to_csv('data/ukb_merged/population_labels.tsv.gz',
                         sep='\t', compression='infer', index=False)
