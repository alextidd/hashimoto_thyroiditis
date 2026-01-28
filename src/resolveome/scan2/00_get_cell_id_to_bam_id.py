#!/usr/bin/env python3

# packages
import pysam
import pandas as pd

# function: get sample name from header
def sample_from_bam(bampath):
    with pysam.AlignmentFile(bampath, 'r') as af:
        sample = af.header.to_dict()['RG'][0]['SM']
    print('Got sample name "%s" from BAM %s' % (sample, bampath))
    return sample

# load samplesheet
ss_path = "/lustre/scratch125/casm/teams/team268/at31/projects/hashimoto_thyroiditis/data/bams/samplesheet_local.csv"
ss = pd.read_csv(ss_path)

# get bam_id for each sample
ss['bam_id'] = ss['bam'].apply(sample_from_bam)

# create dataframe with id and bam_id
df = ss[['id', 'cell_id', 'bam_id']].copy()
print(df)

# save
out_path = "out/resolveome/scan2/id_to_bam_id.csv"
df.to_csv(out_path, index = False)
