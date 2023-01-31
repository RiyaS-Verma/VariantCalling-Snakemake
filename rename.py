#code to rename files to generate matching trio sample key used in the snakemake analysis
#assuming proband, father and mother have unique IDs, run this to generate file names with trio identifying keys

import os
import pandas as pd
import subprocess
from glob import glob
from pathlib import Path
from numpy import unique

df = pd.read_csv("/path/to/csv/trionames.csv") #sample trionames.csv in main file
WORK_DIR='/path/to/working/directory/'
for index, row in df.iterrows():
    trio = row['Trio#']
    probandId = row['Proband_ID']
    fatherId = row['Father_ID']
    motherId = row['Mother_ID']

    proband_dir = WORK_DIR + '/readsfolder/' + str(probandId)
    for filename in os.listdir(proband_dir):
        readnum=filename.split('_',)[1]
        os.rename(proband_dir + '/' + filename, WORK_DIR + '/' + "Trio{}_P_{}.fastq.gz".format(trio, readnum))

    # rename father files
    os.rename("{0}/{1}_1.fastq.gz".format(WORK_DIR, fatherId), "{0}/Trio{2}_F_1.fastq.gz".format(WORK_DIR, fatherId, trio))
    os.rename("{0}/{1}_2.fastq.gz".format(WORK_DIR, fatherId), "{0}/Trio{2}_F_2.fastq.gz".format(WORK_DIR, fatherId, trio))
    
    # rename mother files
    os.rename("{0}/{1}_M_1.fastq.gz".format(WORK_DIR, motherId), "{0}/Trio{2}_M_1.fastq.gz".format(WORK_DIR, motherId, trio))
    os.rename("{0}/{1}_M_2.fastq.gz".format(WORK_DIR, motherId), "{0}/Trio{2}_M_2.fastq.gz".format(WORK_DIR, motherId, trio))
