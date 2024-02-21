#!/usr/bin/env python
# pylint: disable=logging-fstring-interpolation

''' Gets coverage for bam '''

import logging
import pysam
import pandas as pd

def coverage(bam):

    ''' Gets coverage for bam '''

    logging.debug(f'Getting coverage for {bam}')
    covs = pysam.coverage(bam).split('\n')

    header = covs[0].split('\t')

    df = pd.DataFrame(columns = header)
    for cov in covs:
        cov = cov.split('\t')
        if '#rname' not in cov and len(cov) == len(header):
            df.loc[len(df.index)] = cov

    df['endpos']    = pd.to_numeric(df['endpos'], downcast='integer')
    df['meandepth'] = pd.to_numeric(df['meandepth'])

    df = df.sort_values(by=['endpos', '#rname'], ascending= [False, True], ignore_index=True)

    logging.debug(f'The coverage dataframe for {bam}')
    logging.debug(df)

    return df
