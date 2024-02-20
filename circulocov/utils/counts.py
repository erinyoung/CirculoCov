#!/usr/bin/env python
# pylint: disable=R0912
# pylint: disable=R0915

''' Gets coverage for sam '''

import concurrent.futures
import os
import pandas as pd
import logging

from .coverage import coverage
from .depth import depth
from .create_dataframe import create_depth_dataframe

def counts(bam, genome_dict, analysis, args, tmp):
    ''' Gets relevant counts '''

    logging.info(f'Getting coverage for {analysis}')

    df_cov   = coverage(bam)
    cov_cols = ['#rname',
                'startpos',
                'endpos',
                'numreads',
                'covbases',
                'coverage',
                'meandepth',
                'meanbaseq',
                'meanmapq']

    df_cov.to_csv(args.out + '/' + analysis + '_cov.txt',
                  columns=cov_cols,
                  index=False,
                  sep = '\t')

    if args.all:
        logging.info(f'Getting depth for {analysis}')

        # run by contig to speed up processing
        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
            tasks = []
            for contig in genome_dict.keys():
                future = executor.submit(depth, bam, contig, tmp)
                tasks.append(future)

            completed_tasks, _ = concurrent.futures.wait(tasks,return_when=concurrent.futures.ALL_COMPLETED)
            results = [task.result() for task in completed_tasks]

        df_depth = pd.DataFrame(columns = ['contig', 'pos', 'depth'])
        for dep_contig in results:
            if os.stat(dep_contig).st_size > 0:
                df_dep_contig = pd.read_table(dep_contig, names=['contig','pos', 'depth'])
                df_depth = pd.concat([df_depth, df_dep_contig], ignore_index=True)

        df_depth = df_depth.sort_values(['contig', 'pos']).reset_index(drop=True)

        df_depth.to_csv(args.out + '/' + analysis + '_full_depth.txt', index=False, sep = '\t')

        df_window_depth = create_depth_dataframe(df_depth, genome_dict, args)

        df_window_depth.to_csv(args.out + '/' + analysis + '_window_depth.txt', index=False, sep = '\t')

    else:
        df_window_depth = pd.DataFrame()

    return df_window_depth, df_cov
