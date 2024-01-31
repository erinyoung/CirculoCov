#!/usr/bin/env python

""" Gets coverage for sam """

import concurrent.futures
import pandas as pd

from .coverage import coverage
from .depth import depth
from .regions import regions, padded_regions
from .create_dataframe import create_depth_dataframe

def counts(bam, genome_dict, analysis, args):
    """ Gets relevant counts """

    df_cov   = coverage(bam)
    cov_cols = ["#rname",
                "startpos",
                "endpos",
                "numreads",
                "covbases",
                "coverage",
                "meandepth",
                "meanbaseq",
                "meanmapq"]

    df_cov.to_csv(  args.out + '/' + analysis + '_cov.txt'  , columns=cov_cols, index=False, sep = "\t")

    if args.all:
        # dividing up the genome into windows for coverage
        regs  = regions(genome_dict, args)
        pregs = padded_regions(genome_dict, args)

        # TODO : remove these once finished
        regs=['5:2501-3000', '5:3001-3500', '5:3501-4000', '5:4001-4500']
        pregs=['5:11019-11518', '5:11519-12018', '5:12019-12518']

        results = []
        padded_results = []
        # to speed up processing
        with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
            tasks  = []
            ptasks = []
            for reg in regs:
                future = executor.submit(depth, bam, reg)
                tasks.append(future)

            for preg in pregs:
                pfuture = executor.submit(depth, bam, preg)
                ptasks.append(pfuture)

            completed_tasks, _ = concurrent.futures.wait(tasks,return_when=concurrent.futures.ALL_COMPLETED)
            results = [task.result() for task in completed_tasks]

            pcompleted_tasks, _ = concurrent.futures.wait(ptasks,return_when=concurrent.futures.ALL_COMPLETED)
            padded_results = [ptask.result() for ptask in pcompleted_tasks]

        df_depth = create_depth_dataframe(results, padded_results, genome_dict)

        # creating csv of results
        dep_cols = ["contig", "position", "end", "mean_depth"]
        
        df_depth.to_csv(args.out + '/' + analysis + '_depth.txt', columns=dep_cols, index=False, sep = "\t")
    else:
        df_depth = pd.DataFrame()

    return df_depth, df_cov
    
    
    


