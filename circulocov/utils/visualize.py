#!/usr/bin/env python

''' visualize top 10 contigs '''

import concurrent.futures
import logging

from .top_ten import top_ten
from .create_figure import create_circos_figure
from .create_figure import create_linear_figure

def visualize (genome_dict, df, args):

    ''' visualize top 10 contigs '''

    len_threshold = top_ten(genome_dict)

    with concurrent.futures.ThreadPoolExecutor(max_workers = args.threads) as executor:
        tasks = []
        for contig in genome_dict.keys():
            if genome_dict[contig]['length'] >= len_threshold:
                logging.info(f"Creating depth figure for {genome_dict[contig]['name']}")
                if genome_dict[contig]['circ'] :
                    future = executor.submit(create_circos_figure, genome_dict[contig], df, args)
                    tasks.append(future)
                else:
                    future = executor.submit(create_linear_figure, genome_dict[contig], df, args)
                    tasks.append(future)

        ct, _ = concurrent.futures.wait(tasks,return_when=concurrent.futures.ALL_COMPLETED)

        logging.debug(ct)
