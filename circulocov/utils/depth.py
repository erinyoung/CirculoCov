#!/usr/bin/env python
# pylint: disable=logging-fstring-interpolation

""" Gets depth for bam """

import logging
import pysam
import pandas as pd

def depth(bam, reg):
    """ Gets coverage for bam over subrange """

    depths = pysam.depth(bam, "-r", reg).split("\n")

    df = pd.DataFrame(columns = ["ref", "pos", "depth"])
    for dep in depths:
        dept = dep.split("\t")
        if len(dept) == 3:
            df.loc[len(df.index)] = dept

    logging.debug(f"The depth dataframe for {bam}")
    logging.debug(df)

    if df.empty:
        ave = 0.0
    else:
        df["depth"] = pd.to_numeric(df["depth"])

        ave = df["depth"].mean()

    logging.debug(f"The average for {reg} is {ave}")

    return reg, ave
