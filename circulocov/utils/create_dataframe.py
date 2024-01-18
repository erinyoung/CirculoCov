#!/usr/bin/env python

""" Get dataframe of coverage """

import logging
import pandas as pd

def create_depth_dataframe(depths, pdepths, genome_dict):
    """ Creating dataframe of coverages """

    logging.info("Creating dataframe of coverages")
    df = pd.DataFrame(columns = ["contig", "position", "end", "depth"])

    for reg, depth in depths:
        contig = reg.split(":")[0]
        pos    = reg.split(":")[1].split("-")[0]
        end    = reg.split(":")[1].split("-")[1]
        df.loc[len(df.index)] = [contig, pos, end, depth]

    df["position"] = pd.to_numeric(df["position"], downcast="integer")
    df["depth"] = pd.to_numeric(df["depth"])
    df["match"]    = df["contig"] + ":" + df["position"].astype("string")
    df = df.sort_values(by=['contig', 'position'], ascending= [True, True], ignore_index=True)

    # now for the tricky part
    pdf = pd.DataFrame(columns = ["contig", "position", "pdepth"])

    for reg, depth in pdepths:
        contig = reg.split(":")[0]
        pos    = reg.split(":")[1].split("-")[0]
        pdf.loc[len(pdf.index)] = [contig, pos, depth]
        pdf["position"] = pd.to_numeric(pdf["position"], downcast="integer")
        pdf = pdf.sort_values(by=['contig', 'position'], ascending= [True, True], ignore_index=True)

    pdf["length"] = 0
    for contig in genome_dict.keys():
        pdf.loc[pdf["contig"] == contig, "length"] = genome_dict[contig]["length"]

    pdf["pdepth"] = pd.to_numeric(pdf["pdepth"])
    pdf["original"] = pdf["position"] - pdf["length"]
    pdf["match"]    = pdf["contig"] + ":" + pdf["original"].astype("string")
    pdf = pdf[["match", "pdepth"]]

    df = pd.merge(df, pdf, left_on="match", right_on="match", how="left")
    df["pdepth"]     = df[["pdepth"]].fillna(0)
    df["mean_depth"] = df["depth"] + df["pdepth"]

    logging.debug("Coverage dataframe")
    logging.debug(pdf)
    logging.debug(df)

    df = df[["contig", "position", "end", "match", "mean_depth"]]

    return df
