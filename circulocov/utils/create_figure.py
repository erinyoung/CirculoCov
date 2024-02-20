#!/usr/bin/env python

''' Visualize coverage'''

import pycirclize
import pandas as pd
import matplotlib.pyplot as plt

def create_circos_figure(contig_dict, df, args):
    ''' Visualize coverage for circular seqs '''

    df['pos'] = pd.to_numeric(df['pos'])
    chr_df = df[df['contig'].astype(str) == contig_dict['name']].copy()

    if not chr_df.empty:
        circos = pycirclize.Circos(sectors={contig_dict['name']: contig_dict['length']})
        circos.text(f"{args.sample} {contig_dict['name']}", size=15)

        for sector in circos.sectors:
            track1 = sector.add_track((50, 54))
            track1.axis(fc='black')
            major_interval = 250000
            minor_interval = int(major_interval / 10)
            if sector.size > minor_interval:
                track1.xticks_by_interval(major_interval, label_formatter=lambda v: f'{v / 1000:.0f} Kb') # pylint: disable=C0301
                track1.xticks_by_interval(minor_interval, tick_length=1, show_label=False)

            x = chr_df['pos'].to_numpy()
            i = 58

            if 'nanopore_depth' in df.columns:
                y      = chr_df['nanopore_depth'].to_numpy()
                track2 = sector.add_track((i, i + 10), r_pad_ratio=0.1)
                track2.fill_between(x, y, color='#6699CC')
                track2.axis()
                track2.grid()
                i = i + 11

            if 'illumina_depth' in df.columns:
                y      = chr_df['illumina_depth'].to_numpy()
                track3 = sector.add_track((i, i + 10), r_pad_ratio=0.1)
                track3.fill_between(x, y, color='#EECC66')
                track3.axis()
                track3.grid()
                i = i + 11

            if 'pacbio_depth' in df.columns:
                y      = chr_df['pacbio_depth'].to_numpy()
                track4 = sector.add_track((i, i + 10), r_pad_ratio=0.1)
                track4.fill_between(x, y, color='#EE99AA')
                track4.axis()
                track4.grid()

            if sector.name == circos.sectors[0].name:
                if 'nanopore_depth' in df.columns:
                    circos.text('nanopore', r=track2.r_center, deg = 10)
                if 'illumina_depth' in df.columns:
                    circos.text('illumina', r=track3.r_center, deg = 10)
                if 'pacbio_depth' in df.columns:
                    circos.text('pacbio'  , r=track4.r_center, deg = 10)

            # Future plans : add line for average coverage
            # track.line([track.start, track.end], [vmin, vmax], lw=1, ls='dotted')

        circos.savefig(f"{args.out}/{args.sample}_{contig_dict['name']}.png")

    return args.out + '/' + args.sample + '_' + contig_dict['name'] + '.png'

def create_linear_figure(contig_dict, df, args):

    ''' Visualize coverage for non-closed seqs '''

    df['pos'] = pd.to_numeric(df['pos'])
    chr_df = df[df['contig'].astype(str) == contig_dict['name']]

    if not chr_df.empty:

        mean_depth_columns = []
        colors = []
        if 'illumina_depth' in chr_df.columns:
            mean_depth_columns.append('illumina_depth')
            colors.append('#EECC66')
            # dark : 997700
            #linestyles.append('-')

        if 'nanopore_depth' in chr_df.columns:
            mean_depth_columns.append('nanopore_depth')
            colors.append('#6699CC')
            # dark : 004488
            #linestyles.append('--')

        if 'pacbio_depth' in chr_df.columns:
            mean_depth_columns.append('pacbio_depth')
            colors.append('EE99AA')
            # dark : 994455
            #linestyles.append(':')

        cov_plot = chr_df.plot.line(x='pos', y=mean_depth_columns, color=colors)
        cov_plot.plot()
        plt.title(f"{args.sample} {contig_dict['name']} coverage")
        cov_plot.figure.savefig(f"{args.out}/{args.sample}_{contig_dict['name']}.png")
        plt.close()

    return args.out + '/' + args.sample + '_' + contig_dict['name'] + '.png'
