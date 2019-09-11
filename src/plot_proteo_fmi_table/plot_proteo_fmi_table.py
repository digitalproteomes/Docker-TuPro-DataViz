#!/usr/bin/env python

import sys
import os
import math
import argparse
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.patches import Patch
import pandas as pd



def assign_color(value):
  """ Return a color depending on value
  - orange if value > 1
  - aqua if value < -1
  - white if 1 < value < 1
  - gray if NaN"""

  if value < -1:
    color = 'aqua'
  elif value > 1:
    color = 'orange'
  elif -1 <= value <= 1:
    color = 'white'
  else:
    color = 'grey'
  return color


def plot_table(protq, col_count=10, out_file_name='fmi_table.png'):
  """Plots FMI targets table with proteomics colors.

  Params:
  - data: a pandas DataFrame with columns ['Gene', 'log2fc']
  - col_count: the number of columns in the resulting table
  - out_file_name: the name of output table figure.
                   NOTE by changing the extension the output format is
                   also automatically adjusted.
  """
  cell_text = []
  cell_color = []
  row_count = math.ceil(protq.shape[0]/col_count)
  pos = 0

  for _ in range(0, row_count):
    row_text = []
    row_color = []
    for _ in range(0, col_count):
      if pos < protq.shape[0]:
        row_text.append(protq['Gene'][pos])
        row_color.append(assign_color(protq['log2fc'][pos]))
      else:
        row_text.append('')
        row_color.append('white')
      pos += 1
    cell_text.append(row_text)
    cell_color.append(row_color)

  _ = plt.figure(figsize=(15, 1))
  plt.table(cellText=cell_text,
            cellLoc='center',
            cellColours=cell_color,
            loc='top')
  _ = plt.axis('tight')
  _ = plt.axis('off')
  legend_elements = [Patch(facecolor='orange', edgecolor='r',
                           label='More than 2-fold up-regulated in sample'),
                     Patch(facecolor='aqua', edgecolor='r',
                           label='More than 2-fold down-regulated in sample'),
                     Patch(facecolor='white', edgecolor='r',
                           label='Less than 2-fold change'),
                     Patch(facecolor='grey', edgecolor='r',
                           label='Not detected')]
  plt.legend(handles=legend_elements)
  plt.savefig(out_file_name, bbox_inches='tight')


def get_arguments(parser):
  """Set up command line parameters
  """
  parser.add_argument("-i", "--infile",
                      help="The input csv file. This table must have a Gene and a log2fc column.",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-o", "--outfile",
                      help="""The name of the figure to generate. 
                      NOTE the extension controls the format.""",
                      metavar="FILE",
                      default='fmi_table.png')
  parser.add_argument("-c", "--colnum",
                      help="The number of columns in the resulting table.",
                      type=int,
                      metavar="N",
                      default=10)

  return parser.parse_args()


def main():
  """Main
  """
  parser = argparse.ArgumentParser(description=
                                   """Generates a table with proteomics 
                                   values corresponding to the FMI targets""")

  args = get_arguments(parser)
  if not os.path.isfile(args.infile):
    print("ERROR: Cannot read input file")
    exit(1)

  protq = pd.read_csv('{}'.format(args.infile), sep='\t')
  plot_table(protq, args.colnum, '{}'.format(args.outfile))


if __name__ == "__main__":
  main()
