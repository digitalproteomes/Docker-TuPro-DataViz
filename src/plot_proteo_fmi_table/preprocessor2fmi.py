#!/usr/bin/env python

import pandas as pd
import numpy as np
import os
import re
import argparse

def load_fmi_targets(fmi_file):
  """Loads all the FMI gene names into a list
  """
  fmi_targets = []
  with open(fmi_file, 'r') as fin:
    for gene in fin:
      fmi_targets.append(gene.strip())
  return fmi_targets


def is_fmi_target(target, fmi_targets):
  """Checks if target is an fmi_target
  """
  if target in fmi_targets:
    return True
  else:
    return False
      

def assemble_fmi_df(infile, fmifile):
  """Assembles a dataframe for FMI targets only
  """
  # Load FMI targets
  fmi_targets = load_fmi_targets(fmifile)
    
  # Load input data
  input_df = pd.read_csv(infile, sep='\t')

  # Only keep the first name for a gene, to mitigate non-uniqueness
  input_df['Gene'] = input_df.apply(lambda x: x['Gene'].split('|')[0], axis=1)

  # Extract FMI targets
  fmi_df = input_df[input_df.apply(lambda x: is_fmi_target(x['Gene'], fmi_targets), axis=1)]
  fmi_df = fmi_df[['Gene', 'log2fc']]

  # Add all targets that were not seen with log2fc = NaN
  fmi_na_df = pd.DataFrame(list(set(fmi_targets).difference(set(fmi_df['Gene']))))
  fmi_na_df.columns = ['Gene']
  fmi_na_df['log2fc'] = np.nan

  fmi_df = fmi_df.append(fmi_na_df)
  fmi_df.sort_values(['Gene'], inplace=True)

  return fmi_df
  

  
def get_arguments(parser):
  """Set up command line parameters
  """
  parser.add_argument("-i", "--infile",
                      help="The input csv file. This table must have a Gene and a log2fc column.",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-f", "--fmifile",
                      help="A file listing FMI target gene one per line.",
                      metavar="FILE",
                      default='/usr/local/data/fmi_gene_list.txt')
  return parser.parse_args()


def main():
  """Main
  """
  parser = argparse.ArgumentParser(description=
                                   """Extracts FMI targets from preprocessor output table""")

  args = get_arguments(parser)
  if not os.path.isfile(args.infile):
    print("ERROR: Cannot read input file")
    exit(1)

  if not os.path.isfile(args.fmifile):
    print("ERROR: Cannot read FMI targets file")
    exit(1)

  fmi_df = assemble_fmi_df(args.infile, args.fmifile)
  out_filename = re.match(r'(.*)\.tsv', args.infile).group(1) + '_fmi.tsv'
  fmi_df.to_csv(out_filename, sep='\t', index=False)


if __name__ == "__main__":
  main()
