"Runs some basic QC on the input table for TP-Viz processing."

import pandas as pd
import numpy as np
import argparse


def get_arguments(parser):
  """Set up command line parameters
  """
  parser.add_argument("-i", "--infile",
                      help="""The input tsv file. This table must have a 
                      Gene and a log2fc column. A PG.ProteinAcessions column 
                      must be present and contain Uniprot entries.""",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-r", "--referencecols",
                      help="""The first and last reference columns separated 
                      by space (e.g. -r 4 18)""",
                      nargs=2,
                      type=int,
                      metavar=("start_col", "end_col"),
                      required=True)
  parser.add_argument("-s", "--samplecols",
                      help="""The first and last ample columns separated by 
                      space (e.g. -r 19 20).
                      WARNING!! CURRENTLY ONLY ONE SAMPLE IMPLEMENTED""",
                      nargs=2,
                      type=int,
                      metavar=("start_col", "end_col"),
                      required=True)
  return parser.parse_args()


def qc_input(intensity_matrix, ref_start, ref_end, sample_start, sample_end):
  """Perform some basic QC on the input matrix
  """
  dup_protein_entries = [protein for protein, count
                         in Counter(intensity_matrix['PG.ProteinAccessions']).items()
                         if count > 1]
  if len(dup_protein_entries) > 0:
    raise SystemExit("""Found duplicate entries for the following proteins:
    {}
    Please manually fix input file and re-run.""".format(dup_protein_entries))

  reference_matrix = intensity_matrix.iloc[:, ref_start:ref_end]
  for (columnName, columnData) in reference_matrix.iteritems():
    try:
      columnData.astype(np.float64)
    except:
      raise ValueError("Error reference column {} contain non-numeric value!\nPlease fix input matrix.".format(columnName))

  sample_matrix = intensity_matrix.iloc[:, sample_start:sample_end]
  for (columnName, columnData) in sample_matrix.iteritems():
    try:
      columnData.astype(np.float64)
    except:
      raise ValueError("Error sample column {} contain non-numeric value!\nPlease fix input matrix.".format(columnName))


def main():
  """Main
  """
  # Get user parameters
  parser = argparse.ArgumentParser(description="""Runs some basic QC on the input
  table for TP-Viz processing.""")
  args = get_arguments(parser)

  infile = args.infile
  ref_start = args.referencecols[0] - 1
  ref_end = args.referencecols[1] 
  sample_start = args.samplecols[0] - 1
  sample_end = args.samplecols[1] 

  # Preprocess intensity matrix
  intensity_matrix = pd.read_csv(infile, sep='\t')
  qc_input(intensity_matrix, ref_start, ref_end, sample_start, sample_end)
  
if __name__ == "__main__":
  main()
