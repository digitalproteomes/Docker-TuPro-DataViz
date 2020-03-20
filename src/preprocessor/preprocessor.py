"Preprocesses Spectronaut output for TP-Viz processing"

import numpy as np
import pandas as pd
import argparse
from collections import Counter


def get_arguments(parser):
  """Set up command line parameters
  """
  parser.add_argument("-i", "--infile",
                      help="""The input tsv file. This table must have a 
                      Gene and a log2fc column. A PG.ProteinAcessions column 
                      must be present and contain Uniprot entries.""",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-o", "--outfile",
                      help="The output tsv filename root.",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-p", "--peptidefile",
                      help="""The input peptide tsv file.""",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-g", "--uniprot2genefile",
                      help="""The mapping of Uniprot to Gene names. 
                      Required columns Entry (Uniprot Entry) and Gene (Gene name).""",
                      metavar="FILE",
                      default='/usr/local/data/uniprotEntry2gene')
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

      
# TODO: Let user know which entries have been removed
def map_genes(intensity_matrix, uniprot2genefile):
  """Map Uniprot entries to genes.

  CAVE AT: Entries from the intensity_matrix that could not be mapped to a gene are removed!

  Returns: 
  The intensity_matrix DataFrame with an additional column ('Gene') with the gene name.
  """
  uniprot2gene = pd.read_csv(uniprot2genefile, sep='\t', usecols=['Entry', 'Gene'])
  uniprot2gene.set_index('Entry', inplace=True)
  mapped_matrix = intensity_matrix.merge(uniprot2gene,
                                     left_index=True, right_index=True)
  # It is still possible that a Uniprot is entry is in the mapping
  # file, with an actual gene name associated with it. Clean this as well.
  mapped_matrix['Gene'].replace('', np.nan, inplace=True)
  mapped_matrix.dropna(subset=['Gene'], inplace=True)
  return mapped_matrix 


def process_refs(intensity_matrix, ref_start, ref_end, outfile_root):
  """Calculates the median across the reference sample intensities

  - Extracts columns ref_start to end_start
  - Calculates median and adds it to a new 'ref_median' column (NaNs are ignored)

  Returns:
  Dataframe with the median protein abundances across all references
  """
  reference_matrix = intensity_matrix.iloc[:, ref_start:ref_end]
  print('Extracting reference samples:\n{}'.format(reference_matrix.columns))
  print(reference_matrix.loc['A0AVT1'])
  reference_matrix['ref_median'] = reference_matrix.apply(lambda x: np.nanmedian(x), axis=1)
  return pd.DataFrame(reference_matrix['ref_median'])


def process_samples(intensity_matrix, sample_start, sample_end, outfile_root):
  """Extracts and preprocesses the sample intensities

  - Extracts columns sample_start to end_start
    
  Returns:
  Dataframe with the proteins intensities across the samples
  """
  sample_matrix = intensity_matrix.iloc[:, sample_start:sample_end]
  print('Extracting samples:\n{}'.format(sample_matrix.columns))
  return sample_matrix


def fix_log2fc_nan(row, max_log2fc, min_log2fc, offset=20):
  """Reassign log2fc nan values.

  Returns:
  If log2fc is not NaN return log2fc
  If reference is NaN assign max + offset
  If reference is not NaN assign min - offset
  """
  if not np.isnan(row['log2fc']):
    return row['log2fc']
  elif np.isnan(row['ref_median']):
    return max_log2fc + offset
  else:
    return min_log2fc - offset
  
# TODO: Currently only one sample is considered.
#       Once the downstream scripts have been adapted to expect log2_sample_name and
#       to expect more than one sample per analysis, rename log2fc column.
def calculate_log2fc(reference_matrix, sample_matrix, fake_nan=True):
  """Calculates log2fc(sample intensity/median_reference intensity)
  """
  log2fc_df = pd.DataFrame()
  for sample_name, sample_intensities in sample_matrix.iteritems():
    log2fc = np.log2(sample_intensities / reference_matrix['ref_median'])
    if fake_nan:
      ref_log2fc = pd.concat([reference_matrix['ref_median'], log2fc], axis=1)
      ref_log2fc.columns = ['ref_median', 'log2fc']
      max_log2fc = ref_log2fc['log2fc'].max()
      min_log2fc = ref_log2fc['log2fc'].min()
      log2fc = ref_log2fc.apply(lambda x: fix_log2fc_nan(x, max_log2fc, min_log2fc), axis=1)
    log2fc_df['log2fc'] = log2fc
  return reference_matrix.merge(log2fc_df, left_index=True, right_index=True)


def get_protein_with_proteotypic_evidence(peptide_infile, intensity_matrix):
  """Extracts a list of proteins with proteotypic evidence
  """
  peptide_matrix = pd.read_csv(peptide_infile, sep='\t')
  unique_prot_count = peptide_matrix.groupby('EG.PrecursorId')['PG.ProteinAccessions'].nunique()
  proteotypic_peptides =  unique_prot_count[unique_prot_count == 1]
  peptide_matrix[peptide_matrix['EG.PrecursorId'].isin(proteotypic_peptides)]
  proteotypic_proteins = set(peptide_matrix['PG.ProteinAccessions'])
  intensity_matrix = intensity_matrix[intensity_matrix['PG.ProteinAccessions'].isin(proteotypic_proteins)]
  return intensity_matrix


def zero_low_intensity(intensity_matrix, intensity_col_start,
                       intensity_col_end, threshold=3):
  """Set intensity values below threshold to NaNs

  This is a very crude compensation to remove imputed values.

  Args:
  - intensity_col_start: the first column with ref/samples intensities
  - intensity_col_end: the last column with ref/sample intensities
  - threshold: the min intensity threshold
  
  Returns:
  intensity_matrix with values below threshold set as NaNs
  """
  for col in range(intensity_col_start, intensity_col_end+1):
    intensity_matrix.iloc[:,col][intensity_matrix.iloc[:, col] < 3] = np.nan
  return intensity_matrix

  
def main():
  """Main
  """
  # Get user parameters
  parser = argparse.ArgumentParser(description="""Preprocess spectronaut output 
table for TP-Viz processing.
  
Creates 2 files:
  - _cleaned.tsv:       same as input protein matrix, with:
    - 'CONT' and 'iRT' entries removed.
    - proteins with no proteotypic evidence removed (NOTE: protein ratio are not 
      recalculated based on proteotypic peptides).
  - _preprocessed.tsv:  starting from the cleaned matrix:
    - sets PG.ProteinAccessions as Entry column
    - adds a ref_median column (with the median of the reference intensities)
    - adds a Gene column (if a protein cannot be mapped it is removed)
    - adds a log2fc column with log2(ref_median/sample) (NOTE: currently, only 
      one sample per matrix is considered.)
  """)
  args = get_arguments(parser)

  infile = args.infile
  outfile_root = args.outfile
  pepfile = args.peptidefile
  uniprot2genefile = args.uniprot2genefile
  ref_start = args.referencecols[0] - 1
  ref_end = args.referencecols[1] 
  sample_start = args.samplecols[0] - 1
  sample_end = args.samplecols[1] 

  # Preprocess intensity matrix
  intensity_matrix = pd.read_csv(infile, sep='\t')
  intensity_matrix = intensity_matrix[~(intensity_matrix['PG.ProteinAccessions'].str.contains('CONT') |
                                        (intensity_matrix['PG.ProteinAccessions'].str.contains('iRT')))]
  intensity_matrix = get_protein_with_proteotypic_evidence(pepfile, intensity_matrix)
  intensity_matrix = zero_low_intensity(intensity_matrix,
                                        min(ref_start, sample_start),
                                        min(ref_end, sample_end))
  intensity_matrix.to_csv(outfile_root + '_cleaned.tsv', sep='\t', index=False, na_rep='NaN')
  intensity_matrix.set_index('PG.ProteinAccessions', inplace=True)
  intensity_matrix = map_genes(intensity_matrix, uniprot2genefile)
  
  # Extract median reference intensities
  # We have to subtract 1 from start and end since we have removed the ProteinAccessions column
  reference_matrix = process_refs(intensity_matrix, ref_start - 1, ref_end - 1, outfile_root)

  # Extract sample intensities
  # We have to subtract 1 from start and end since we have removed the ProteinAccessions column
  sample_matrix = process_samples(intensity_matrix, sample_start - 1, sample_end - 1, outfile_root)

  # Calculate log2FC
  print(reference_matrix.head())
  print(sample_matrix.head())
  log2fc_df = calculate_log2fc(reference_matrix.merge(intensity_matrix['Gene'], left_index=True, right_index=True), sample_matrix)
  log2fc_df.to_csv(outfile_root + '_preprocessed.tsv', sep='\t', index_label='Entry', na_rep='NaN')

if __name__ == "__main__":
  main()
