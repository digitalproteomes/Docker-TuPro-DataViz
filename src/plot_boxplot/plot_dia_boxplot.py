from os.path import basename
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.switch_backend('agg')
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
  parser.add_argument("-o", "--outfile",
                      help="The output tsv filename root.",
                      metavar="FILE",
                      required=True)
  parser.add_argument("-m", "--markersfile",
                      help="""The list of clinical markers.
                      Required columns Uniprot and Gene""",
                      metavar="FILE",
                      default='/usr/local/data/Melanoma_marker_Anja_IDs.csv')
  parser.add_argument("-r", "--referencecols",
                      help="""The first and last reference columns separated 
                      by space (e.g. -r 4 18)""",
                      nargs=2,
                      type=int,
                      metavar=("start_col", "end_col"),
                      required=True)
  parser.add_argument("-s", "--samplecols",
                      help="""The first and last ample columns separated by 
                      space (e.g. -r 19 20).""",
                      nargs=2,
                      type=int,
                      metavar=("start_col", "end_col"),
                      required=True)
  return parser.parse_args()


def process_refs(intensity_matrix, ref_start, ref_end, outfile_root):
  """Calculates the median across the reference sample intensities

  - Extracts columns ref_start to end_start
  - Calculates log2 of intensities

  Returns:
  Dataframe with the log2 reference intensities
  """
  reference_matrix = intensity_matrix.iloc[:, ref_start:ref_end]
  print('Extracting reference samples:\n{}'.format(reference_matrix.columns))
  reference_matrix = np.log2(reference_matrix)
  return pd.DataFrame(reference_matrix)


def process_samples(intensity_matrix, sample_start, sample_end, outfile_root):
  """Extracts and preprocesses the sample intensities

  - Extracts columns sample_start to end_start
  - Calculates log2 of intensities
    
  Returns:
  Dataframe with log2 proteins intensities across the samples
  """
  sample_matrix = intensity_matrix.iloc[:, sample_start:sample_end]
  print('Extracting samples:\n{}'.format(sample_matrix.columns))
  sample_matrix = np.log2(sample_matrix)
  sample_matrix.fillna(0, inplace=True)
  print(sample_matrix.head())
  return sample_matrix


def plot_boxplot(reference_matrix, sample_matrix, outfile_root):
  """Plot the boxplot
  """
  sns.set_style("whitegrid")
  plt.subplots(figsize=(15,10))

  data = reference_matrix.transpose()
  print(data)
  # Boxplots
  g = sns.boxplot(data=data, showmeans=False, color="skyblue")
  # Convert to black and white
  for i, box in enumerate(g.artists):
    box.set_edgecolor('black')
    box.set_facecolor('white')

    # iterate over whiskers and median lines
    for j in range(6*i,6*(i+1)):
      g.lines[j].set_color('black')

  # Add scatterpoints for the sample values
  g = sns.scatterplot(data=sample_matrix, color='red', s=100)

  _ = g.set_xticklabels(reference_matrix.index, rotation=90)
  _ = g.set(xlabel='gene name', ylabel='log2(intensity)', title='Marker intensities in sample vs. reference population')
  _ = g.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0.)

  g.get_figure().savefig(outfile_root + "_DIA-BoxPlot.pdf", bbox_inches='tight')
  return

  
def main():
  """Main
  """
  # Get user parameters
  parser = argparse.ArgumentParser(description="""Creates DIA boxplots
for TP-Viz processing.

This script supports an arbitrary number of samples.
""")
  args = get_arguments(parser)

  infile = args.infile
  outfile_root = args.outfile
  markers_file = args.markersfile
  ref_start = args.referencecols[0] - 1
  ref_end = args.referencecols[1]
  sample_start = args.samplecols[0] -1 
  sample_end = args.samplecols[1]

  # Preprocess intensity matrix
  intensity_matrix = pd.read_csv(infile, sep='\t')
  intensity_matrix.set_index('PG.ProteinAccessions', inplace=True)
  intensity_matrix.index.names = ['Uniprot']

  # Load markers df
  markers_df = pd.read_csv(markers_file, sep=',', index_col='Uniprot')

  # Extract intensity of markers
  markers_intensity = intensity_matrix.merge(markers_df, left_index=True, right_index=True)
  markers_intensity.set_index('Gene', inplace=True)
  markers_intensity.sort_index(inplace=True)
  markers_intensity.to_csv(outfile_root + '_DIA-Marker.tsv', sep='\t')
  
  # Extract markers intensities from reference samples
  # We have to subtract 1 from start and end since we have removed the ProteinAccessions column
  reference_matrix = process_refs(markers_intensity, ref_start - 1, ref_end - 1, outfile_root)

  # Extract markers intensities from actual samples
  # We have to subtract 1 from start and end since we have removed the ProteinAccessions column
  sample_matrix = process_samples(markers_intensity, sample_start - 1, sample_end - 1, outfile_root)

  # Do the actual plot
  plot_boxplot(reference_matrix, sample_matrix, outfile_root)
  

if __name__ == "__main__":
  main()
