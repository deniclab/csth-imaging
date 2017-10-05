import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'Merge analysis output \
                                 csv files.')
parser.add_argument('-d', '--csv_dir', required = True, 
                    help = 'directory containing analysis output csvs.')
args = parser.parse_args()
print(args)
csv_dir = args.csv_dir

os.chdir(csv_dir)
csv_list = [f for f in os.listdir() if f.endswith('.csv') and not f.startswith('.')]
raw_csvs = [f for f in csv_list if 'raw' in f]
summary_csvs = [f for f in csv_list if 'summary' in f]
merged_raw_df = pd.DataFrame()
merged_summary_df = pd.DataFrame()
for c in raw_csvs:
    print('current raw csv: ' + c)
    c_df = pd.read_csv(c)
    merged_raw_df = merged_raw_df.append(c_df)
print('merging of raw csvs complete. moving on to summaries.')
for c in summary_csvs:
    print('current summary csv: ' + c)
    c_df = pd.read_csv(c)
    merged_summary_df = merged_summary_df.append(c_df)
merged_raw_df.to_csv(csv_dir + '/merged_raw_analysis_output.csv')
merged_summary_df.to_csv(csv_dir + '/merged_summary_analysis_output.csv')
