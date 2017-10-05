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
merged_df = pd.DataFrame()
for c in csv_list:
    print('current csv: ' + c)
    c_df = pd.read_csv(c)
    merged_df = merged_df.append(c_df)
merged_df = merged_df.to_csv(csv_dir + '/merged_analysis_output.csv')
