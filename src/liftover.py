# File: liftover.py
# Author: Tomas Bencomo
# Description: Liftover genomic coordinates to hg38 for the Durinck/Cho
# and Pickering mutation callsets. The Durinck/Cho callset is aligned
# to hg18 and the genomic positions are 1 indexed (confirmed by checking some reference alleles against UCSC Genome Browser).
# Pickering is aligned to hg19 and is also 1 indexed based.
# WARNING: It is important to note pyliftover treats coordinates as 0 indexed,
# so it is critical to subtract 1 from coordinates before using pyliftover. 
#After coordinates have been lifted over, they'll need to be converted back to 1 based

import os
import numpy as np
from pyliftover import LiftOver
import pandas as pd

def concat_cho_data():
    filepath = os.path.join('../data/original-callsets', 'stab_1_8WES.xls')
    sheets = pd.read_excel(filepath, sheet_name = None)
    mutation_tables = []
    for sheet in sheets:
        sheets[sheet]['patient'] = sheet
        mutation_tables.append(sheets[sheet])
    cho = pd.concat(mutation_tables, ignore_index=True)
    return cho

def liftover_cho(df):
    lo = LiftOver('hg18', 'hg38')
    def lift_coord(row):
        chrom = 'chr' + str(row['Chromosome'])
        pos = row['Genomic position'] - 1
        result = lo.convert_coordinate(chrom, pos)
        if len(result) == 0:
            print(f"Didn't find hg38 coordinate for {row['Chromosome']}:{row['Genomic position']}")
            return 'NA'
        return result[0][1] + 1
    df['Genomic position'] = df.apply(lift_coord, axis=1)
    return df

def load_pickering():
    filepath = os.path.join('../data/original-callsets', 'NIHMS635298-supplement-3_39WES.xlsx')
    df = pd.read_excel(filepath, header = 2)
    return df

def liftover_pickering(df):
    lo = LiftOver('hg19', 'hg38')
    chroms = ('chr' + df['Chromosome'].astype(str)).tolist()
    startpos = df['Start_position'].tolist()
    endpos = df['End_position'].tolist()
    rows = zip(chroms, startpos, endpos)
    new_startpos = []
    new_endpos = []
    for row in rows:
        new_start = lo.convert_coordinate(row[0], row[1]-1)
        if len(new_start) == 0:
            print(f"Didn't find hg38 coordinate for {row[0]}:{row[1]}-{row[2]}")
            new_startpos.append('NA')
        else:
            new_startpos.append(new_start[0][1]+1)
        new_end = lo.convert_coordinate(row[0], row[2]-1)
        if len(new_end) == 0:
            print(f"Didn't find hg38 coordinate for {row[0]}:{row[1]}-{row[2]}")
            new_endpos.append('NA')
        else:
            new_endpos.append(new_end[0][1] + 1)
    df['Start_position'] = new_startpos
    df['End_position'] = new_endpos
    return df
    
def main():
    cho = concat_cho_data()
    print("Lifting over Durinck coordinates")
    cho = liftover_cho(cho)

    pickering = load_pickering()
    print("Lifting over Pickering coordinates")
    pickering = liftover_pickering(pickering)

    print("Saving updated files...")
    cho.to_csv('../data/durinck_mutations_hg38.csv', index=False, header=True)
    pickering.to_csv('../data/pickering_mutations_hg38.csv', index=False, header=True) 

if __name__ == '__main__':
    main()


