import pandas as pd
import os
import gzip
import io
import numpy as np
import time
import sys
from pathlib import Path

start_time = time.time()
print(start_time)

def load_refseq_data(refseq_file):
    try:
        refseq_df = pd.read_csv(refseq_file, sep='\t',
                                converters={'chrom': str, 'exonStarts': lambda x: x.strip(','),
                                            'exonEnds': lambda x: x.strip(',')})
        print("Loaded refseq data with columns:", refseq_df.columns)
        necessary_columns = ['chrom', 'txStart', 'txEnd', 'gene_symbol']
        if not all(col in refseq_df.columns for col in necessary_columns):
            raise ValueError(f"Required columns {necessary_columns} not all present in the DataFrame")
        return refseq_df
    except Exception as e:
        print(f"Error loading reference data: {e}")
        return None

def parse_refseq_for_boundaries_corrected(refseq_file):
    refseq_df = pd.read_csv(refseq_csv_path, sep='\t',
                            converters={'chrom': str, 'exonStarts': lambda x: x.strip(','),
                                        'exonEnds': lambda x: x.strip(',')})
    refseq_df.columns = [col.strip() for col in refseq_df.columns]
    parsed_data = []

    for _, row in refseq_df.iterrows():
        gene_data = {
            'chrom': row['chrom'],
            'gene_symbol': row['gene_symbol'],
            'strand': row['strand'],
            'exonCount': row['exonCount']  
        }
        exon_starts = list(map(int, row['exonStarts'].split(',')))
        exon_ends = list(map(int, row['exonEnds'].split(',')))
        if row['strand'] == '+':
            gene_data['5\'utr_region'] = f"{row['txStart']}-{exon_starts[0]}"
            gene_data['3\'utr_region'] = f"{exon_ends[-1]}-{row['txEnd']}"
            exons = list(zip(exon_starts, exon_ends))
        else:
            gene_data['5\'utr_region'] = f"{exon_ends[-1]}-{row['txEnd']}"
            gene_data['3\'utr_region'] = f"{row['txStart']}-{exon_starts[0]}"
            exons = list(zip(exon_starts, exon_ends))
        for i, (start, end) in enumerate(exons, 1):
            gene_data[f'exon{i}'] = f'{start}-{end}'
        introns = []
        for i in range(len(exons) - 1):
            start = exons[i][1]
            end = exons[i + 1][0]
            introns.append((start, end))

        for i, (start, end) in enumerate(introns, 1):
            gene_data[f'intron{i}'] = f'{start}-{end}'

        parsed_data.append(gene_data)

    parsed_df = pd.DataFrame(parsed_data)
    return parsed_df

def reverse_exon_intron_order(refseq_df):
    corrected_data = []
    for _, row in refseq_df.iterrows():
        gene_data = row.to_dict()  
        if row['strand'] == '-':
            exon_columns = [f'exon{i}' for i in range(1, int(row['exonCount']) + 1)]
            intron_columns = [f'intron{i}' for i in range(1, int(row['exonCount']))]  
            reversed_exons = {f'exon{i}': gene_data[col] for i, col in enumerate(reversed(exon_columns), 1)}
            reversed_introns = {f'intron{i}': gene_data[col] for i, col in enumerate(reversed(intron_columns), 1) if col in gene_data}
            gene_data.update(reversed_exons)
            gene_data.update(reversed_introns)
        corrected_data.append(gene_data)
    return pd.DataFrame(corrected_data)


def read_vcf(path):
    open_func = gzip.open if path.endswith('.gz') else open
    with open_func(path, 'rt') as f:
        lines = [l for l in f if not l.startswith('##')]
    vcf_df = pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
    sample_start_idx = vcf_df.columns.get_loc('FORMAT') + 1
    
    vcf_df['GT'] = vcf_df.apply(lambda row: extract_field(row, 'GT', sample_start_idx), axis=1)
    vcf_df['PR'] = vcf_df.apply(lambda row: extract_field(row, 'PR', sample_start_idx), axis=1)
    vcf_df['SR'] = vcf_df.apply(lambda row: extract_field(row, 'SR', sample_start_idx), axis=1)

    return vcf_df

def extract_field(row, field_name, sample_start_idx):
    format_fields = row['FORMAT'].split(':')
    sample_fields = row.iloc[sample_start_idx].split(':')
    try:
        field_index = format_fields.index(field_name)
        return sample_fields[field_index]
    except ValueError:
        print(f"Field '{field_name}' not found in FORMAT column.")
        return np.nan  
    except IndexError:
        print(f"Sample data does not contain enough fields for '{field_name}'")
        return np.nan 


def prefilter_variants(df):
    def get_alt_supporting_count(value):
        try:
            wt, alt = map(int, value.split(','))
            return alt
        except (ValueError, AttributeError):
            return np.nan  

    if 'PR' in df.columns:
        df['PR_alt'] = df['PR'].apply(get_alt_supporting_count)
    else:
        df['PR_alt'] = np.nan  # Set to NaN if 'PR' is absent

    if 'SR' in df.columns:
        df['SR_alt'] = df['SR'].apply(get_alt_supporting_count)
    else:
        df['SR_alt'] = np.nan  # Set to NaN if 'SR' is absent

    filtered_df = df[(df['PR_alt'].isna() | (df['PR_alt'] >= 20)) &
                     (df['SR_alt'].isna() | (df['SR_alt'] >= 20))]

    return filtered_df

def preprocess_sv_df(sv_df):
    def extract_end(info_str):
        for field in info_str.split(';'):
            if field.startswith("END="):
                return int(field.split('=')[1])
        return None

    sv_df['end'] = sv_df['INFO'].apply(extract_end)
    sv_df['end'] = sv_df['end'].fillna(sv_df['POS']).astype(int)
    sv_df = sv_df.rename(columns={'CHROM': 'chrom', 'POS': 'start'})
    sv_df = sv_df[['chrom', 'start', 'end', 'ALT', 'INFO', 'GT', 'PR', 'SR']]
    
    return sv_df

def find_genes_for_sv(sv_df, refseq_df):
    sv_df['gene_symbols'] = pd.NA
    sv_df['mane_selects'] = pd.NA

    for sv_index, sv in sv_df.iterrows():
        if 'SVTYPE=BND' in sv['INFO']:
            matching_genes_1 = refseq_df[(refseq_df['chrom'] == sv['chrom']) &
                                         (refseq_df['txStart'] <= sv['end']) &
                                         (refseq_df['txEnd'] >= sv['start'])]
            print(matching_genes_1)
            matching_genes_2 = refseq_df[(refseq_df['chrom'] == sv['temp_chr']) &
                                         (refseq_df['txStart'] <= sv['temp_end']) &
                                         (refseq_df['txEnd'] >= sv['temp_start'])]
            print(matching_genes_2)
            combined_genes = set(matching_genes_1['gene_symbol'].unique()).union(
                set(matching_genes_2['gene_symbol'].unique())
            )
            combined_mane = set(matching_genes_1['name'].unique()).union(
                set(matching_genes_2['name'].unique())
            )

            if combined_genes:
                sv_df.at[sv_index, 'gene_symbols'] = ','.join(combined_genes)
                sv_df.at[sv_index, 'mane_selects'] = ','.join(combined_mane)
        else:
            matching_genes = refseq_df[(refseq_df['chrom'] == sv['chrom']) &
                                       (refseq_df['txStart'] <= sv['end']) &
                                       (refseq_df['txEnd'] >= sv['start'])]
            if not matching_genes.empty:
                sv_df.at[sv_index, 'gene_symbols'] = ','.join(matching_genes['gene_symbol'].unique())
                sv_df.at[sv_index, 'mane_selects'] = ','.join(matching_genes['name'].unique())

    return sv_df




def parse_refseq_for_boundaries_corrected(refseq_file):
    refseq_df = pd.read_csv(refseq_file, sep='\t',
                            converters={'chrom': str, 'exonStarts': lambda x: x.strip(','),
                                        'exonEnds': lambda x: x.strip(',')})
    refseq_df.columns = [col.strip() for col in refseq_df.columns]
    parsed_data = []

    for _, row in refseq_df.iterrows():
        gene_data = {
            'chrom': row['chrom'],
            'gene_symbol': row['gene_symbol'],
            'strand': row['strand'],
            'exonCount': row['exonCount']
        }
        exon_starts = list(map(int, row['exonStarts'].split(',')))
        exon_ends = list(map(int, row['exonEnds'].split(',')))

        if row['strand'] == '+':
            gene_data['5\'utr_region'] = f"{row['txStart']}-{row['cdsStart']}"
            gene_data['3\'utr_region'] = f"{row['cdsEnd']}-{row['txEnd']}"
        else:
            gene_data['5\'utr_region'] = f"{row['cdsEnd']}-{row['txEnd']}"
            gene_data['3\'utr_region'] = f"{row['txStart']}-{row['cdsStart']}"

        exons = list(zip(exon_starts, exon_ends))
        for i, (start, end) in enumerate(exons, 1):
            gene_data[f'exon{i}'] = f'{start}-{end}'

        introns = []
        for i in range(len(exons) - 1):
            intron_start = exons[i][1]
            intron_end = exons[i + 1][0]
            introns.append((intron_start, intron_end))

        for i, (start, end) in enumerate(introns, 1):
            gene_data[f'intron{i}'] = f'{start}-{end}'

        parsed_data.append(gene_data)

    parsed_df = pd.DataFrame(parsed_data)
    return parsed_df


def annotate_sv_positions(sv_df, refseq_df):
    print("Columns in sv_df:", sv_df.columns)

    if 'gene_symbols' not in sv_df.columns:
        print("Error: 'gene_symbols' column not found in sv_df")
        return sv_df

    # Initialize new columns
    sv_df['breakend1_pos'] = '.'
    sv_df['breakend2_pos'] = '.'
    sv_df['breakend1_flag'] = '.'
    sv_df['breakend2_flag'] = '.'

    def find_feature_position(sv_row, refseq_df, chrom, pos):
        matches = []
        gene_symbols = sv_row['gene_symbols']

        if pd.notna(gene_symbols):  # Check if gene_symbols is not NaN/NA
            gene_list = gene_symbols.split(',') if ',' in gene_symbols else [gene_symbols]

            for gene_symbol in gene_list:
                gene_symbol = gene_symbol.strip()
                gene_rows = refseq_df[(refseq_df['gene_symbol'] == gene_symbol) & (refseq_df['chrom'] == chrom)]
                if gene_rows.empty:
                    continue

                for _, gene_row in gene_rows.iterrows():
                    features = [col for col in gene_row.index if 'exon' in col or 'intron' in col]
                    for feature in features:
                        feature_range = str(gene_row[feature])
                        if '-' in feature_range:
                            start_feature, end_feature = map(int, feature_range.split('-'))
                            if start_feature <= pos <= end_feature:
                                matches.append(f"{gene_symbol}:{feature}")
        return matches if matches else ['.']

    def flag_utr_position(pos, gene_row, strand, gene_symbol):
        five_utr_range = str(gene_row["5'utr_region"])
        three_utr_range = str(gene_row["3'utr_region"])

        flags = set()  # Use a set to avoid duplicate flags

        if '-' in five_utr_range:
            start_utr, end_utr = map(int, five_utr_range.split('-'))
            if strand == '+':
                if start_utr - 20 <= pos <= start_utr + 20:
                    flags.add(f"{gene_symbol}:5' UTR +- 20")
                if start_utr <= pos <= end_utr:
                    flags.add(f"{gene_symbol}:5' UTR")
            else:
                if end_utr - 20 <= pos <= end_utr + 20:
                    flags.add(f"{gene_symbol}:5' UTR +- 20")
                if start_utr <= pos <= end_utr:
                    flags.add(f"{gene_symbol}:5' UTR")

        if '-' in three_utr_range:
            start_utr, end_utr = map(int, three_utr_range.split('-'))
            if strand == '+':
                if end_utr - 20 <= pos <= end_utr + 20:
                    flags.add(f"{gene_symbol}:3' UTR +- 20")
                if start_utr <= pos <= end_utr:
                    flags.add(f"{gene_symbol}:3' UTR")
            else:
                if start_utr - 20 <= pos <= start_utr + 20:
                    flags.add(f"{gene_symbol}:3' UTR +- 20")
                if start_utr <= pos <= end_utr:
                    flags.add(f"{gene_symbol}:3' UTR")

        return flags  # Return flags as a set

    for sv_index, sv_row in sv_df.iterrows():
        if 'SVTYPE=BND' in sv_row['INFO']:
            # Handling for BND variants
            breakend1_positions = find_feature_position(sv_row, refseq_df, sv_row['chrom'], sv_row['start'])
            breakend2_positions = find_feature_position(sv_row, refseq_df, sv_row['temp_chr'], sv_row['temp_start'])

            sv_df.at[sv_index, 'breakend1_pos'] = ';'.join(breakend1_positions) if breakend1_positions else '.'
            sv_df.at[sv_index, 'breakend2_pos'] = ';'.join(breakend2_positions) if breakend2_positions else '.'

            if pd.notna(sv_row['gene_symbols']):
                gene_list = sv_row['gene_symbols'].split(',') if ',' in sv_row['gene_symbols'] else [sv_row['gene_symbols']]
                breakend1_flags = set()
                breakend2_flags = set()

                for gene_symbol in gene_list:
                    gene_symbol = gene_symbol.strip()
                    gene_rows = refseq_df[(refseq_df['gene_symbol'] == gene_symbol) & (refseq_df['chrom'] == sv_row['chrom'])]
                    if gene_rows.empty:
                        continue

                    for _, gene_row in gene_rows.iterrows():
                        strand = gene_row['strand']
                        breakend1_flags.update(flag_utr_position(sv_row['start'], gene_row, strand, gene_symbol))

                    gene_rows_temp = refseq_df[(refseq_df['gene_symbol'] == gene_symbol) & (refseq_df['chrom'] == sv_row['temp_chr'])]
                    for _, gene_row in gene_rows_temp.iterrows():
                        strand = gene_row['strand']
                        breakend2_flags.update(flag_utr_position(sv_row['temp_start'], gene_row, strand, gene_symbol))

                breakend1_flags.discard('.')
                breakend2_flags.discard('.')
                sv_df.at[sv_index, 'breakend1_flag'] = ';'.join(breakend1_flags) if breakend1_flags else '.'
                sv_df.at[sv_index, 'breakend2_flag'] = ';'.join(breakend2_flags) if breakend2_flags else '.'
        else:
            # Existing logic for non-BND variants
            sv_row['pos'] = sv_row['start']
            breakend1_positions = find_feature_position(sv_row, refseq_df, sv_row['chrom'], sv_row['start'])
            sv_row['pos'] = sv_row['end']
            breakend2_positions = find_feature_position(sv_row, refseq_df, sv_row['chrom'], sv_row['end'])

            sv_df.at[sv_index, 'breakend1_pos'] = ';'.join(breakend1_positions) if breakend1_positions else '.'
            sv_df.at[sv_index, 'breakend2_pos'] = ';'.join(breakend2_positions) if breakend2_positions else '.'

            # Apply UTR flagging for non-BND variants
            if pd.notna(sv_row['gene_symbols']):
                gene_list = sv_row['gene_symbols'].split(',') if ',' in sv_row['gene_symbols'] else [sv_row['gene_symbols']]
                breakend1_flags = set()
                breakend2_flags = set()

                for gene_symbol in gene_list:
                    gene_symbol = gene_symbol.strip()
                    gene_rows = refseq_df[(refseq_df['gene_symbol'] == gene_symbol) & (refseq_df['chrom'] == sv_row['chrom'])]
                    if gene_rows.empty:
                        continue

                    for _, gene_row in gene_rows.iterrows():
                        strand = gene_row['strand']
                        breakend1_flags.update(flag_utr_position(sv_row['start'], gene_row, strand, gene_symbol))
                        breakend2_flags.update(flag_utr_position(sv_row['end'], gene_row, strand, gene_symbol))

                breakend1_flags.discard('.')
                breakend2_flags.discard('.')
                sv_df.at[sv_index, 'breakend1_flag'] = ';'.join(breakend1_flags) if breakend1_flags else '.'
                sv_df.at[sv_index, 'breakend2_flag'] = ';'.join(breakend2_flags) if breakend2_flags else '.'

    return sv_df


def handle_bnd_variants(sv_df):
    sv_df['temp_chr'] = sv_df['chrom']
    sv_df['temp_start'] = sv_df['start']
    sv_df['temp_end'] = sv_df['start']

    for sv_index, sv_row in sv_df.iterrows():
        if 'SVTYPE=BND' in sv_row['INFO']:
            alt = sv_row['ALT']
            if pd.isna(alt):
                print(f"ALT field is missing at index {sv_index}. Skipping this variant.")
                continue
            if ']' in alt or '[' in alt:
                mate_info = alt.split('[' if '[' in alt else ']')
                mate_chr_pos = mate_info[1] if len(mate_info) > 1 else mate_info[0]
                parts = mate_chr_pos.split(':')
                if len(parts) >= 2:
                    mate_chr = parts[0]
                    mate_pos = parts[1]
                else:
                    # Handle cases where the split does not result in expected format
                    print(f"Warning: Unexpected format in mate_chr_pos: {mate_chr_pos}")
                    mate_chr, mate_pos = None, None

                #mate_chr, mate_pos = mate_chr_pos.split(':')
                mate_pos = int(mate_pos)

                sv_df.at[sv_index, 'temp_chr'] = mate_chr.lower()
                sv_df.at[sv_index, 'temp_start'] = mate_pos
                sv_df.at[sv_index, 'temp_end'] = mate_pos

                sv_df.at[sv_index, 'end'] = sv_row['start']

                end_field_1 = f"END={mate_pos}"
                sv_df.at[sv_index, 'INFO'] = sv_row['INFO'] + ';' + end_field_1

    return sv_df

def count_affected_exons(sv_row, refseq_df):
    gene_symbols = sv_row['gene_symbols'].split(',')
    sv_start = sv_row['start']
    sv_end = sv_row['end']
    exons_affected_list = []
    for gene in gene_symbols:
        gene_data = refseq_df[refseq_df['gene_symbol'] == gene]
        num_affected_exons = 0
        if not gene_data.empty:
            for exon_col in [col for col in gene_data.columns if 'exon' in col]:
                exon_ranges = str(gene_data.iloc[0][exon_col]).split(',')
                for exon_range in exon_ranges:
                    if '-' in exon_range:
                        exon_start, exon_end = map(int, exon_range.split('-'))
                        if (exon_start <= sv_start <= exon_end) or (exon_start <= sv_end <= exon_end) or (sv_start <= exon_start and sv_end >= exon_end):
                            num_affected_exons += 1
        exons_affected_list.append(num_affected_exons)
    return exons_affected_list

def extract_svlen_abs(info):
    for part in info.split(';'):
        if part.startswith('SVLEN='):
            try:
                return abs(int(part.split('=')[1]))
            except ValueError:
                return np.nan
    return np.nan

def filter_variants(sv_df):
    bnd_mask = sv_df['INFO'].str.contains('SVTYPE=BND')
    non_bnd_mask = ~bnd_mask & (sv_df['max_exon_affected'] != 0)
    combined_mask = bnd_mask | non_bnd_mask
    filtered_df = sv_df[combined_mask]
    filtered_df = filtered_df.drop(columns=['temp_chr', 'temp_start', 'temp_end'])
    return filtered_df

def merge_vcfs(manta_df, cnvkit_df):
    merged_df = pd.concat([manta_df, cnvkit_df], axis=0, ignore_index=True)
    return merged_df

def process_file(sv_path, output_dir, refseq_df, processed_refseq_df, sample_id, cnvkit_vcf):
    print(f"Starting to process: {sv_path}")

    sv_df = read_vcf(sv_path)
    print(sv_df)
    cnvkit_df = read_vcf(cnvkit_vcf)
    sv_df = pd.concat([sv_df, cnvkit_df], ignore_index = True)
    sv_df.to_csv('merged_test.tsv', index = False)
    print('saved_df')
    print(sv_df)
    print("VCF loaded into DataFrame.")
    print("Columns in sv_df:", sv_df.columns)  
    #cnvkit_df = read_vcf(cnvkit_vcf)
    #cnvkit_deletions_df = filter_pms2_chek2_deletions(cnvkit_df)
    #cnvkit_df.to_csv('cnvkit.csv')
    #cnvkit_deletions_df.to_csv('cnvkit_deletion.csv')
    #print("CNVkit deletions for PMS2 and CHEK2:")
    #print(cnvkit_deletions_df)

    #merged_df = merge_vcfs(sv_df, cnvkit_deletions_df)
    #merged_df.to_csv('merged_df_test.csv')

    #sv_df = prefilter_variants(merged_df)
    ##########sv_df = prefilter_variants(sv_df)
    print("filtered sv")
    sv_df.to_csv('temp.tsv')
    
    sv_df = preprocess_sv_df(sv_df)
    print("Structural variants preprocessed.")
    print(sv_df)
    sv_df.to_csv('preprocessed.csv')
    sv_df = handle_bnd_variants(sv_df)
    print("BND variants processed and END positions annotated.")
    sv_df.to_csv('bnd.tsv')

    sv_df = find_genes_for_sv(sv_df,refseq_df)
    print("Genes associated with SVs found.")
    sv_df.to_csv('found_genes.tsv')

    sv_df = annotate_sv_positions(sv_df, processed_refseq_df)
    sv_df.to_csv('test.tsv')
    print("SV positions annotated.")
    print(sv_df.columns)
    sv_df = sv_df[sv_df['gene_symbols'].notna() & (sv_df['gene_symbols'] != '')]
    print("Filtered entries with gene.")

    sv_df['exons_affected_each_gene'] = sv_df.apply(lambda row: count_affected_exons(row, processed_refseq_df), axis=1)
    print("Counted exons affected by each SV.")

    sv_df['exons_affected_each_gene'] = sv_df['exons_affected_each_gene'].apply(lambda x: ','.join(map(str, x)))
    print("Formatted the count of affected exons.")

    sv_df['max_exon_affected'] = sv_df['exons_affected_each_gene'].apply(lambda x: max([int(i) for i in x.split(',')]))
    print("Determined the maximum number of exons affected.")

    sv_df['SV_length'] = sv_df['INFO'].apply(extract_svlen_abs)
    print("Extracted SV length from INFO.")

    filtered_df = filter_variants(sv_df)

    filtered_df['sample_id'] = sample_id

    #filtered_df = sv_df[sv_df['max_exon_affected'] != 0]
    #filtered_df = filtered_df.drop(columns=['ALT', 'temp_chr', 'temp_start', 'temp_end'])
    print("Filtered SVs where exons are affected.")
    sample_id = os.path.basename(sv_path).split('.')[0]
    file_type = 'str' if 'str' in sv_path else 'sv'

    output_path = os.path.join(output_dir, f'{sample_id}.{file_type}.tsv')
    filtered_df.to_csv(output_path, sep='\t', index=False)
    print(f"Processed data saved to: {output_path}")
    end_time = time.time()
    print(end_time)
    elapse_time = end_time - start_time

    
def main():
    # Configuration
    input_sv_vcf = sys.argv[1]  
    sample_id = sys.argv[2]
    cnvkit_vcf = sys.argv[3]
    output_dir = sys.argv[4]
    refseq_csv_path = "./refseq_original.tsv"

    # Construct the VCF file path using the sample ID and template
    #sv_path = sv_file_template.format(sample_id=sample_id)
    sv_path=input_sv_vcf
    print(f"VCF path for sample {sample_id}: {sv_path}")

    # Load the reference data
    refseq_df = pd.read_csv(refseq_csv_path, sep='\t')
    processed_refseq_df = parse_refseq_for_boundaries_corrected(refseq_csv_path)
    processed_refseq_df = reverse_exon_intron_order(processed_refseq_df)

    # Check if the VCF file exists before processing
    if os.path.exists(sv_path):
        process_file(sv_path, output_dir, refseq_df, processed_refseq_df,sample_id, cnvkit_vcf)
    else:
        print(f"Warning: VCF file for sample {sample_id} not found at {sv_path}")

if __name__ == "__main__":
    main()

