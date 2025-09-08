import pandas as pd
import matplotlib

def read_gene_list(genes_file):
    with open(genes_file, 'r') as f:
        return set(line.strip() for line in f if line.strip())

def parse_svtype(info_str, alt_str):
    if isinstance(info_str, str):
        info_upper = info_str.upper()
        if "SVTYPE=DEL" in info_upper:
            return "DEL"
        elif "SVTYPE=DUP" in info_upper:
            return "DUP"

    alt_upper = str(alt_str).upper()
    if "<DEL>" in alt_upper:
        return "DEL"
    elif "<DUP" in alt_upper:  
        return "DUP"

    return "UNK"  

def is_within_gene_position(breakend_pos: str) -> bool:
    """
    Returns True if the breakend annotation suggests it is within a gene
    (in an intron, exon, or UTR). Otherwise False.
    Example breakend_pos: 'SNRNP25:intron1', 'KANK1:5' UTR', '.', etc.
    """
    if not isinstance(breakend_pos, str):
        return False
    if breakend_pos == '.':
        return False

    lower = breakend_pos.lower()
    return ('intron' in lower) or ('exon' in lower) or ('utr' in lower)
def is_exon_position(breakend_pos: str) -> bool:
    if not isinstance(breakend_pos, str):
        return False
    if breakend_pos == '.':
        return False
    
    return 'exon' in breakend_pos.lower()

def is_coding_exon(pos: str, flag: str) -> bool:
    if not isinstance(pos, str) or not isinstance(flag, str):
        return False
    
    pos_lower = pos.lower()
    flag_lower = flag.lower()

    if "exon" not in pos_lower:
        return False

    if "utr" in flag_lower:
        return False

    return True

def calculate_cnv_impact(preloaded_gene_df, cnv_start: int, cnv_end: int, gene_name: str):
    try:
        gene_row = preloaded_gene_df[preloaded_gene_df['gene_symbol'] == gene_name].iloc[0]
    except IndexError:
        return {}

    strand = gene_row['strand']
    exon_count = int(gene_row['exonCount'])
    utr_5_region = gene_row["5'utr_region"]
    utr_3_region = gene_row["3'utr_region"]

    try:
        utr_5_start, utr_5_end = map(int, utr_5_region.split('-'))
        utr_3_start, utr_3_end = map(int, utr_3_region.split('-'))
    except (AttributeError, ValueError):
        utr_5_start, utr_5_end = 0, 0
        utr_3_start, utr_3_end = 0, 0
    exon_columns = [f"exon{i+1}" for i in range(exon_count)]
    exons = gene_row[exon_columns].values.flatten()

    total_protein_coding_bases = 0
    affected_exon_bases = 0

    for exon in exons:
        if pd.isna(exon):
            continue

        exon_start, exon_end = map(int, exon.split('-'))
        if strand == '+':
            if utr_5_start < exon_end and utr_5_end >= exon_start:
                if exon_start >= utr_5_start and exon_end <= utr_5_end:
                    continue
                elif exon_start < utr_5_start < exon_end:
                    exon_start = utr_5_end
                elif exon_start <= utr_5_end < exon_end:
                    exon_start = max(exon_start, utr_5_end)

            if utr_3_start < exon_end and utr_3_end >= exon_start:
                if exon_start >= utr_3_start and exon_end <= utr_3_end:
                    continue
                elif exon_start < utr_3_start < exon_end:
                    exon_end = utr_3_start
                elif exon_start <= utr_3_end < exon_end:
                    exon_end = min(exon_end, utr_3_start)
        else:  # strand == '-'
            if utr_5_start < exon_end and utr_5_end >= exon_start:
                if exon_start >= utr_5_start and exon_end <= utr_5_end:
                    continue
                elif exon_start < utr_5_start < exon_end:
                    exon_end = utr_5_start
                elif exon_start <= utr_5_end < exon_end:
                    exon_end = min(exon_end, utr_5_start)

            if utr_3_start < exon_end and utr_3_end >= exon_start:
                if exon_start >= utr_3_start and exon_end <= utr_3_end:
                    continue
                elif exon_start < utr_3_start < exon_end:
                    exon_start = utr_3_end
                elif exon_start <= utr_3_end < exon_end:
                    exon_start = max(exon_start, utr_3_end)

        coding_bases = max(0, exon_end - exon_start)
        total_protein_coding_bases += coding_bases

        overlap_start = max(cnv_start, exon_start + 1)
        overlap_end = min(cnv_end, exon_end)

        if overlap_start <= overlap_end:
            affected_exon_bases += (overlap_end - overlap_start + 1)

    frame = "inframe" if affected_exon_bases % 3 == 0 else "frameshift"

    return {
        "frame": frame,
        "affected_exon_bases": affected_exon_bases,
        "total_protein_coding_bases": total_protein_coding_bases
    }

def process_dataframe(df, preloaded_gene_df, gene_list):
  
    df['frame'] = None
    df['percent_removed'] = None
    df['type'] = None
    df['remark'] = None
    
    for index, row in df.iterrows():
        cnv_start = row['start']
        cnv_end = row['end']
        
        sv_type_detected = parse_svtype(row.get('INFO', ''), row.get('ALT', ''))
        is_del = (sv_type_detected == "DEL")
        is_dup = (sv_type_detected == "DUP")
        if pd.isna(row['gene_symbols']):
            continue
        genes_in_this_row = row['gene_symbols'].split(',')
        exons_affected_str = row.get('exons_affected_each_gene', '')
        exons_affected_list = [int(x.strip()) for x in exons_affected_str.split(',')] if exons_affected_str else []
        
        panel_genes = set(genes_in_this_row).intersection(gene_list)
        frame_results, percent_removed_results, type_results, remark_results = [], [], [], []
        for gene_name in panel_genes:
            try:
                gene_index = genes_in_this_row.index(gene_name)
            except ValueError:
                continue
            
            exons_affected_for_this_gene = exons_affected_list[gene_index] if 0 <= gene_index < len(exons_affected_list) else 0
            
            gene_data = preloaded_gene_df[preloaded_gene_df['gene_symbol'] == gene_name]
            if gene_data.empty:
                continue
            strand = gene_data.iloc[0]['strand']
            coding_exons = gene_data.iloc[0].get('coding_exons', '')
            if not coding_exons:
                coding_exons = []
            else:
                try:
                    coding_exons = [int(x.strip()) for x in coding_exons.split(',')]
                except ValueError:
                    coding_exons = []
            
            try:
                f_utr_start, f_utr_end = map(int, gene_data.iloc[0]["5'utr_region"].split('-'))
                t_utr_start, t_utr_end = map(int, gene_data.iloc[0]["3'utr_region"].split('-'))
            except (ValueError, AttributeError, IndexError):
                continue
            gene_numeric_start = min(f_utr_start, f_utr_end, t_utr_start, t_utr_end)
            gene_numeric_end = max(f_utr_start, f_utr_end, t_utr_start, t_utr_end)
            cnv_left, cnv_right = cnv_start, cnv_end
            
            if (cnv_left <= gene_numeric_start) and (cnv_right >= gene_numeric_end):
                sv_region = "Whole gene affected"
            else:
                if strand == '+':
                    if cnv_left <= gene_numeric_start and cnv_right < gene_numeric_end:
                        sv_region = "5' UTR to body"
                    elif cnv_left > gene_numeric_start and cnv_right >= gene_numeric_end:
                        sv_region = "3' UTR to body"
                    else:
                        sv_region = "Within gene"
                else:
                    if cnv_left <= gene_numeric_start and cnv_right < gene_numeric_end:
                        sv_region = "3' UTR to body"
                    elif cnv_left > gene_numeric_start and cnv_right >= gene_numeric_end:
                        sv_region = "5' UTR to body"
                    else:
                        sv_region = "Within gene"
            
            result = calculate_cnv_impact(preloaded_gene_df, cnv_left, cnv_right, gene_name)
            if not result:
                continue
            total_coding_bases = result['total_protein_coding_bases']
            affected_bases = result['affected_exon_bases']
            frame = result['frame']
            pct_removed = (affected_bases / total_coding_bases) * 100 if total_coding_bases > 0 else 0.0
            
            coding_exons_affected = count_coding_exons_outside_utr(gene_data, cnv_left, cnv_right, gene_name)
            #coding_exons_affected = [exon for exon in coding_exons if cnv_left <= exon <= cnv_right]
            
            frame_results.append(f"{gene_name}: {frame}")
            percent_removed_results.append(f"{gene_name}: {pct_removed:.2f}%")
            type_results.append(f"{gene_name}: {sv_region}")

            if is_del and sv_region == "3' UTR to body" and pct_removed > 0 and coding_exons_affected >= 2:
                remark_results.append("pathogenic; 3' UTR to body with multiple coding exons affected")
                
            if (is_dup or is_del) and (sv_region == "Within gene") and (frame == "frameshift"):
                remark_results.append("pathogenic; frameshift")

            elif is_del and (sv_region == "Whole gene affected"):
                remark_results.append("pathogenic; whole gene deletion")

            elif is_del and (sv_region == "5' UTR to body") and (pct_removed > 0):
                remark_results.append("pathogenic; 5' UTR to body") 
            
        df.at[index, 'frame'] = "; ".join(frame_results)
        df.at[index, 'percent_removed'] = "; ".join(percent_removed_results)
        df.at[index, 'type'] = "; ".join(type_results)
        df.at[index, 'remark'] = "; ".join(set(remark_results))
    
    return df

def count_coding_exons_outside_utr(gene_data, cnv_start, cnv_end, gene_name):
    try:
        strand = gene_data.iloc[0]['strand']
        exon_count = int(gene_data.iloc[0]['exonCount'])
        utr_5_region = gene_data.iloc[0]["5'utr_region"]
        utr_3_region = gene_data.iloc[0]["3'utr_region"]

        # Parse UTR coordinates
        utr_5_start, utr_5_end = map(int, utr_5_region.split('-'))
        utr_3_start, utr_3_end = map(int, utr_3_region.split('-'))
    except (AttributeError, ValueError, IndexError):
        return 0

    exon_columns = [f"exon{i+1}" for i in range(exon_count)]
    exons = gene_data[exon_columns].values.flatten()

    coding_exons_affected = 0
    for exon in exons:
        if pd.isna(exon):
            continue

        exon_start, exon_end = map(int, exon.split('-'))

        if strand == '+':
            if utr_5_start < exon_end and utr_5_end >= exon_start:
                if exon_start >= utr_5_start and exon_end <= utr_5_end:
                    continue  # Entire exon in 5'UTR
                elif exon_start < utr_5_start < exon_end:
                    exon_start = utr_5_end
                elif exon_start <= utr_5_end < exon_end:
                    exon_start = max(exon_start, utr_5_end)

            if utr_3_start < exon_end and utr_3_end >= exon_start:
                if exon_start >= utr_3_start and exon_end <= utr_3_end:
                    continue  # Entire exon in 3'UTR
                elif exon_start < utr_3_start < exon_end:
                    exon_end = utr_3_start
                elif exon_start <= utr_3_end < exon_end:
                    exon_end = min(exon_end, utr_3_end)
        elif strand == '-':
            # Trim 5'UTR (on - strand, 3'UTR is at the start)
            if utr_5_start < exon_end and utr_5_end >= exon_start:
                if exon_start >= utr_5_start and exon_end <= utr_5_end:
                    continue
                elif exon_start < utr_5_start < exon_end:
                    exon_end = utr_5_start
                elif exon_start <= utr_5_end < exon_end:
                    exon_end = min(exon_end, utr_5_start)

            # Trim 3'UTR (on - strand, 5'UTR is at the end)
            if utr_3_start < exon_end and utr_3_end >= exon_start:
                if exon_start >= utr_3_start and exon_end <= utr_3_end:
                    continue
                elif exon_start < utr_3_start < exon_end:
                    exon_start = utr_3_end
                elif exon_start <= utr_3_end < exon_end:
                    exon_start = max(exon_start, utr_3_end)

        # Check if the adjusted exon overlaps the CNV
        overlap_start = max(cnv_start, exon_start + 1)
        overlap_end = min(cnv_end, exon_end)

        if overlap_start <= overlap_end:
            coding_exons_affected += 1

    return coding_exons_affected



if __name__ == "__main__":
    gene_definition_file = "parsed_refseq.csv"
    input_dataframe_file= "merged_SF_dominant.tsv"
    genes_txt_file = "SF_dominant_LoF_genes.txt"
    output_file = "../summary/SF_dominant_full.tsv"
    # Read gene definitions
    preloaded_gene_df = pd.read_csv(gene_definition_file)
    input_dataframe = pd.read_csv(input_dataframe_file, sep='\t', encoding='ISO-8859-1')

    gene_list = read_gene_list(genes_txt_file)

    output_dataframe = process_dataframe(input_dataframe, preloaded_gene_df, gene_list)
    invalid_mask = ~output_dataframe["GT"].str.match(r'^\d[|/]\d$')
    output_dataframe.loc[invalid_mask, "GT"] = "1/1"
    output_dataframe.to_csv(output_file,sep='\t', index=False)
    print(output_dataframe)


