import sys
import csv
import gzip
from collections import defaultdict

def parse_manta_tsv(manta_file):
    manta_variants = []
    with open(manta_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader) 
        for row in reader:
            chrom = row[0]  
            start = int(float(row[1]))
            end = int(float(row[2]))
            manta_variants.append((chrom, start, end, row))
    return header, manta_variants

def parse_cnvkit_vcf(cnvkit_file):
    cnvkit_variants = defaultdict(list)
    with gzip.open(cnvkit_file, 'rt') as f:  
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos = fields[0], int(float(fields[1]))

            end = None
            for info_field in fields[7].split(';'):
                if info_field.startswith('END='):
                    end = int(float(info_field.split('=')[1]))
                    break

            if end is not None:
                cnvkit_variants[chrom].append((pos, end))
            else:
                print(f"Warning: No END field found in INFO column for position {pos} on chromosome {chrom}")

    return cnvkit_variants


def calculate_overlap(start1, end1, start2, end2):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start <= overlap_end:
        return overlap_end - overlap_start + 1
    return 0

def find_overlaps(manta_variants, cnvkit_variants):
    overlapping_variants = []
    for chrom, start, end, row in manta_variants:
        manta_length = end - start + 1
        max_overlap_fraction = 0
        for cnvkit_start, cnvkit_end in cnvkit_variants.get(chrom, []):
            overlap_length = calculate_overlap(start, end, cnvkit_start, cnvkit_end)
            if overlap_length > 0:
                overlap_fraction = overlap_length / manta_length
                max_overlap_fraction = max(max_overlap_fraction, overlap_fraction)
        
        row.append(f"{max_overlap_fraction:.4f}")
        overlapping_variants.append(row)
    
    return overlapping_variants

def write_output(header, variants, output_file, size_threshold=25000, overlap_threshold=0.7):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        header.append("Manta_Cnvkit_Overlap_Fraction")
        writer.writerow(header)
        for row in variants:
            variant_size = int(float(row[2])) - int(float(row[1])) + 1
            overlap_fraction = float(row[-1])
            if variant_size <= size_threshold or overlap_fraction >= overlap_threshold:
                writer.writerow(row)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py manta.tsv cnvkit.vcf output.tsv")
        sys.exit(1)

    manta_file = sys.argv[1]
    cnvkit_file = sys.argv[2]
    output_file = sys.argv[3]

    header, manta_variants = parse_manta_tsv(manta_file)
    cnvkit_variants = parse_cnvkit_vcf(cnvkit_file)
    overlapping_variants = find_overlaps(manta_variants, cnvkit_variants)
    write_output(header, overlapping_variants, output_file)

print("Overlap analysis complete. Results written to", output_file)
