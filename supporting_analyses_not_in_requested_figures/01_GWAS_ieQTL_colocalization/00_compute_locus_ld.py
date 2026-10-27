#!/usr/bin/env python3
"""Export lead-SNP LD separately from plotting (adapted from plot_gwas_ieqtl_ld_batch.py).

Requires PLINK 1.9 BED/BIM/FAM, a lead variant present in BIM, and a region
on the same genome assembly. Positions are 1-based inclusive. No allele or
assembly harmonization is performed. Missing LD stays missing, never r2=0.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plink', required=True, help='PLINK 1.9 executable')
    parser.add_argument('--bfile', type=Path, required=True, help='BED/BIM/FAM prefix')
    parser.add_argument('--lead', required=True, help='Exact variant ID in BIM')
    parser.add_argument('--chrom', required=True)
    parser.add_argument('--start', type=int, required=True)
    parser.add_argument('--end', type=int, required=True)
    parser.add_argument('--chr-set', type=int, default=26)
    parser.add_argument('--out-dir', type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.start <= args.end:
        parser.error('Require 1 <= start <= end (1-based inclusive).')
    for suffix in ('.bed', '.bim', '.fam'):
        if not Path(str(args.bfile) + suffix).is_file():
            parser.error('Missing genotype file: ' + str(args.bfile) + suffix)
    # Reject a missing/misplaced lead before starting the external calculation.
    with Path(str(args.bfile) + '.bim').open() as handle:
        leads = [line.split() for line in handle if len(line.split()) >= 4 and line.split()[1] == args.lead]
    if len(leads) != 1 or leads[0][0] != args.chrom or not args.start <= int(leads[0][3]) <= args.end:
        parser.error('Lead must occur exactly once in BIM on the requested chromosome/region.')
    args.out_dir.mkdir(parents=True, exist_ok=False)
    prefix = args.out_dir / 'lead_ld'
    # Set window size to cover the requested region, rather than a fixed 1 Mb.
    window_kb = (args.end - args.start + 999) // 1000 + 1
    command = [args.plink, '--bfile', str(args.bfile), '--chr-set', str(args.chr_set),
               '--chr', args.chrom, '--from-bp', str(args.start), '--to-bp', str(args.end),
               '--ld-snp', args.lead, '--r2', 'gz', '--ld-window', '99999999',
               '--ld-window-kb', str(window_kb), '--ld-window-r2', '0', '--out', str(prefix)]
    version = subprocess.run([args.plink, '--version'], check=True, capture_output=True, text=True).stdout.strip()
    (args.out_dir / 'run.json').write_text(json.dumps({'command': command, 'plink_version': version,
        'coordinate_convention': '1-based inclusive', 'lead': args.lead}, indent=2) + '\n')
    with (args.out_dir / 'process.log').open('w') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    with gzip.open(str(prefix) + '.ld.gz', 'rt') as handle:
        header = next(handle).split()
        if not {'SNP_B', 'R2'}.issubset(header):
            raise ValueError('Expected PLINK 1.9 LD columns SNP_B and R2')
        with (args.out_dir / 'ld.tsv').open('w', newline='') as output:
            writer = csv.writer(output, delimiter='\t')
            writer.writerow(['variant_id', 'r2'])
            for line in handle:
                row = dict(zip(header, line.split()))
                writer.writerow([row['SNP_B'], row['R2']])


if __name__ == '__main__':
    main()
