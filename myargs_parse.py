import argparse

def get_args():
    p = argparse.ArgumentParser(
        description='Endmotifs stat required arguments!')
    p.add_argument('-i', '--indir', required=True, help='indir')
    p.add_argument(
        '-sl',
        '--samples',
        required=True,
        help='Sample name list:sample_name.txt')
    p.add_argument(
        '-n',
        '--num_base',
        type=int,
        default=4,
        help='Base num, one base motif: A, T, G, C.')
    p.add_argument(
        '-r',
        '--ref',
        default=None,
        help='Chromosome,for example, chrM, chr1,chr2...')
    p.add_argument(
        '-se',
        '--start_end',
        type=int,
        nargs=2,
        help=
        'Start position: end position. If not provided, the default analysis will be performed on the entire reference.',
        default=[None, None])
    p.add_argument(
        '-o',
        '--outdir',
        required=True,
        help='outdir,for example:/mnt/data/Seq_data/20220726')
    p.add_argument(
        '-fa',
        '--fasta',
        default=None,
        help=
        'FASTA path, for example: /mnt/data/Seq_data/20220726/hg38.fa. Default is None. If not provided, endmotifs counts will not be standardized.'
    )
    p.add_argument(
        '-suff',
        '--suffix',
        default="bam",
        help='bam files suffix,if bam format,index before run this script.')
    p.add_argument(
        '-v',
        '--vref',
        action='store_true',
        default=False,
        help='Analyze chromosomes other than those specified')
    args = p.parse_args()
    return args