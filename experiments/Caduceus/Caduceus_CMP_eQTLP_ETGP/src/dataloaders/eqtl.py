import os,sys
import torch
import tabix
import numpy as np
import pandas as pd
from tqdm import tqdm
import pyfaidx
import kipoiseq
from kipoiseq import Interval


data_path = "data_long_range_dna/eQTL"

HYENADNA_ALPHABET = {'A': 7, 'C': 8, 'G': 9, 'T': 10, 'N': 11}

_DNA_TOKEN_LOOKUP = np.full(256, 6, dtype=np.int64)  # UNK
for _base, _token_id in HYENADNA_ALPHABET.items():
    _DNA_TOKEN_LOOKUP[ord(_base)] = _token_id
    _DNA_TOKEN_LOOKUP[ord(_base.lower())] = _token_id

def encode_hyenadna(sequence):
    raw = np.frombuffer(
        sequence.encode('ascii', errors='replace'),
        dtype=np.uint8,
    )
    return _DNA_TOKEN_LOOKUP[raw]


class FastaStringExtractor:
    """
    Extract sequences from a fasta file.
    """

    def __init__(self, fasta_file):
        self.fasta = pyfaidx.Fasta(fasta_file)
        self._chromosome_sizes = {k: len(v) for k, v in self.fasta.items()}

    def extract(self, interval: Interval, **kwargs) -> str:
        # Truncate interval if it extends beyond the chromosome lengths.
        chromosome_length = self._chromosome_sizes[interval.chrom]
        trimmed_interval = Interval(interval.chrom,
                                    max(interval.start, 0),
                                    min(interval.end, chromosome_length),
                                    )
        # pyfaidx wants a 1-based interval
        sequence = str(self.fasta.get_seq(trimmed_interval.chrom,
                                          trimmed_interval.start + 1,
                                          trimmed_interval.stop).seq).upper()
        # Fill truncated values with N's.
        pad_upstream = 'N' * max(-interval.start, 0)
        pad_downstream = 'N' * max(interval.end - chromosome_length, 0)
        return pad_upstream + sequence + pad_downstream

    def close(self):
        return self.fasta.close()


class EQTLseqDataSet(torch.utils.data.IterableDataset):

    def __init__(self, config_file, subset):
        #########################
        # load config
        #########################
        self.config = parse_config(config_file)
        # init
        self.organism = 'Unknown'
        if self.config.get('organism', None) is not None:
            self.organism = self.config['organism']
        self.subset = subset
        if self.subset is None:
            if self.config.get('subset', None) is not None:
                self.subset = self.config['subset']
            else:
                print("subset is required either in the config file or passed as a parameter", file = sys.stderr)
                exit(1)
        try:
            assert self.subset in ['train', 'valid', 'test']
        except AssertionError:
            print("subset must be either train, valid, or test", file = sys.stderr)
            exit(1)
        # verbose level
        self.verbose = 0
        if self.config.get('verbose', None) is not None:
            self.verbose = self.config['verbose']
        # show config
        if self.verbose > 0:
            print_config(self.config)
        print("> load config done")

        ############################
        # load fasta file
        ############################
        if self.config.get('genome_fa', None) is None:
            print('genome_fa is required in the config file', file = sys.stderr)
            exit(1)
        self.fasta_reader = FastaStringExtractor(os.path.join(data_path, self.config['genome_fa']))
        print("> init fasta extractor done")

        ############################
        # load dataset, ie eQTL data
        ############################
        if self.config.get('seq_len_cutoff', None) is None:
            print('seq_len_cutoff is required in the config file', file = sys.stderr)
            print('will use default value 450000', file = sys.stderr)
            self.config['seq_len_cutoff'] = 450000
        if self.config.get('eQTL_tabix_file', None) is None:
            print('eQTL_tabix_file is required in the config file', file = sys.stderr)
            exit(1)
        self.dataset = parse_eQTL(self.config['eQTL_file'], self.config['eQTL_tabix_file'], self.fasta_reader, self.subset, self.config['seq_len_cutoff'], self.config['tss_flank_upstream'], self.config['tss_flank_downstream'], self.config['region_flank_upstream'], self.config['region_flank_downstream'])
        if self.verbose > 0:
            print(f'Loaded {len(self.dataset)} eQTL records')

    def __iter__(self):
        # iterate over the data
        target2int = {'positive':1, 'negative':0}
        for index, record in enumerate(self.dataset):
            sequence_ref = encode_hyenadna(record[0])
            sequence_alt = encode_hyenadna(record[1])
            target = np.array(target2int[record[2]]).astype('int8')
            yield (sequence_ref, sequence_alt, target)
            # yield {
            #     'sequence_ref': sequence_ref,
            #     'sequence_alt': sequence_alt,
            #     'target': target
            # }


def rcDNA(seq):
    """Get reverse complemented of a DNA sequence"""
    diccomseq = {'A':'T', 'T':'A', 'G':'C', 'C':'G', 'N':'N'}
    return "".join(diccomseq[nu.upper()] for nu in seq[::-1])


def one_hot_encode(sequence):
    return kipoiseq.transforms.functional.one_hot_dna(sequence).astype(np.int8)


def parse_config(config_file):
    """
    load config as a dictionary
    """
    config_table = {}
    with open(config_file, 'r') as fin:
        for line in fin:
            row = line.strip().split()
            key = row[0]
            value = row[1]
            data_type = row[2]
            if data_type == 'int':
                value = int(value)
            elif data_type == 'float':
                value = float(value)
            elif data_type == 'bool':
                value = bool(value)
            elif data_type == 'str':
                value = str(value)
            elif data_type == 'list':
                value_list = value.split(',')
                # get the data type of the list
                new_value_list = []
                for value in value_list:
                    if value.isnumeric():
                        value = int(value)
                    elif value.replace('.','',1).isdigit():
                        value = float(value)
                    else:
                        value = str(value)
                    new_value_list.append(value)
                value = new_value_list
            else:
                value = str(value)
            #
            config_table[key] = value
    return config_table


def print_config(config_table):
    """
    print config
    """
    for key, value in config_table.items():
        print(f'{key}: {value}')


def parse_eQTL(eQTL_file, blacklist_tabix_file, fasta_reader, subset_group_name, seq_len_cutoff = 450000, tss_flank_upstream = 3000, tss_flank_downstream = 3000, region_flank_upstream = 500, region_flank_downstream = 500):
    """
    this function will parse the following file format to reresent an variant-gene interaction entry
        > gene_chrom  gene_start gene_end region_chrom region_start region_end gene_id gene_strand region_id target subset allele1 allele2
        > other columns may exists but will be ignored
    then, the function will search for any variants between the tested region and the gene tss the sequence in these variants will be masked N. This is crucial because we want to ensure that these sequences don’t influence the model training process.
    """
    # check subset
    if subset_group_name not in ['train', 'valid', 'test']:
        print(f'Error: unknown subset group {subset_group_name}')
        exit(1)
    # load records as dataframe
    df = pd.read_csv(os.path.join(data_path, eQTL_file), sep='\t', header=0)
    # load the tabix
    blacklist_tabix = tabix.open(os.path.join(data_path, blacklist_tabix_file))
    # parse the records and convert it into a list of DNA sequence and target (ie, true or false eQTL)
    dataset = []
    print("> Start parsing eQTL records to build the dataset %s" % (subset_group_name))
    N_total = df.shape[0]
    N_pass = 0
    N_skip_for_chrom_diff = 0 # skip if variant and gene tss are in different chromosomes
    N_skip_for_len = 0 # skip if the distance between variant and gene tss is too long
    N_skip_for_strand = 0 # skip if gene strand is not + or - or unknown
    for index, row in tqdm(df.iterrows(), total = N_total):
        gene_chrom = row['gene_chrom']
        gene_start = row['gene_start']
        gene_end = row['gene_end']
        region_chrom = row['region_chrom']
        region_start = row['region_start']
        region_end = row['region_end']
        gene_id = row['gene_id']
        gene_strand = row['gene_strand']
        region_id = row['region_id']
        target = row['target']
        subset_group = row['subset']
        allele1 = row['allele1']
        allele2 = row['allele2']
        # skip if in different subset group
        if subset_group != subset_group_name:
            continue
        # skip if in different chromosomes
        if gene_chrom != region_chrom:
            N_skip_for_chrom_diff += 1
            continue
        # get tss
        if gene_strand == '+':
            tss_start = gene_start
        elif gene_strand == '-':
            tss_start = gene_end - 1
        else:
            print(f'Warning: unknown strand {gene_strand}')
            N_skip_for_strand += 1
            continue
        tss_end = tss_start + 1
        # add flanking sequence
        tss_start = tss_start - tss_flank_upstream
        tss_end = tss_end + tss_flank_downstream
        variant_start = region_start
        variant_end = region_end
        region_start = region_start - region_flank_upstream
        region_end = region_end + region_flank_downstream
        # calculate the distance between tss and region
        distance = max(0, max(tss_start, region_start) - min(tss_end, region_end))
        # check if the sequence length is too long
        if distance > seq_len_cutoff:
            N_skip_for_len += 1
            continue
        # get the interverl between the region and tss
        sequence_start = min(tss_start, region_start)
        sequence_end = max(tss_end, region_end)
        # get the sequence
        region_interval = Interval(region_chrom, sequence_start, sequence_end)
        region_seq = fasta_reader.extract(region_interval)
        # check if the variant sequence matches the reference sequence
        variant_rev_start = variant_start - sequence_start
        variant_rev_end = variant_end - sequence_start
        assert region_seq[variant_rev_start:variant_rev_end] == row['allele1']
        # mask any eQTL sequence
        if distance > 0:
            if tss_start > region_end:
                query_start = region_end
                query_end = tss_start
            else:
                query_start = tss_end
                query_end = region_start
            query_result = blacklist_tabix.query(region_chrom, query_start, query_end)
            for overlap in query_result:
                o_start = int(overlap[1])
                o_end = int(overlap[2])
                # mask the correspoinding DNA sequence as 'N'
                rev_start = o_start - sequence_start
                rev_end = o_end - sequence_start
                if rev_start > 0 and rev_start < len(region_seq) and rev_end > 0 and rev_end < len(region_seq):
                    region_seq = region_seq[:rev_start] + 'N' * (rev_end - rev_start) + region_seq[rev_end:]
        # create variant sequence
        region_seq_var =  region_seq[:variant_rev_start] + row['allele2'] + region_seq[variant_rev_end:]
        # flip (reverse complement) the sequence if gene is on the downstream of region (ie, genomic coordinate is larger)
        if gene_start > region_end:
            region_seq = rcDNA(region_seq)
            region_seq_var = rcDNA(region_seq_var)
        # padding N to the sequence if the sequence
        if len(region_seq) <= seq_len_cutoff:
            region_seq = region_seq + 'N' * (seq_len_cutoff - len(region_seq))
        else:
            region_seq = region_seq[:seq_len_cutoff]
        if len(region_seq_var) <= seq_len_cutoff:
            region_seq_var = region_seq_var + 'N' * (seq_len_cutoff - len(region_seq_var))
        else:
            region_seq_var = region_seq_var[:seq_len_cutoff]
        N_pass +=1
        # append to the dataset
        
        dataset.append((region_seq, region_seq_var, target))
    # report the statistics
    print("# Finish parsing eQTL records")
    print("# Total records: ", N_total)
    print("# Skipped records due to different chromosomes: ", N_skip_for_chrom_diff)
    print("# Skipped records due to distance cutoff: ", N_skip_for_len)
    print("# Skipped records due to unknown strand: ", N_skip_for_strand)
    print("# Select records %s with subset %s " % (N_pass, subset_group_name))
    return dataset

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# load dataset
# change to train, valid, or test to load different dataset
train_dataset = EQTLseqDataSet("{}/config/gtex_hg38.Adipose_Subcutaneous.config".format(data_path), 'valid')
train_loader = torch.utils.data.DataLoader(train_dataset, num_workers=0, batch_size=1)
for i, batch in enumerate(train_loader):
    seq_ref, seq_alt, target = batch
        # batch_gpu = {k:v.to(device) for k,v in batch.items()}
        # seq_ref = batch_gpu['sequence_ref']
        # seq_alt = batch_gpu['sequence_alt']
        # target = batch_gpu['target']
    print(seq_ref.shape)
    print(seq_ref[:, 0])
    print(seq_alt.shape)
    print(target.shape)
    print(target)
    break
# should show the shape of the sequence and target
# torch.Size([5, 4507000, 4])
# torch.Size([5])
