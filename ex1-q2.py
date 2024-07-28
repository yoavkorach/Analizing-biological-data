import re


def is_valid_DNA_sequence(seq):
    """
    :param seq: a sequence to test
    :return: True if a valid DNA sequence, else False
    """
    pattern = r"^[ATCG]*$"
    return bool(re.match(pattern, seq))


def find_start_codons(seq):
    """
    :param seq: a sequence to test
    :return: a list of all the start codons positions
    """
    start = r'(ATG)'
    return [iter.start(0) for iter in re.finditer(start, seq)]


def find_stop_codons(seq):
    """
    :param seq: a sequence to test
    :return: a list of all the stop codons positions
    """
    patterns = r'T(AG|AA|GA)'
    lst = [m.start() for m in re.finditer(patterns, seq)]
    return lst


def extract_coding_sequences(seq):
    """
        :param seq: a sequence to test
        :return: a list of all the coding sequences in the sequences
    """
    pattern = re.compile(r'(?=(ATG(?:[ACGT]{3})*(?:TAA|TAG|TGA)))')
    matches = pattern.findall(seq)
    return [match[3:-3] for match in matches if len(match) % 3 == 0]


def translate_dna_to_protein(seq):
    """
        :param seq: a sequence to test
        :return: a string, the resulting aa sequence
    """
    ans = ""
    for i in range(0, len(seq), 3):
        ans += codon_to_aa[seq[i:i + 3]]
    return ans


def find_nglycosylation_sites(seq):
    """
        :param seq: a sequence to test
        :return: True if contains nglycosylation_sites, False otherwise
    """
    pattern = re.compile(r'N[^P][ST][^P]')
    return bool(pattern.search(seq))


def dna_to_nglycosylation(seq):
    """
        :param seq: a sequence to test
        :return: number of protein in the sequence that contain N-glycosylation sites
    """
    if (not is_valid_DNA_sequence(seq)):
        raise ValueError
    lst = extract_coding_sequences(seq)
    cnt = 0
    for dna in lst:
        pro = translate_dna_to_protein(dna)
        if find_nglycosylation_sites(pro):
            cnt += 1
    return cnt


codon_to_aa = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}

if __name__ == '__main__':
    # you can write whatever you want here
    pass


