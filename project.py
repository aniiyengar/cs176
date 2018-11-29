""" 
    RNA Alignment Assignment
    
    Implement each of the functions below using the algorithms covered in class.
    You can construct additional functions and data structures but you should not
    change the functions' APIs.

    You will be graded on the helper function implementations as well as the RNA alignment, although
    you do not have to use your helper function.
    
    *** Make sure to comment out any print statement so as not to interfere with the grading script
"""

import sys # DO NOT EDIT THIS
from shared import *

import time
import multiprocessing as mp

ALPHABET = [TERMINATOR] + BASES

def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K(object):
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0  
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

def get_suffix_array(s):
    """
    Naive implementation of suffix array generation (0-indexed). You do not have to implement the
    KS Algorithm. Make this code fast enough so you have enough time in Aligner.__init__ (see bottom).

    Input:
        s: a string of the alphabet ['A', 'C', 'G', 'T'] already terminated by a unique delimiter '$'
    
    Output: list of indices representing the suffix array

    >>> get_suffix_array('GATAGACA$')
    [8, 7, 5, 3, 1, 6, 4, 0, 2]
    """
    def compare(i, j):
        c = 0
        while s[i + c] == s[j + c]:
            c += 1
        if s[i + c] > s[j + c]:
            return 1
        elif s[i + c] < s[j + c]:
            return -1
        else:
            return 0
    
    buckets = {}
    for i in range(len(s)):
        if s[i:i+100] in buckets:
            buckets[s[i:i+100]].append(i)
        else:
            buckets[s[i:i+100]] = [i]

    keys = list(buckets.keys())
    for i in range(len(keys)):
        item = keys[i]
        if len(buckets[item]) > 1:
            buckets[item].sort(key=cmp_to_key(compare))

    result = []
    for item in sorted(keys):
        result.extend(buckets[item])

    return result

def get_bwt(s, sa):
    """
    Input:
        s: a string terminated by a unique delimiter '$'
        sa: the suffix array of s

    Output:
        L: BWT of s as a string
    """
    return ''.join([s[i-1] for i in sa])

def get_F(L):
    """
    Input: L = get_bwt(s)

    Output: F, first column in Pi_sorted
    """
    return ''.join(sorted(list(L)))

def get_M(F):
    """
    Returns the helper data structure M (using the notation from class). M is a dictionary that maps character
    strings to start indices. i.e. M[c] is the first occurrence of "c" in F.

    If a character "c" does not exist in F, you may set M[c] = -1
    """
    M = {}
    for i in range(len(F)):
        char = F[i]
        if char not in M:
            M[char] = i
    for char in ALPHABET:
        if char not in M:
            M[char] = -1
    return M

def get_occ(L):
    """
    Returns the helper data structure OCC (using the notation from class). OCC should be a dictionary that maps 
    string character to a list of integers. If c is a string character and i is an integer, then OCC[c][i] gives
    the number of occurrences of character "c" in the bwt string up to and including index i
    """
    occ = {i: [0] for i in ALPHABET}
    for i in L:
        for j in occ:
            if i != j:
                occ[j].append(occ[j][-1])
            else:
                occ[i].append(occ[i][-1]+1)
    for i in occ:
        occ[i].pop(0)
    return occ

def construct_L(M, occ):
    length = 0
    for char in occ:
        length += occ[char][-1]
    L = ['']*length
    for char in occ:
        for i in range(occ[char][-1]):
            L[occ[char].index(i+1)] = char
    return ''.join(L)

def exact_suffix_matches(p, M, occ):
    """
    Find the positions within the suffix array sa of the longest possible suffix of p 
    that is a substring of s (the original string).
    
    Note that such positions must be consecutive, so we want the range of positions.

    Input:
        p: the pattern string
        M, occ: buckets and repeats information used by sp, ep

    Output: a tuple (range, length)
        range: a tuple (start inclusive, end exclusive) of the indices in sa that contains
            the longest suffix of p as a prefix. range=None if no indices matches any suffix of p
        length: length of the longest suffix of p found in s. length=0 if no indices matches any suffix of p

        An example return value would be ((2, 5), 7). This means that p[len(p) - 7 : len(p)] is
        found in s and matches positions 2, 3, and 4 in the suffix array.

    >>> s = 'ACGT' * 10 + '$'
    >>> sa = get_suffix_array(s)
    >>> sa
    [40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3]
    >>> L = get_bwt(s, sa)
    >>> L
    'TTTTTTTTTT$AAAAAAAAAACCCCCCCCCCGGGGGGGGGG'
    >>> F = get_F(L)
    >>> F
    '$AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTT'
    >>> M = get_M(F)
    >>> sorted(M.items())
    [('$', 0), ('A', 1), ('C', 11), ('G', 21), ('T', 31)]
    >>> occ = get_occ(L)
    >>> type(occ) == dict, type(occ['$']) == list, type(occ['$'][0]) == int
    (True, True, True)
    >>> occ['$']
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    >>> exact_suffix_matches('ACTGA', M, occ)
    ((1, 11), 1)
    >>> exact_suffix_matches('$', M, occ)
    ((0, 1), 1)
    >>> exact_suffix_matches('AA', M, occ)
    ((1, 11), 1)
    """
    
    # initialize sp, ep
    length = len(p)
    last_char = p[-1]
    sp = M[last_char]
    if sp == -1:
        return (None, 0)
    
    # find next(last_char)
    nxt = float('inf')
    nxt_item = None
    for item in M:
        if nxt > M[item] > M[last_char] and item != last_char:
            nxt = M[item]
            nxt_item = item
    if nxt_item == None:
        ep = len(occ['$']) - 1
    else:
        ep = nxt - 1
    
    # changed for loop a bit, only works on strings len >= 2
    # for strings of len 1, it skips to the end and works
    for i in range(length-2,-1,-1):
        sp_ph = M[p[i]] + occ[p[i]][sp-1]
        ep_ph = M[p[i]] + occ[p[i]][ep]-1
        if sp_ph > ep_ph:
            return ((sp,ep+1),length-i-1)
        sp = sp_ph
        ep = ep_ph
    return ((sp,ep+1), length)

def bowtie(p, M, occ):
    num_mismatches = 0
    length = len(p)
    last_char = p[-1]
    sp = M[last_char]
    if sp == -1:
        return (None, 0)
    
    # find next(last_char)
    nxt = float('inf')
    nxt_item = None
    for item in M:
        if nxt > M[item] > M[last_char] and item != last_char:
            nxt = M[item]
            nxt_item = item
    if nxt_item == None:
        ep = len(occ['$']) - 1
    else:
        ep = nxt - 1
    
    first_mismatch = -1
    should_replace = False
    # changed for loop a bit, only works on strings len >= 2
    # for strings of len 1, it skips to the end and works
    for i in range(length-2,-1,-1):
        sp_ph = M[p[i]] + occ[p[i]][sp-1]
        ep_ph = M[p[i]] + occ[p[i]][ep]-1

        poss_switches = ALPHABET.copy()
        poss_switches.remove(p[i])
        is_mismatch = False
        while sp_ph > ep_ph and poss_switches:
            new_char = poss_switches.pop()
            sp_ph = M[new_char] + occ[new_char][sp-1]
            ep_ph = M[new_char] + occ[new_char][ep]-1
            is_mismatch = True

            if first_mismatch == -1:
                first_mismatch = i

        num_mismatches += is_mismatch
        sp = sp_ph
        ep = ep_ph

    if num_mismatches > 6:
        q = mp.Queue()
        best_align_score = num_mismatches
        best_align = ((sp,ep+1), num_mismatches)
        best_ix = None
        ps = [mp.Process(target=bowtie_offset_queue, args=(p, M, occ, q, i, best_align, best_align_score)) for i in range(8)]
        for proc in ps:
            proc.start()
        for i in range(8):
            ix, j = q.get()
            if j[1] < best_align_score:
                best_align_score = j[1]
                best_align = j
                best_ix = ix
        return best_align
    else:
        return ((sp,ep+1), num_mismatches)

def bowtie_offset(p, M, occ, mismatches):
    num_mismatches = 0
    length = len(p)
    last_char = p[-1]
    sp = M[last_char]
    if sp == -1:
        return (None, 0)
    
    # find next(last_char)
    nxt = float('inf')
    nxt_item = None
    for item in M:
        if nxt > M[item] > M[last_char] and item != last_char:
            nxt = M[item]
            nxt_item = item
    if nxt_item == None:
        ep = len(occ['$']) - 1
    else:
        ep = nxt - 1

    # changed for loop a bit, only works on strings len >= 2
    # for strings of len 1, it skips to the end and works
    for i in range(length-2,-1,-1):
        if i in mismatches:
            # force a mismatch here
            poss_switches = list('$CGAT')
            poss_switches.remove(p[i])
            sp_ph = 1
            ep_ph = 0
            while sp_ph > ep_ph and poss_switches:
                new_char = poss_switches.pop()
                sp_ph = M[new_char] + occ[new_char][sp-1]
                ep_ph = M[new_char] + occ[new_char][ep]-1

            num_mismatches += 1
        else:
            sp_ph = M[p[i]] + occ[p[i]][sp-1]
            ep_ph = M[p[i]] + occ[p[i]][ep]-1

            poss_switches = ALPHABET.copy()
            poss_switches.remove(p[i])
            is_mismatch = False
            while sp_ph > ep_ph and poss_switches:
                new_char = poss_switches.pop()
                sp_ph = M[new_char] + occ[new_char][sp-1]
                ep_ph = M[new_char] + occ[new_char][ep]-1
                is_mismatch = True

            num_mismatches += is_mismatch
            sp = sp_ph
            ep = ep_ph

    return ((sp,ep+1), num_mismatches)

def suffix_replacer(hits):
    (start, end), length = hits
    genome_seqs = []
    for i in range(start,end):
        index = a.genome_sa[i]
        genome_seqs.append((index,index+length))
    return genome_seqs

MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 10000

class Aligner: 
    def __init__(self, genome_sequence, known_genes):
        """
        Initializes the aligner. Do all time intensive set up here. i.e. build suffix array.

        genome_sequence: a string (NOT TERMINATED BY '$') representing the bases of the of the genome
        known_genes: a python set of Gene objects (see shared.py) that represent known genes. You can get the isoforms 
                     and exons from a Gene object

        Time limit: 500 seconds maximum on the provided data. Note that our server is probably faster than your machine, 
                    so don't stress if you are close. Server is 1.25 times faster than the i7 CPU on my computer

        """
        with open('genome.fa') as f:
            genome_sequence = f.read().split()[1].strip() + '$'
        last_time = time.time()
        self.genome_sa = get_suffix_array(genome_sequence)
        self.genome_L = get_bwt(genome_sequence, self.genome_sa)
        self.genome_M = get_M(get_F(self.genome_L))
        self.genome_occ = get_occ(self.genome_L)
        
        print('FM index of genome: ' + str(time.time() - last_time))
        last_time = time.time()
        
        isos_sa = {}
        isos_L = {}
        isos_M = {}
        isos_occ = {}
        isos = {}
        iso_exon_starts = {}
        iso_names = []
        known_genes = []
        with open('genes.tab') as f:
            curr = f.readline().split()
            gene_id = curr[1]
            isoform_id = ''
            isoforms = []
            exons = []
            first_gene = -1
            while curr:
                if curr[0] == 'gene':
                    if exons:
                        isoforms.append(Isoform(isoform_id,exons))
                    if first_gene != -1:
                        known_genes.append(Gene(gene_id, isoforms))
                    else:
                        first_gene = 0
                    gene_id = curr[1]
                    isoforms = []
                    first_iso = -1
                elif curr[0] == 'isoform':
                    if first_iso != -1:
                        isoforms.append(Isoform(isoform_id, exons))
                    else:
                        first_iso = 0
                    isoform_id = curr[1]
                    exons = []
                elif curr[0] == 'exon':
                    exons.append(Exon(curr[1],int(curr[2]),int(curr[3])))
                curr = f.readline().split()
        
        for gene in known_genes:
            for isoform in gene.isoforms:
                iso_names.append(isoform.id)
                spliced_isoform = ''
                i = isoform.id
                isos[i] = isoform
                prev = 0
                iso_exon_starts[i] = []
                curr = 0
                for exon in isoform.exons:
                    exon_str = genome_sequence[exon.start:exon.end]
                    spliced_isoform += exon_str
                    curr += len(exon_str)
                    iso_exon_starts[i].append(curr)
                spliced_isoform += '$'
                isos_sa[i] = get_suffix_array(spliced_isoform)
                isos_L[i] = get_bwt(spliced_isoform, isos_sa[i])
                isos_M[i] = get_M(get_F(isos_L[i]))
                isos_occ[i] = get_occ(isos_L[i])
                print('FM index of isoform, length %d: ' % len(spliced_isoform) + str(time.time() - last_time))
                last_time = time.time()
        
        self.isos_sa = isos_sa
        self.isos_L = isos_L
        self.isos_M = isos_M
        self.isos_occ = isos_occ
        self.iso_names = iso_names
        self.isos = isos
        self.iso_exon_starts = iso_exon_starts

    def align(self, read_sequence):
        """
        Returns an alignment to the genome sequence. An alignment is a list of pieces. 
        Each piece consists of a start index in the read, a start index in the genome, and a length 
        indicating how many bases are aligned in this piece. Note that mismatches are count as "aligned".

        Note that <read_start_2> >= <read_start_1> + <length_1>. If your algorithm produces an alignment that 
        violates this, we will remove pieces from your alignment arbitrarily until consecutive pieces 
        satisfy <read_start_2> >= <read_start_1> + <length_1>

        Return value must be in the form (also see the project pdf):
        [(<read_start_1>, <reference_start_1, length_1), (<read_start_2>, <reference_start_2, length_2), ...]

        If no good matches are found: return the best match you can find or return []

        Time limit: 0.5 seconds per read on average on the provided data.
        """
        
        # bowtie, read to isoforms
        # first see if the read matches a known isoform for under 6 mismatches
        def suffix_replacer(a, hits):
            (start, end), length = hits
            genome_seqs = []
            for i in range(start,end):
                index = self.genome_sa[i]
                genome_seqs.append((index,index+length))
            return genome_seqs

        min_mismatches = float('inf')
        best_match = None
        for iso_name in self.iso_names:
            match_range, num_mismatches = bowtie(read_sequence, self.isos_M[iso_name], self.isos_occ[iso_name])
            if num_mismatches < min_mismatches:
                min_mismatches = num_mismatches
                best_match = (iso_name, match_range)
        
        alignment = []
        seqlength = len(read_sequence)
        
        #Transcriptome Alignment
        if min_mismatches <= 6 and best_match:
            # we have a match, direct to transcriptome
            # if we have multiple matches, it'll be weird but just take the first one
            iso_name, match_range = best_match
            curr = self.isos_sa[iso_name][match_range[0]] # start pt in spliced isoform
            print(curr)
            exon_starts = self.iso_exon_starts[iso_name]
            exon_num = 0
            curr_pt = 0
            for i in range(len(exon_starts)):
                if curr < exon_starts[i]:
                    exon_num = i
                    break
            exons = self.isos[iso_name].exons
            print(exons)
            print(exon_starts)
            read_pt = 0
            curr_pt = exons[exon_num].start + curr
            while seqlength:
                offset = exons[exon_num].end-curr_pt
                diff = min(offset, seqlength)
                alignment.append((read_pt, curr_pt, diff))
                seqlength -= diff
                curr_pt += diff
                read_pt += diff
                exon_num += 1
            return alignment
        #Genome Alignment
        else:
            min_mis = 7
            best_hits = []
            hits = bowtie(read, self.genome_M, self.genome_occ)
            curr_mis = hits[1]
            if curr_mis < min_mis:
                min_mis = curr_mis
                best_hits = [(0, suffix_replacer(hits)[0][0], 50)]
            for offset in range(10, 41):
                seed1, seed2 = read[:offset], read[offset:]
                hits1 = bowtie(seed1, self.genome_M, self.genome_occ)
                hits2 = bowtie(seed2, self.genome_M, self.genome_occ)
                curr_mis = hits1[1] + hits2[1]
                if curr_mis < min_mis:
                    hits1 = suffix_replacer(hits1)
                    hits2 = suffix_replacer(hits2)
                    for seq1 in hits1:
                        for seq2 in hits2:
                            if MIN_INTRON_SIZE < seq2[0] - seq1[1] < MAX_INTRON_SIZE:
                                min_mis = curr_mis
                                seq1len = len(seq1)
                                seq2len = len(seq2)
                                best_hits = [(0, seq1[0], seq1len), (seq1len, seq2[0], seq2len)]
                                break
            for offset1 in range(10,31):
                for offset2 in range(10+offset1, 41):
                    seed1, seed2,seed3 = read[:offset1], read[offset1:offset2], read[offset2:]
                    hits1 = bowtie(seed1, self.genome_M, self.genome_occ)
                    hits2 = bowtie(seed2, self.genome_M, self.genome_occ)
                    hits3 = bowtie(seed3, self.genome_M, self.genome_occ)
                    curr_mis = hits1[1] + hits2[1] + hits3[1]
                    if curr_mis < min_mis:
                        hits1 = suffix_replacer(hits1)
                        hits2 = suffix_replacer(hits2)
                        hits3 = suffix_replacer(hits3)
                        for seq1 in hits1:
                            for seq2 in hits2:
                                for seq3 in hits3:
                                    if MIN_INTRON_SIZE < seq2[0] - seq1[1] < MAX_INTRON_SIZE and MIN_INTRON_SIZE < seq3[0] - seq2[1] < MAX_INTRON_SIZE:
                                        min_mis = curr_mis
                                        seq1len = len(seq1)
                                        seq2len = len(seq2)
                                        seq3len = len(seq3)
                                        best_hits = [(0, seq1[0], seq1len), (seq1len, seq2[0], seq2len), (seq1len+seq2len, seq3[0], seq3len)]
                                        break
            return best_hits
        return []
