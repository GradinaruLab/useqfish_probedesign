import subprocess
import os

from Bio import Align, SeqIO, Entrez, SearchIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from primer3 import calc_hairpin     #dependency: primer3-py

import numpy as np
import pandas as pd

import nupack

def designHCR3Probes(gene_id="", gene_name="", hairpin_id=None, email=None, 
                sequence="", db=os.getcwd(), result_path=os.getcwd(), 
                prb_length=20, gc_range=[40, 60], prb_space=1, dg_thresh=-9, 
                spacer = ['ta','at'], to_excel=False):
    """ Creates an Excel file in result_path containing probe designs for a given gene.
    # TODO: fill in descriptions of arguments
    Args:
        gene_id:
        gene_name:
        hairpin_id:
        db:
        result_path:
        prb_space:
        dg_thresh:
        spacer:
        to_excel: bool, whether or not save excel file in result path

    """
    if email:
        Entrez.email=email

    if not sequence:
        # retrieve target sequence from genbank by using accession id
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype = "gb", retmode = "text")
        target = SeqIO.read(handle, "genbank")
        if gene_name.lower() not in target.description.lower():
            print('Target name and accession number are not matched!')
            handle.close()
            exit()
        handle.close()

        # find the coding region of the target
        cds_start = 0
        cds_end = 0
        for feature in target.features:
            if feature.type == 'CDS':
                cds_start = feature.location._start.position
                cds_end = feature.location._end.position
    else:
        target = SeqRecord(
            Seq(sequence.upper()),
            name=gene_name,
        )

        cds_start = 0
        cds_end = len(sequence)

    # create a fasta file including all candidates
    prbs = findAllCandidates(target, prb_length, result_path, gene_name)
    num_prbs = len(prbs)

    ## bowtie2 alignment
    ProbeBowtie2(os.path.join(result_path, f"{gene_name}_prbs_candidates.fasta"),db=db,
                 result_path=os.path.join(result_path, f"{gene_name}_prbs_candidates_alignment_results.sam"))

    ## parse sam file to get mapq
    print(" 1. filtering GC contents, Tm, repeats, dG ...")

    ## basic filtering
    bad_gc, bad_repeats, bad_dg = basicFilter(prbs, num_prbs, prb_length=prb_length, gc_range=gc_range, dg_thresh=dg_thresh)

    # Get full probe sequences and align
    init_seq = GetInitiatorSeq(hairpin_id)
    prbs_full = []
    for i in range(num_prbs):
        full_rec = SeqRecord(init_seq[0:int(len(init_seq)/2)]+spacer[0]+prbs[i].seq[prb_length:prb_length*2], f'{i+1}-1', '', '')
        prbs_full.append(full_rec)
        full_rec = SeqRecord(prbs[i].seq[0:prb_length]+spacer[1]+init_seq[int(len(init_seq)/2):len(init_seq)], f'{i+1}-2', '', '')
        prbs_full.append(full_rec) 
    count = SeqIO.write(prbs_full, os.path.join(result_path, f"{gene_name}_prbs_candidates_full.fasta"), "fasta")
    print("Converted %i records" % count)
    ProbeBowtie2(os.path.join(result_path, f"{gene_name}_prbs_candidates_full.fasta"), db=db,
                 result_path=os.path.join(result_path, f"{gene_name}_prbs_candidates_full_alignment_results.sam"))

    ## Find basic bad probes
    bad_inds = bad_gc + bad_repeats + bad_dg

    # Select probe sequences that are apart
    prb_pos = np.argwhere(np.logical_not(bad_inds))
    prb_pos = np.array(prb_pos.ravel())   
    # prb_pos = [ind for sublist in prb_pos for ind in sublist]
    # print(prb_pos)
    # print(len(prb_pos))

    # Check NCBI databank for description matching
    print(' 2. aligning probe sequences on refseq database using bowtie2')
    bad_unique = IsUnique(os.path.join(result_path, f"{gene_name}_prbs_candidates_alignment_results.sam"), gene_name, prb_pos)
    bad_unique_each = IsUnique(os.path.join(result_path, f"{gene_name}_prbs_candidates_full_alignment_results.sam"),
                               gene_name, prb_pos, full=True)
    prb_pos = prb_pos[~(bad_unique | bad_unique_each)]


    # only select sequences within the coding region
    # and sequences that are apart from the adjacent one
    prb_final_pos = []
    for pos in prb_pos:
        if (pos > cds_start) & (pos+prb_length*2 < cds_end):
            if len(prb_final_pos) == 0:
                prb_final_pos.append(pos)
            elif pos > prb_final_pos[-1] + prb_length*2 + prb_space:
                prb_final_pos.append(pos)
    print(prb_final_pos)
    print('- done! # of final probes: %i' % len(prb_final_pos))

    # recall probe pairs
    prb_final_A = []
    prb_final_B = []
    for pos in prb_final_pos:
        prb_final_A.append(str(init_seq[0:int(len(init_seq)/2)]+spacer[0]+prbs[pos].seq[prb_length:prb_length*2]))
        prb_final_B.append(str(prbs[pos].seq[0:prb_length]+spacer[1]+init_seq[int(len(init_seq)/2):len(init_seq)]))

    result = {'name': gene_name, \
        'accession': gene_id, \
        'position':prb_final_pos, \
        'probe A':prb_final_A, \
        'probe B':prb_final_B}
    resultdf = pd.DataFrame(result)
    if to_excel:
        resultdf.to_excel(excel_writer = os.path.join(result_path, f"{gene_name}_probes.xlsx"))

    return resultdf


def designUSeqFISHProbes(gene_id="", gene_name="", email=None,
                sequence="", db=os.getcwd(), ugi_path=os.getcwd(), 
                ugi_num=1, ugi="", result_path=os.getcwd(), 
                prb_length=20, gc_range=[40, 60], primer_end="TAATGTTATCTT",
                padlock_start="ACATTA", padlock_end="AAGATA", spacer1="attta",
                spacer2 = "atta", prb_space=1, dg_thresh=-9, 
                to_excel=False):
    """ Creates an Excel file in result_path containing probe designs for a given gene.
    # TODO: fill in descriptions of arguments
    Args:
        gene_id:
        gene_name:
        email:
        sequence:
        db:
        result_path:
        prb_space:
        dg_thresh:
        spacer:
        to_excel: bool, whether or not save excel file in result path

    """
    if email:
        Entrez.email = email

    if not ugi:
        ugi_df = pd.read_excel(ugi_path, index_col=0)
        ugi_db = ugi_df['ugi'].values.tolist()
        ugi = ugi_db[ugi_num-1]

    # Get Sequence
    if not sequence:
        # retrieve target sequence from genbank by using accession id
        print("Looking up gene with accession id")
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype = "gb", retmode = "text")
        target = SeqIO.read(handle, "genbank")
        if gene_name.lower() not in target.description.lower():
            print('Target name and accession number are not matched!')
            handle.close()
            exit()
        handle.close()
        # find the coding region of the target
        cds_start = 0
        cds_end = 0
        for feature in target.features:
            if feature.type == 'CDS':
                cds_start = feature.location._start.position
                cds_end = feature.location._end.position
    else:
        target = SeqRecord(
            Seq(sequence.upper()),
            name=gene_name,
        )

        cds_start = 0
        cds_end = len(sequence)

    # Find all probe candidates, create a fasta file
    print("- finding all potential probes ...")
    prbs = findAllCandidates(target, prb_length, result_path, gene_name)
    num_prbs = len(prbs)

    # Bowtie2 alignment
    ProbeBowtie2(os.path.join(result_path, f"{gene_name}_prbs_candidates.fasta"), db=db,
                 result_path=os.path.join(result_path, f"{gene_name}_prbs_candidates_alignment_result.sam"), score_min="G,10,4")

    # basic filtering
    print(" 1. filtering GC contents, Tm, repeats, dG ...") 
    bad_gc, bad_repeats, bad_dg = basicFilter(prbs, num_prbs, prb_length=prb_length, gc_range=gc_range, dg_thresh=dg_thresh)

    # Get full probe sequences (primer and padlock) and align
    primers = []
    padlocks = []
    for i in range(num_prbs):
        primer_rec = SeqRecord(prbs[i].seq[0:prb_length] + primer_end, '%i' % (i+1), '', '')
        padlock_rec = SeqRecord(padlock_start + prbs[i].seq[prb_length:prb_length*2] \
            + spacer1 + ugi + spacer2 + padlock_end, '%i' % (i+1), '', '')
        primers.append(primer_rec)
        padlocks.append(padlock_rec)
    count = SeqIO.write(primers, os.path.join(result_path, f"{gene_name}_primers.fasta"), "fasta")
    count = SeqIO.write(padlocks, os.path.join(result_path, f"{gene_name}_padlocks.fasta"), "fasta")
    print("Converted %i records" % count)

    ProbeBowtie2(os.path.join(result_path, f"{gene_name}_primers.fasta"), db=db,
                 result_path=os.path.join(result_path, f"{gene_name}_primers_candidates_alignment_results.sam"), score_min='G,20,8')
    ProbeBowtie2(os.path.join(result_path, f"{gene_name}_padlocks.fasta"), db=db,
                 result_path=os.path.join(result_path, f"{gene_name}_padlocks_candidates_alignment_results.sam"), score_min='G,20,8')
    
    ## nupack secondary structure analysis
    print(" 2. secondary structure modeling ...")
    bad_secondary = []
    bond_count = []
    for primer, padlock in zip(primers, padlocks):
        bad_primer, bond_count_primer = secondaryFilter(str(primer.seq), part='primer', linker_length=len(primer_end))
        bad_padlock, bond_count_padlock = secondaryFilter(str(padlock.seq), part='padlock', linker_length=len(padlock_start))
        bad_secondary.append(bad_primer | bad_padlock)
        bond_count.append((bond_count_primer, bond_count_padlock))

    # Find bad probes
    bad_inds = bad_gc + bad_repeats + bad_dg + bad_secondary

    # Select probe sequences that are apart
    prb_pos = np.argwhere(np.logical_not(bad_inds))
    prb_pos = np.array(prb_pos.ravel())

    # Check NCBI databank for description matching
    print(' 3. aligning probe sequences on refseq database using bowtie2')
    bad_unique = IsUnique(os.path.join(result_path, f"{gene_name}_prbs_candidates_alignment_result.sam"), gene_name, prb_pos)
    bad_unique_primers = IsUnique(os.path.join(result_path, f"{gene_name}_primers_candidates_alignment_results.sam"), gene_name, prb_pos)
    bad_unique_padlocks = IsUnique(os.path.join(result_path, f"{gene_name}_padlocks_candidates_alignment_results.sam"), gene_name, prb_pos)
    bad_unique_full = bad_unique_primers | bad_unique_padlocks
    prb_pos = prb_pos[~(bad_unique | bad_unique_full)]

    # only select sequences within the coding region
    # and sequences that are apart from the adjacent one
    prb_final_pos = []
    for pos in prb_pos:
        if (pos > cds_start) & (pos+prb_length*2 < cds_end):
            # if len(prb_final_pos) == 0:
            prb_final_pos.append(pos)
            # elif pos > prb_final_pos[-1] + prb_length*2 + prb_space:
            #     prb_final_pos.append(pos)
    print(prb_final_pos)
    print('- done! # of final probes: %i' % len(prb_final_pos))


    # recall probe pairs
    primers_final = []
    padlocks_final = []
    bond_count_final = []
    for pos in prb_final_pos:
        primers_final.append(primers[pos].seq)
        padlocks_final.append('/5Phos/'+padlocks[pos].seq)
        bond_count_final.append(bond_count[pos])

    result = {'name': gene_name, \
        'accession': gene_id, \
        'position':prb_final_pos, \
        'primer':primers_final, \
        'padlock':padlocks_final, \
        'bonds (primer, padlock)':bond_count_final}
    resultdf = pd.DataFrame(result)
    if to_excel:
        resultdf.to_excel(excel_writer = os.path.join(result_path, f"{gene_name}_probes.xlsx"))
    return resultdf


##### HELPER FUNCTIONS #####


def ProbeBowtie2(fastafile, db=os.path.join(os.getcwd(), 'mouse_refseq_rna'), 
                 result_path="prbs_candidates_alignment_result.sam", 
                 score_min='G,20,8'):
    """
    runs bowtie2 to create alignment results from a fasta file
    """
    call = ['bowtie2', '--very-sensitive-local', '-f', '--no-sq', '--no-hd', '--reorder', '--score-min', score_min, \
        '-x', db, '-U', fastafile, '-S', result_path]
    subprocess.check_call(call)
    return


def GetInitiatorSeq(hairpin_id=2, I_id=2):
    init_seqs = [['gAggAgggCAgCAAACgggAAgAgTCTTCCTTTACg', 'gCATTCTTTCTTgAggAgggCAgCAAACgggAAgAg'], \
        ['CCTCgTAAATCCTCATCAATCATCCAgTAAACCgCC', 'AgCTCAgTCCATCCTCgTAAATCCTCATCAATCATC'], \
            ['gTCCCTgCCTCTATATCTCCACTCAACTTTAACCCg', 'AAAgTCTAATCCgTCCCTgCCTCTATATCTCCACTC'], \
                ['CCTCAACCTACCTCCAACTCTCACCATATTCgCTTC', 'CACATTTACAgACCTCAACCTACCTCCAACTCTCAC'], \
                    ['CTCACTCCCAATCTCTATCTACCCTACAAATCCAAT', 'CACTTCATATCACTCACTCCCAATCTCTATCTACCC']]
    return init_seqs[hairpin_id-1][I_id-1]


def IsUnique(samfile_path, gene_name, prb_pos, full=False):
    # find hits (ids) from alignment 
    num_prbs = len(prb_pos)
    hits = [[] for _ in range(num_prbs)]
    if full:
        hits = [[] for _ in range(num_prbs*2)]
    # print(hits)
    with open(samfile_path, "rb") as samfile:
        i = 0
        for line in samfile:
            info = line.decode().split('\t')
            if '-' in info[0]:
                index = int(info[0][:-2])
            else:
                index = int(info[0])-1
            if index not in prb_pos:
                continue
            hits[i].append(info[2])
            i += 1
    
    # search hitted ids if it is variants of the same gene
    variants = ['*']
    bad_unique = np.zeros((num_prbs,), dtype=bool)
    for i, hits_for_oneprb in enumerate(hits):
        for hit_id in hits_for_oneprb:
            if hit_id not in variants:
                handle = Entrez.efetch(db="nucleotide", id=hit_id, rettype="gb", retmode="text")
                search_result = SeqIO.read(handle, "genbank")
                if gene_name.lower() in search_result.description.lower():
                    variants.append(hit_id)
                else:
                    bad_unique[i] = 1

    return bad_unique


def findAllCandidates(target, prb_length, result_path, gene_name):
    """ finds all candidate probes
    """
    prbs = []
    limit = len(target.seq)
    for i in range(0, limit-prb_length*2+1):
        start = i
        end = start + prb_length*2
        temp = target.seq[start:end].reverse_complement()
        temp_rec = SeqRecord(temp, '%i' % (i+1), '', '')
        prbs.append(temp_rec)
    count = SeqIO.write(prbs, os.path.join(result_path, f"{gene_name}_prbs_candidates.fasta"), "fasta")
    print("Converted %i records" % count)

    return prbs


def basicFilter(prbs, num_prbs, prb_length=20, gc_range=[40,60], dg_thresh=-9):
    GC = np.zeros((2, num_prbs))
    repeats = np.zeros((2, num_prbs))
    dg = np.zeros((2, num_prbs))

    for i in range(num_prbs):
        GC[1,i] = (prbs[i].seq[0:prb_length].count("G") \
            + prbs[i].seq[0:prb_length].count("C"))/prb_length *100
        GC[0,i] = (prbs[i].seq[prb_length:prb_length*2].count("G") \
            + prbs[i].seq[prb_length:prb_length*2].count("C"))/prb_length *100

        repeats[1,i] = prbs[i].seq[0:prb_length].count("AAAA") \
            + prbs[i].seq[0:prb_length].count("CCC") \
                + prbs[i].seq[0:prb_length].count("GGG") \
                    + prbs[i].seq[0:prb_length].count("TTTT")
        repeats[0,i] = prbs[i].seq[prb_length:prb_length*2].count("AAAA") \
            + prbs[i].seq[prb_length:prb_length*2].count("CCC") \
                + prbs[i].seq[prb_length:prb_length*2].count("GGG") \
                    + prbs[i].seq[prb_length:prb_length*2].count("TTTT")

        dg[1,i] = calc_hairpin(str(prbs[i].seq[0:prb_length])).dg/1000
        dg[0,i] = calc_hairpin(str(prbs[i].seq[prb_length:prb_length*2])).dg/1000

    bad_gc = (GC < gc_range[0]) | (GC > gc_range[1])
    bad_gc = bad_gc[0,:] | bad_gc[1,:]
    bad_repeats = repeats > 0
    bad_repeats = bad_repeats[0,:] | bad_repeats[1,:]
    bad_dg = dg <= dg_thresh
    bad_dg = bad_dg[0,:] | bad_dg[1,:]
    
    return bad_gc, bad_repeats, bad_dg


def secondaryFilter(seq, part='primer', linker_length=6):
    prb_seq = nupack.Strand(seq, name='prb')
    model = nupack.Model(material='dna', celsius=37, sodium=0.39, magnesium=0.0)
    set = nupack.ComplexSet(strands=[prb_seq])
    result = nupack.complex_analysis(complexes=set, model=model, compute=['mfe'])
    secondstruct = str(result['(prb)'].mfe[0].structure)

    bad = True
    bond_count = secondstruct.count('(')
    if bond_count == 0:
            bad = False
    else:
        secondstruct = secondstruct.replace(')','(').replace('+','(')
        splitted = secondstruct.split('(')
        if (part=='primer') and (splitted[-1].count('.')>linker_length):
            bad = False
        if (part=='padlock') and (splitted[0].count('.')>linker_length) and (splitted[-1].count('.')>linker_length):
            bad = False

    return bad, bond_count

