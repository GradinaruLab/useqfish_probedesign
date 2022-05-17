import subprocess
import os

from Bio import Align, SeqIO, Entrez, SearchIO, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from primer3 import calcHairpin     #dependency: primer3-py

import numpy as np
import pandas as pd
from sqlalchemy.orm.query import AliasOption


def designHCR3Probes(gene_id="", gene_name="", hairpin_id=None, email=None,
                db=os.getcwd(), result_path=os.getcwd(), 
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
    # retrieve target sequence from genbank by using accession id
    handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype = "gb", retmode = "text")
    target = SeqIO.read(handle, "genbank")
    if gene_name not in target.description.lower():
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

    # create a fasta file including all candidates
    prbs = findAllCandidates(target, prb_length, result_path)
    num_prbs = len(prbs)

    ## bowtie2 alignment
    ProbeBowtie2(os.path.join(result_path, "prbs_candidates.fasta"),db=db, result_path=os.path.join(result_path, "alignment_results.sam"))

    ## parse sam file to get mapq and
    ## find only unique probe sequences
    print(" 0. aligning probe sequences on refseq database using bowtie2")
    is_unique = IsUnique(os.path.join(result_path, "alignment_results.sam"), gene_name, num_prbs)
    bad_unique = np.logical_not(is_unique)


    ## basic filtering
    bad_gc, bad_repeats, bad_dg = basicFilter(prbs, num_prbs, prb_length=prb_length, gc_range=gc_range, dg_thresh=dg_thresh)


    # Get full probe sequences and align
    init_seq = GetInitiatorSeq(hairpin_id)
    prbs_full = []
    for i in range(num_prbs):
        full_rec = SeqRecord(init_seq[0:int(len(init_seq)/2)]+spacer[0]+prbs[i].seq[prb_length:prb_length*2], '%i-1' % (i+1), '', '')
        prbs_full.append(full_rec)
        full_rec = SeqRecord(prbs[i].seq[0:prb_length]+spacer[1]+init_seq[int(len(init_seq)/2):len(init_seq)], '%i-2' % (i+1), '', '')
        prbs_full.append(full_rec) 
    count = SeqIO.write(prbs_full, os.path.join(result_path, "prbs_candidates_full.fasta"), "fasta")
    print("Converted %i records" % count)
    ProbeBowtie2(os.path.join(result_path, "prbs_candidates_full.fasta"), db=db, result_path=os.path.join(result_path, "prbs_candidates_full_alignment_results.sam"))
    is_unique_full = IsUnique(os.path.join(result_path, "prbs_candidates_full_alignment_results.sam"), gene_name, num_prbs*2)
    bad_unique_full = np.zeros_like(bad_unique)
    for i in range(num_prbs):
        if is_unique_full[2*i] == 0 | is_unique_full[2*i+1] == 0:
            bad_unique_full[i] = 1

    ## Find bad probes
    bad_inds = bad_gc + bad_repeats + bad_dg + bad_unique + bad_unique_full

    # Select probe sequences that are apart
    prb_pos = np.argwhere(np.logical_not(bad_inds))
    prb_pos = np.array(prb_pos.ravel())   
    # prb_pos = [ind for sublist in prb_pos for ind in sublist]
    # print(prb_pos)
    # print(len(prb_pos))

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
        resultdf.to_excel(excel_writer = os.path.join(result_path, "probes.xlsx"))

    return resultdf


def designuseqFISHProbes(gene_id="", gene_name="", gene_host="", email=None,
                sequence="", db=os.getcwd(), barcode_path=os.getcwd(), 
                barcode_num=44, barcode="", result_path=os.getcwd(), 
                prb_length=20, gc_range=[40, 60], primer_end="TAATGTTATCTT",
                padlock_start="ACATTA", padlock_end="AAGATA", spacer1="attta",
                spacer2 = "atta", prb_space=1, dg_thresh=-9, 
                to_excel=False):
    """ Creates an Excel file in result_path containing probe designs for a given gene.
    # TODO: fill in descriptions of arguments
    Args:
        gene_id:
        gene_name:
        gene_host:
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

    if not barcode:
        barcode_df = pd.read_excel(barcode_path, index_col=0)
        barcode_db = barcode_df['barcode'].values.tolist()
        barcode = barcode_db[barcode_num-1]

    # Get Sequence
    if not sequence:
        # retrieve target sequence from genbank by using accession id
        print("Looking up gene with accession id")
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype = "gb", retmode = "text")
        target = SeqIO.read(handle, "genbank")
        if gene_name not in target.description.lower():
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
    prbs = findAllCandidates(target, prb_length, result_path)
    num_prbs = len(prbs)

    # Bowtie2 alignment
    ProbeBowtie2(os.path.join(result_path, "prbs_candidates.fasta"), db=db, result_path=os.path.join(result_path, "prbs_candidates_alignment_result.sam"), score_min="G,10,4")

    # parse sam file to get mapq and
    # find only unique probe sequences
    print(" 0. aligning probe sequences on refseq database using bowtie2")
    bad_unique = IsUnique(os.path.join(result_path, "prbs_candidates_alignment_result.sam"), gene_name, num_prbs)

    # basic filtering
    bad_gc, bad_repeats, bad_dg = basicFilter(prbs, num_prbs, prb_length=prb_length, gc_range=gc_range, dg_thresh=dg_thresh)

    # Get full probe sequences (primer and padlock) and align
    primers = []
    padlocks = []
    for i in range(num_prbs):
        primer_rec = SeqRecord(prbs[i].seq[0:prb_length] + primer_end, '%i' % (i+1), '', '')
        padlock_rec = SeqRecord(padlock_start + prbs[i].seq[prb_length:prb_length*2] \
            + spacer1 + barcode + spacer2 + padlock_end, '%i' % (i+1), '', '')
        primers.append(primer_rec)
        padlocks.append(padlock_rec)
    count = SeqIO.write(primers, os.path.join(result_path, "primers.fasta"), "fasta")
    count = SeqIO.write(padlocks, os.path.join(result_path, "padlocks.fasta"), "fasta")
    print("Converted %i records" % count)

    ProbeBowtie2(os.path.join(result_path, "primers.fasta"), db=db, result_path=os.path.join(result_path, "primers_candidates_alignment_results.sam"), score_min='G,20,8')
    ProbeBowtie2(os.path.join(result_path, "padlocks.fasta"), db=db, result_path=os.path.join(result_path, "padlocks_candidates_alignment_results.sam"), score_min='G,20,8')
    bad_unique_primers = IsUnique(os.path.join(result_path, "primers_candidates_alignment_results.sam"), gene_name, num_prbs)
    bad_unique_padlocks = IsUnique(os.path.join(result_path, "padlocks_candidates_alignment_results.sam"), gene_name, num_prbs)
    bad_unique_full = bad_unique_primers | bad_unique_padlocks

    # Find bad probes
    bad_inds = bad_gc + bad_repeats + bad_dg + bad_unique + bad_unique_full

    # Select probe sequences that are apart
    prb_pos = np.argwhere(np.logical_not(bad_inds))
    prb_pos = np.array(prb_pos.ravel())   

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
    primers_final = []
    padlocks_final = []
    for pos in prb_final_pos:
        primers_final.append(primers[pos].seq)
        padlocks_final.append('/5Phos/'+padlocks[pos].seq)

    result = {'name': gene_name, \
        'accession': gene_id, \
        'position':prb_final_pos, \
        'primer':primers_final, \
        'padlock':padlocks_final}
    resultdf = pd.DataFrame(result)
    if to_excel:
        resultdf.to_excel(excel_writer = os.path.join(result_path, "probes.xlsx"))
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

def IsUnique(samfile_path, gene_name, num_prbs=0):
    # find hits (ids) from alignment 
    hits = [[] for _ in range(num_prbs)]
    # print(hits)
    with open(samfile_path, "rb") as samfile:
        for line in samfile:
            info = line.decode().split('\t')
            hits[int(info[0])-1].append(info[2])
    
    # search hitted ids if it is variants of the same gene
    variants = ['*']
    bad_unique = np.zeros((num_prbs,), dtype=bool)
    for i, hits_for_oneprb in enumerate(hits):
        for hit_id in hits_for_oneprb:
            if hit_id not in variants:
                handle = Entrez.efetch(db="nucleotide", id=hit_id, rettype="gb", retmode="text")
                search_result = SeqIO.read(handle, "genbank")
                if gene_name in search_result.description.lower():
                    variants.append(hit_id)
                else:
                    bad_unique[i] = 1

    return bad_unique


def findAllCandidates(target, prb_length, result_path):
    """ finds all candidate probes
    """
    print("- finding all potential probes ...")
    prbs = []
    limit = len(target.seq)
    for i in range(0, limit-prb_length*2+1):
        start = i
        end = start + prb_length*2
        temp = target.seq[start:end].reverse_complement()
        temp_rec = SeqRecord(temp, '%i' % (i+1), '', '')
        prbs.append(temp_rec)
    count = SeqIO.write(prbs, os.path.join(result_path, "prbs_candidates.fasta"), "fasta")
    print("Converted %i records" % count)

    return prbs


def basicFilter(prbs, num_prbs, prb_length=20, gc_range=[40,60], dg_thresh=-9):
    print(" 1. filtering GC contents, Tm, repeats, dG ...") 
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

        dg[1,i] = calcHairpin(str(prbs[i].seq[0:prb_length])).dg/1000
        dg[0,i] = calcHairpin(str(prbs[i].seq[prb_length:prb_length*2])).dg/1000

    bad_gc = (GC < gc_range[0]) | (GC > gc_range[1])
    bad_gc = bad_gc[0,:] | bad_gc[1,:]
    bad_repeats = repeats > 0
    bad_repeats = bad_repeats[0,:] | bad_repeats[1,:]
    bad_dg = dg <= dg_thresh
    bad_dg = bad_dg[0,:] | bad_dg[1,:]
    
    return bad_gc, bad_repeats, bad_dg
