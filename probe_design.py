import subprocess
import os

from Bio import Align, SeqIO, Entrez, SearchIO, SeqFeature
from Bio.SeqRecord import SeqRecord

from primer3 import calcHairpin     #dependency: primer3-py

import numpy as np
import pandas as pd
from sqlalchemy.orm.query import AliasOption


def designProbes(gene_id="", gene_name="", hairpin_id=None, email=None,
                db=os.getcwd(), result_path=os.getcwd(), 
                prb_length=20, gc_range=[40, 60], 
                hits_thresh=5, num_offtar_thresh=0, prb_space=1, dg_thresh=-9, 
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

    ## probe design
    # create a fasta file including all candidates
    print("- finding all potential probes ...")
    prbs = []
    limit = len(target.seq)
    for i in range(0, limit-prb_length*2+1):
        start = i
        end = start + prb_length*2
        temp = target.seq[start:end].reverse_complement()
        temp_rec = SeqRecord(temp, '%i' % (i+1), '', '')
        prbs.append(temp_rec)

    num_prbs = len(prbs)
    count = SeqIO.write(prbs, os.path.join(result_path, "prbs_candidates.fasta"), "fasta")
    print("Converted %i records" % count)

    ## bowtie2 alignment
    ProbeBowtie2(os.path.join(result_path, "prbs_candidates.fasta"),db=db, result_path=os.path.join(result_path, "alignment_results.sam"))

    ## parse sam file to get mapq and
    ## find only unique probe sequences
    print(" 0. aligning probe sequences on refseq database using bowtie2")
    is_unique = IsUnique(os.path.join(result_path, "alignment_results.sam"))
    bad_unique = np.logical_not(is_unique)


    ## basic filtering
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
    is_unique_full = IsUnique(os.path.join(result_path, "prbs_candidates_full_alignment_results.sam"))
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


##### HELPER FUNCTIONS #####


def ProbeBowtie2(fastafile, db=os.path.join(os.getcwd(), 'mouse_refseq_rna'), result_path="prbs_candidates_alignment_result.sam"):
    subprocess.check_call(['bowtie2', '--very-sensitive-local', '-f', '--no-sq', '--no-hd', '--reorder', '--score-min', 'G,10,4', \
        '-x', db, '-U', fastafile, '-S', result_path])
    return


def GetInitiatorSeq(hairpin_id=2, I_id=2):
    init_seqs = [['gAggAgggCAgCAAACgggAAgAgTCTTCCTTTACg', 'gCATTCTTTCTTgAggAgggCAgCAAACgggAAgAg'], \
        ['CCTCgTAAATCCTCATCAATCATCCAgTAAACCgCC', 'AgCTCAgTCCATCCTCgTAAATCCTCATCAATCATC'], \
            ['gTCCCTgCCTCTATATCTCCACTCAACTTTAACCCg', 'AAAgTCTAATCCgTCCCTgCCTCTATATCTCCACTC'], \
                ['CCTCAACCTACCTCCAACTCTCACCATATTCgCTTC', 'CACATTTACAgACCTCAACCTACCTCCAACTCTCAC'], \
                    ['CTCACTCCCAATCTCTATCTACCCTACAAATCCAAT', 'CACTTCATATCACTCACTCCCAATCTCTATCTACCC']]
    return init_seqs[hairpin_id-1][I_id-1]

def IsUnique(samfile_path):
    is_unique = []
    # with open("probe_design/prbs_alignment_result.sam", 'rb') as samfile:
    with open(samfile_path, "rb") as samfile:
        for line in samfile:
            if line.decode().find('XS') > -1:
                is_unique.append(0)
            else:
                is_unique.append(1)
    return np.array(is_unique, dtype=bool)