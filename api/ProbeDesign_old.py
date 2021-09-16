
# import ray
import time
import pandas as pd

from Bio.Blast.Applications import  NcbiblastnCommandline
from Bio.Blast import NCBIXML, Record
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_nucleotide
from Bio import SeqIO, Entrez, SearchIO, SeqFeature
from Bio.SeqUtils import MeltingTemp as mt

from primer3 import calcHairpin     #dependency: primer3-py

import numpy as np
import csv
import re
import matplotlib.pyplot as plt

# ray.init()


## input
Entrez.email = "mjjang@caltech.edu"
gene_id = "NM_182993.2"
gene_name = "slc17a7"
gene_synonym = []
barcode_num = 31
# barcode = "TTACTATCTCAATAAATAT"


## parameters
prb_length = 20             # default: 20
gc_range = [40, 60]         # [40, 60]
tm_diff = 1                 # 1
hits_thresh = 5             # 5
num_offtar_thresh = 0       # 0
prb_space = 1               # 1
dg_thresh = -9 # kcal/mol   # 9

primer_end = "TAATGTTATCTT"
padlock_start = "ACATTA"
padlock_end = "AAGATA"
spacer1 = "ATTTA"
spacer2 = "ATTA"

barcode_df = pd.read_excel("barcodes.xlsx", index_col=0)
barcode_db = barcode_df['barcode'].values.tolist()
barcode = barcode_db[barcode_num-1]


def ProbeBlast(fastafile):
    cmd = NcbiblastnCommandline(query=fastafile, db="mouse/mouse_transcript", \
        outfmt=5, evalue=1000, max_target_seqs=5000, \
            task="blastn", \
                out=fastafile+"_blast_results.xml", \
                    num_threads="3")
    stdout, stderr = cmd()

    result_handle = open(fastafile+"_blast_results.xml")
    blast_results = NCBIXML.parse(result_handle)
    # blast_result = next(blast_results)
    # print(blast_result.descriptions[1].title)

    return blast_results



## retrieve target sequence from genbank by using accession id
handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype = "gb", retmode = "text")
#print(handle.read())
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
prb_list = []
limit = len(target.seq)
for i in range(0, limit-prb_length*2+1):
    start = i
    end = start + prb_length*2
    temp = target.seq[start:end]
    temp_rec = SeqRecord(temp, '%i' % (i+1), '', '')
    prb_list.append(temp_rec)

num_prbs = len(prb_list)
SeqIO.write(prb_list, "prblist.fasta", "fasta")

# filter GC contents, tm, repeats
print(" 0. filtering GC contents, Tm, repeats, dG ...") 
GC = np.zeros((2, num_prbs))
tm  = np.zeros((2, num_prbs))
repeats = np.zeros((2, num_prbs))
dg = np.zeros((2, num_prbs))


for i in range(0, num_prbs):
    GC[0,i] = (prb_list[i].seq[0:prb_length].count("G") \
        + prb_list[i].seq[0:prb_length].count("C"))/prb_length *100
    GC[1,i] = (prb_list[i].seq[prb_length:prb_length*2].count("G") \
        + prb_list[i].seq[prb_length:prb_length*2].count("C"))/prb_length *100

    tm[0,i] = mt.Tm_NN(prb_list[i].seq[0:prb_length], nn_table=mt.DNA_NN1)
    tm[1,i] = mt.Tm_NN(prb_list[i].seq[prb_length:prb_length*2], nn_table=mt.DNA_NN1)

    repeats[0,i] = prb_list[i].seq[0:prb_length].count("AAAA") \
        + prb_list[i].seq[0:prb_length].count("CCCC") \
            + prb_list[i].seq[0:prb_length].count("GGGG") \
                + prb_list[i].seq[0:prb_length].count("TTTT")
    repeats[1,i] = prb_list[i].seq[prb_length:prb_length*2].count("AAAA") \
        + prb_list[i].seq[prb_length:prb_length*2].count("CCCC") \
            + prb_list[i].seq[prb_length:prb_length*2].count("GGGG") \
                + prb_list[i].seq[prb_length:prb_length*2].count("TTTT")

    dg[0,i] = calcHairpin(str(prb_list[i].seq[0:prb_length])).dg/1000
    dg[1,i] = calcHairpin(str(prb_list[i].seq[prb_length:prb_length*2])).dg/1000

bad_gc = (GC < gc_range[0]) | (GC > gc_range[1])
bad_gc = bad_gc[0,:] | bad_gc[1,:]
bad_tm = np.absolute(tm[0,:]-tm[1,:]) > tm_diff
bad_repeats = repeats > 0
bad_repeats = bad_repeats[0,:] | bad_repeats[1,:]
bad_dg = dg <= dg_thresh
bad_dg = bad_dg[0,:] | bad_dg[1,:]


## blast target sequences
print(" 1. blasting target sequences ...")
t = time.time()
blast_results = ProbeBlast("prblist.fasta")

hits_name_list = []
hits_pos_list = []
bad_100hits = np.zeros(num_prbs, dtype=bool)
num_hits = np.zeros(num_prbs)
for i, blast_result in enumerate(blast_results):
    hits_name_list_fori = []
    hits_pos_list_fori = []
    for alignment in blast_result.alignments:
        isoff = (alignment.title.lower().find('predicted') == -1) & \
            (alignment.title.lower().find(gene_name) == -1)
        if len(gene_synonym) > 0:
            for name in gene_synonym:
                isoff = isoff & (alignment.title.lower().find(name) == -1)
        if isoff:
            for hsp in alignment.hsps:
                hits_name = re.split(r'\(|\)', alignment.title.lower())[-2]
                qcovs = hsp.align_length / blast_result.query_length * 100
                pident = hsp.identities/ hsp.align_length * 100
                
                if hits_name not in hits_name_list_fori:
                    hits_pos_list_fori.append(i)
                    hits_name_list_fori.append(hits_name)
                    
                if (qcovs == 100) & (pident == 100):
                    # print(hits_name)
                    bad_100hits[i] = True
    hits_name_list.append(hits_name_list_fori)
    hits_pos_list.append(hits_pos_list_fori)
    num_hits[i] = len(hits_name_list_fori)


print('    Elapsed time: %.2f sec' % (time.time() - t))

# np.asarray(bad_100hits, dtype=bool)
# np.asarray(num_hits)

## create full probe sequences
primer_list = []
padlock_list = []
for i in range(0, num_prbs):
    seqrc = prb_list[i].seq.reverse_complement()
    primer_rec = SeqRecord(seqrc[0:prb_length] + primer_end, '%i' % (i+1), '', '')
    padlock_rec = SeqRecord(padlock_start + seqrc[prb_length:prb_length*2] \
        + spacer1 + barcode + spacer2 + padlock_end, '%i' % (i+1), '', '')
    primer_list.append(primer_rec)
    padlock_list.append(padlock_rec)
    # print(prb_list[i])
    # print(seqrc)
    # print(primer_list)
    # print(padlock_list)

SeqIO.write(primer_list, "primerlist.fasta", "fasta")
SeqIO.write(padlock_list, "padlocklist.fasta", "fasta")

## blast full sequences - if hit other than the target, the full seq dropped
print(" 2. blasting full sequences ...")
t = time.time()
prm_blast_results = ProbeBlast("primerlist.fasta")
pdl_blast_results = ProbeBlast("padlocklist.fasta")

bad_fullhits = np.zeros(num_prbs, dtype=bool)
for i, (prm_blast_result, pdl_blast_result) in enumerate(zip(prm_blast_results, pdl_blast_results)):
    for prm_alignment, pdl_alignment in zip(prm_blast_result.alignments, pdl_blast_result.alignments):
        prm_isoff = (prm_alignment.title.lower().find('predicted') == -1) & \
            (prm_alignment.title.lower().find(gene_name) == -1)
        pdl_isoff = (pdl_alignment.title.lower().find('predicted') == -1) & \
            (pdl_alignment.title.lower().find(gene_name) == -1)
        if len(gene_synonym) > 0:
            for name in gene_synonym:
                prm_isoff = prm_isoff & (prm_alignment.title.lower().find(name) == -1)
                pdl_isoff = pdl_isoff & (pdl_alignment.title.lower().find(name) == -1)
        if prm_isoff | pdl_isoff:
            # print(prm_alignment)
            # print(prm_alignment.hsps[0].query[0:32])
            # print(prm_alignment.hsps[0].match[0:32])
            # print(prm_alignment.hsps[0].sbjct[0:32])
            # print(pdl_alignment)
            bad_fullhits[i] = True
print('    Elapsed time: %.2f sec' % (time.time() - t))


## final decision
bad_tarhits = num_hits > hits_thresh

# bad_inds = np.array([bad_gc, bad_repeats, bad_tarhits, bad_100hits, bad_fullhits])
bad_inds = bad_gc + bad_tm + bad_repeats + bad_dg + bad_tarhits + bad_100hits + bad_fullhits

howmany = np.array([bad_gc, bad_tm, bad_repeats, bad_dg, bad_tarhits, bad_100hits, bad_fullhits])
howmany = howmany*1
df = pd.DataFrame(howmany).T
df.to_excel(excel_writer = "howmany.xlsx")
# howmany = np.sum(howmany, axis=1)
# print(howmany)

max_num_offtar = num_prbs
while max_num_offtar > num_offtar_thresh:
    prb_pos = np.argwhere(np.logical_not(bad_inds))
    prb_pos = np.array(prb_pos.ravel())    
    print('# of probes passed: %i' % prb_pos.size)
    if prb_pos.size == 0:
        print('==> No probe passed filters!')
        exit()

    # collecting the name of off-targets
    offtar_hits_names = []
    for i in range(prb_pos.size):
        offtar_hits_names.extend([hits_name_list[prb_pos[i]]])
    offtar_hits_names = [name for sublist in offtar_hits_names for name in sublist]

    # find the off-target hitted by the max number of probes
    offtar_hits_names_uniq = []
    offtar_hits_count = []

    for name in offtar_hits_names:
        if len(name) > 0:
            if name not in offtar_hits_names_uniq:
                offtar_hits_names_uniq.append(name)
                offtar_hits_count.append(1)    
            else:
                hi = offtar_hits_names_uniq.index(name)
                offtar_hits_count[hi] += 1
    
    offtar_hits_count = np.array(offtar_hits_count)
    print(offtar_hits_names_uniq)
    print(offtar_hits_count)

    if offtar_hits_count.size == 0:
        break 

    max_num_offtar = np.amax(offtar_hits_count)
    max_idx = int(np.argwhere(offtar_hits_count == max_num_offtar)[0])
    max_offtar_hits_name = offtar_hits_names_uniq[max_idx]
    
    # print(max_num_offtar)
    # print(max_offtar_hits_name)
    
    # drop probes that hit the off-target above
    for hi, hname in enumerate(hits_name_list):
        if max_offtar_hits_name in hname:
            # print(hname)
            # print(max_offtar_hits_name in hname)
            bad_inds[hi] = True

# Select probe sequences that are apart
prb_pos = np.argwhere(np.logical_not(bad_inds))
prb_pos = np.array(prb_pos.ravel())   
# prb_pos = [ind for sublist in prb_pos for ind in sublist]
print(prb_pos)
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
# print(prb_final_pos)
print('- done. # of final probes: %i' % len(prb_final_pos))

# recall primer and padlock sequences for each position
prb_final_primer = []
prb_final_padlock = []
for pos in prb_final_pos:
    prb_final_primer.append(str(primer_list[pos].seq))
    prb_final_padlock.append("/5Phos/"+str(padlock_list[pos].seq))

result = {'name': gene_name, \
    'accession': gene_id, \
    'position':prb_final_pos, \
    'primer':prb_final_primer.seq, \
    'padlock':prb_final_padlock.seq}
resultdf = pd.DataFrame(result)
resultdf.to_excel(excel_writer = "probes.xlsx")
# print(prb_final_primer)

fig, ax = plt.subplots()
ax.plot(np.linspace(0, num_hits.size, num_hits.size), num_hits)
ax.grid()
plt.show()