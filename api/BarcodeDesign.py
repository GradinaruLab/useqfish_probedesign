import random
import re
import numpy as np
import pandas as pd

from Bio import Align, SeqIO, pairwise2
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment

def SeqRandom(length=int(), letters="CGTA"):
    return ''.join(random.choices(letters, k=length))

def SWAlign(seq1=str(), seq2=str()):
    # a = pairwise2.align.globalxx(seq1, seq2)
    # print(format_alignment(*a[0]))
    score = pairwise2.align.globalxx(seq1, seq2, score_only=True)
    return score
    
def ProbeBlast(fastafile):
    cmd = NcbiblastnCommandline(query=fastafile, db="probe_design/mouse/mouse_transcript", \
        outfmt=5, evalue=1, \
            task="blastn", \
                out=fastafile+"_blast_results.xml")
    stdout, stderr = cmd()

    result_handle = open(fastafile+"_blast_results.xml")
    blast_results = NCBIXML.parse(result_handle)
    # blast_result = next(blast_results)
    # print(blast_result.descriptions[1].title)

    return blast_results


## parameters #############
# gc_range = [10, 20]
# tm_thresh = 40
# scr_thresh = 17

# barcode_length = 19
# barcode_num = 100

gc_range = [10, 40]
tm_thresh = 25
scr_thresh = 17
barcode_length = 17
barcode_num = 50

# gc_range = [40, 60]
# tm_thresh = 70
# # scr_thresh = 17
# scr_thresh = 3
# barcode_length = 20
# barcode_num = 30
###########################

# read barcode database (barcodes.xlsx)
barcodedf = pd.read_excel("probe_design/barcodes2.xlsx", index_col=0)
barcode_db = barcodedf['barcode'].values.tolist()
barcode_gc = barcodedf['GC'].values.tolist()
barcode_tm = barcodedf['Tm'].values.tolist()

# barcode_db = []
# barcode_gc = []
# barcode_tm = []
num = 0
while num < barcode_num:

    # generate random sequence composed of "C", "T", and "A" 
    barcode_seq = SeqRandom(barcode_length, letters="CTAG")
    # print(barcode_seq)

    # filter repeats
    repeats = barcode_seq.count("AAAA") \
        + barcode_seq.count("TTTT") \
            + barcode_seq.count("CCCC") \
                + barcode_seq.count("GGGG")
    if repeats > 0:
        continue

    # filter GC and Tm
    # gc = barcode_seq.count("C")/barcode_length*100
    gc = (barcode_seq.count("C")+barcode_seq.count("G"))/barcode_length*100
    tm = mt.Tm_NN(barcode_seq, nn_table=mt.DNA_NN1)
    
    if (gc < gc_range[0]) | (gc > gc_range[1]) | (tm > tm_thresh):
        continue

    # filter barcode that aligns too much with previous ones
    score = []
    if len(barcode_db) > 0:
        for barcode in barcode_db:
            # score.append(SWAlign(barcode_seq, barcode))
            # score.append(pairwise2.align.localms(barcode_seq, barcode, 1,-5,-5,-5, score_only=True))
            alignment = pairwise2.align.localms(barcode_seq, barcode, 1, -5, -5, -5)
            align1, align2, s, begin, end = alignment[0]
            start1 = str(len(align1[:begin]) - align1[:begin].count("-") + 1) + " "
            start2 = str(len(align2[:begin]) - align2[:begin].count("-") + 1) + " "
            # print(start1, start2)
            # print(pairwise2.format_alignment(*alignment[0]))            
            if start1 == start2:
                score.append(s)
            else:
                score.append(0)

        score_above = np.array(score) > scr_thresh
        if np.sum(score_above) > 0:
            continue
        

    # blast
    barcode_rec = SeqRecord(Seq(barcode_seq), '1', '', '')
    SeqIO.write(barcode_rec, "probe_design/barcode.fasta", "fasta")
    blast_result = ProbeBlast("probe_design/barcode.fasta")
    bresult = next(blast_result)

    offtarget_exist = False
    for alignment in bresult.alignments:
        offtarget_exist = alignment.title.lower().find('predicted)') == -1
    
    if offtarget_exist:
        continue
    

    barcode_gc.append(gc)
    barcode_tm.append(tm)
    barcode_db.append(barcode_seq)
    num += 1
    
# print barcodes and readers generated
result = {'barcode':barcode_db, \
    'GC':barcode_gc, \
        'Tm':barcode_tm}
resultdf = pd.DataFrame(result)
resultdf.to_excel(excel_writer = "probe_design/barcodes2.xlsx")
    
    