import probe_design
import os

RESULT_PATH = os.path.join(os.getcwd(), "probe_design_files")
if not os.path.isdir(RESULT_PATH):
    os.mkdir(RESULT_PATH)

## input
email = "jvendemiatti@hmc.edu"
gene_id = "NM_009891.2"
gene_name = "chat"
gene_host = "mus musculus"
sequence = ""

# if retrieve target sequence from genbank make this empty
# gene_name = "mneongreen"
# sequence = "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCTCTCTCCCAGCGACACATGAGTTACACATCTTTGGCTCCATCAACGGTGTGGACTTTGACATGGTGGGTCAGGGCACCGGCAATCCAAATGATGGTTATGAGGAGTTAAACCTGAAGTCCACCAAGGGTGACCTCCAGTTCTCCCCCTGGATTCTGGTCCCTCATATCGGGTATGGCTTCCATCAGTACCTGCCCTACCCTGACGGGATGTCGCCTTTCCAGGCCGCCATGGTAGATGGCTCCGGATACCAAGTCCATCGCACAATGCAGTTTGAAGATGGTGCCTCCCTTACTGTTAACTACCGCTACACCTACGAGGGAAGCCACATCAAAGGAGAGGCCCAGGTGAAGGGGACTGGTTTCCCTGCTGACGGTCCTGTGATGACCAACTCGCTGACCGCTGCGGACTGGTGCAGGTCGAAGAAGACTTACCCCAACGACAAAACCATCATCAGTACCTTTAAGTGGAGTTACACCACTGGAAATGGCAAGCGCTACCGGAGCACTGCGCGGACCACCTACACCTTTGCCAAGCCAATGGCGGCTAACTATCTGAAGAACCAGCCGATGTACGTGTTCCGTAAGACGGAGCTCAAGCACTCCAAGACCGAGCTCAACTTCAAGGAGTGGCAAAAGGCCTTTACCGATGTGATGGGCATGGACGAGCTGTACAAGTGAG"


barcode_num = 44

## parameters
prb_length = 20
gc_range = [40, 60]
prb_space = 1
dg_thresh = -9

primer_end = "TAATGTTATCTT"
padlock_start = "ACATTA"
padlock_end = "AAGATA"
spacer1 = "attta"
spacer2 = "atta"

resultdf = probe_design.designuseqFISHProbes(gene_id=gene_id, 
                                gene_name=gene_name, 
                                gene_host=gene_host,
                                email=email,
                                sequence=sequence,
                                db="db/mouse/mouse_refseq_rna",
                                barcode_path="db/mouse/barcodes/barcodes.xlsx",
                                barcode_num=44,
                                result_path=RESULT_PATH,
                                prb_length=prb_length,
                                gc_range=gc_range,
                                primer_end=primer_end,
                                padlock_start=padlock_start,
                                padlock_end=padlock_end,
                                spacer1=spacer1,
                                spacer2=spacer2,
                                prb_space=prb_space,
                                dg_thresh=dg_thresh,
                                to_excel=True)