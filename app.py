import re
import probe_design

from flask import Flask, render_template, url_for, request, redirect, flash, session

import os
import traceback


app = Flask(__name__)


app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///test.db'

RESULT_PATH = os.path.join(os.getcwd(), "probe_design_files")
if not os.path.isdir(RESULT_PATH):
    os.mkdir(RESULT_PATH)
app.config['RESULT_PATH'] = RESULT_PATH
app.secret_key = "glab"

@app.route('/', methods=['GET','POST'])
def index():
    if request.method == 'GET':
        return render_template("index.html")
    else:
        prb_type = request.form['prb_type']
        if prb_type == 'useqFISH':
            return render_template('useqFISH.html')
        elif prb_type == 'HCR3':
            return render_template('HCR3.html')
        else:
            return 'there was an error redirecting to correct prb_type page'


@app.route('/useqFISH', methods=['POST', 'GET'])
def useqFISH():
    if request.method == 'POST':
        species = request.form.get("species")
        gene_name = request.form.get('gene_name')
        gene_id = request.form.get('gene_id')
        sequence = request.form.get('sequence')
        barcode = request.form.get('barcode')
        gene_host = request.form.get('gene_host')
        email = request.form.get('email')

        primer_end = request.form.get('primer_end')
        padlock_start = request.form.get('padlock_start')
        padlock_end = request.form.get('padlock_end')
        prb_length = request.form.get('prb_length', type=int)
        gc_range = [request.form.get('min_gc', type=int), request.form.get('max_gc', type=int)]
        prb_space = request.form.get('prb_space', type=int)
        dg_thresh = request.form.get('dg_thresh', type=int)

        db = getdb(species)

        try:

            resultdf = probe_design.designuseqFISHProbes(gene_id=gene_id, 
                        gene_name=gene_name, 
                        gene_host=gene_host, 
                        email=email, 
                        sequence=sequence,
                        primer_end=primer_end,
                        padlock_start=padlock_start,
                        padlock_end=padlock_end,
                        barcode=barcode,
                        db=db, 
                        result_path=app.config['RESULT_PATH'],
                        prb_length=prb_length,
                        gc_range=gc_range,
                        prb_space=prb_space,
                        dg_thresh=dg_thresh)
    
           
            return render_template('useqFISH.html', result=resultdf.to_html())

        except Exception as e: 
            print(traceback.format_exc())
            return 'There was an issue running the Probe Designer.'

    else:
        return render_template("useqFISH.html")



@app.route('/HCR3', methods=['POST', 'GET'])
def HCR3():
    if request.method == 'POST':
        species = request.form.get("species")
        gene_name = request.form.get('gene_name')
        gene_id = request.form.get('gene_id')
        hairpin_id = request.form.get('hairpin_id', type=int)
        email = request.form.get('email')

        prb_length = request.form.get('prb_length', type=int)
        gc_range = [request.form.get('min_gc', type=int), request.form.get('max_gc', type=int)]
        prb_space = request.form.get('prb_space', type=int)
        dg_thresh = request.form.get('dg_thresh', type=int)

        db = getdb(species)

        try:

            resultdf = probe_design.designHCR3Probes(gene_id=gene_id, 
                        gene_name=gene_name, 
                        hairpin_id=hairpin_id, 
                        email=email, 
                        db=db, 
                        result_path=app.config['RESULT_PATH'],
                        prb_length=prb_length,
                        gc_range=gc_range,
                        prb_space=prb_space,
                        dg_thresh=dg_thresh)

            return render_template('HCR3.html', result=resultdf.to_html())

        except Exception as e: 
            print(traceback.format_exc())
            return 'There was an issue running the Probe Designer.'

    else:
        return render_template("HCR3.html")


def getdb(species):
    if species == "mouse":
        return os.path.join(os.getcwd(), "db/mouse/mouse_refseq_rna")
    elif species =="macaque":
        return os.path.join(os.getcwd(), "db/macaque/macaque_refseq_rna")
    else:
        raise ValueError("No dataase for selected species")


if __name__ == "__main__":
    app.run(debug=True)