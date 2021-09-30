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


@app.route('/', methods=['POST', 'GET'])
def index():
    if request.method == 'POST':
        gene_name = request.form['gene_name']
        gene_id = request.form['gene_id']
        hairpin_id = int(request.form['hairpin_id'])
        email = request.form['email']

        session['gene_id'] = gene_id
        session['gene_name'] = gene_name
        session['hairpin_id'] = hairpin_id
        session['email'] = email
        try:
            prb_length =  session['prb_length'] if 'prb_length' in session else 20
            gc_range = session['gc_range'] if 'gc_range' in session else [40, 60]
            hits_thresh = session['hits_thresh'] if 'hits_thresh' in session else 5
            num_offtar_thresh = session['num_offtar_thresh'] if 'num_offtar_thresh' in session else 0
            prb_space= session['prb_space'] if 'probe_space' in session else 1
            dg_thresh = session['dg_thresh'] if 'dg_thresh' in session else -9

            resultdf = probe_design.designProbes(gene_id=gene_id, gene_name=gene_name, 
                        hairpin_id=hairpin_id, email=email, 
                        db=os.path.join(os.getcwd(), "db/mouse/mouse_refseq_rna"), 
                        result_path=app.config['RESULT_PATH'],
                        prb_length=prb_length,
                        gc_range=gc_range,
                        hits_thresh=hits_thresh,
                        num_offtar_thresh=num_offtar_thresh,
                        prb_space=prb_space,
                        dg_thresh=dg_thresh)
            
            session['result'] = resultdf.to_html()
           
            return redirect('/results')

        except Exception as e: 
            print(traceback.format_exc())
            return 'There was an issue running the Probe Designer.'

    else:
        return render_template("index.html")


@app.route('/parameters', methods=['GET', 'POST'])
def update_parameters():
    if request.method == 'POST':

        try:
            session['prb_length'] = request.form.get('prb_length', type=int)
            session['gc_range'] = [request.form.get('min_gc', type=int), request.form.get('max_gc', type=int)]
            session['hits_thresh'] = request.form.get('hits_thresh', type=int)
            session['num_offtar_thresh'] = request.form.get('num_offtar_thresh', type=int)
            session['prb_space'] = request.form.get('prb_space', type=int)
            session['dg_thresh'] = request.form.get('dg_thresh', type=int)
            return redirect('/')
        
        except:
            return 'There was an issue updating the parameters'
            # return f'type of prb_length is {type(session["prb_legnth"])}'
    else:
        return render_template('parameters.html')


@app.route('/results', methods=['GET', 'POST'])
def get_results():
    if request.method == 'GET':
        return render_template('results.html', result = session['result'] )
    else:
        return redirect('/results')

if __name__ == "__main__":
    app.run(debug=True)