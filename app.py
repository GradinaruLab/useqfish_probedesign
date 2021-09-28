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
app.secret_key = os.urandom(24)


@app.route('/', methods=['POST', 'GET'])
def index():
    if request.method == 'POST':
        gene_name = request.form['gene_name']
        gene_id = request.form['gene_id']
        hairpin_id = int(request.form['hairpin_id'])
        email = request.form['email']

        try:
            if 'prb_length' in session:
                prb_length = session['prb_length']
                gc_range = session['gc_range']
                hits_thresh = session['hits_thresh']
                num_offtar_thresh = session['num_offtar_thresh']
                prb_space= session['prb_space']
                dg_thresh = session['dg_thresh']

            else:
                prb_length = 20
                gc_range = [40, 60]
                hits_thresh = 5
                num_offtar_thresh = 0
                prb_space = 1
                dg_thresh = -9

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
            session['prb_length'] = request.form['prb_length']
            session['gc_range'] = [request.form['min_gc'], request.form['max_gc']]
            session['hits_thresh'] = request.form['hits_thresh']
            session['num_offtar_thresh'] = request.form['num_offtar_thresh']
            session['prb_space'] = request.form['prb_space']
            session['dg_thresh'] = request.form['dg_thresh']
            return redirect('/')
        
        except:
            return 'There was an issue updating the parameters'
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