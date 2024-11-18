from fileinput import filename 
import FQ_VSearch.find_variable_sites as FVS

import os

from flask import *
from pyngrok import ngrok
from flask_cors import CORS
app = Flask(__name__) 
CORS(app)

@app.route('/') 
def main(): 
	return render_template("index.html") 

@app.route('/success', methods = ['POST']) 
def success(): 
	if request.method == 'POST': 
		fasta_file = request.files['file1']
		fastq_file = request.files['file2'] 
		out_path = request.form.get('file_name')
		phread_score = int(request.form.get('phred_score'))
		distance_5_prime = int(request.form.get('distance_5_prime'))
		distance_3_prime = int(request.form.get('nucleotide_match'))
		
		fasta_file.save(fasta_file.filename) 
		fastq_file.save(fastq_file.filename) 

		variable_sites_number = int(request.form.get('variable_sites_number'))
	
		out_path = FVS.main(fasta_file.filename,fastq_file.filename,out_path,True,True,phread_score,distance_5_prime,distance_3_prime,variable_sites_number)
		
		return send_file(out_path, as_attachment=True)



#public_url = ngrok.connect(5000)
#print(f"public_url:\t{public_url}")

if __name__ == '__main__': 
	app.run(port=5000,debug=True)
