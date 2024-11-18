from fileinput import filename 
import find_variable_sites as FVS

import os

from flask import *
from pyngrok import ngrok
from flask_cors import CORS
app = Flask(__name__) 
CORS(app)

@app.route('/',methods=['GET', 'POST']) 
def main(): 
	if request.method == 'GET':
		return render_template("index.html")
	if request.method == 'POST': 
		try:
			try:
				fasta_file = request.files['file1']
				fastq_file = request.files['file2'] 
				fasta_file.save(fasta_file.filename) 
				fastq_file.save(fastq_file.filename) 
			except:
				return redirect(url_for('Error' , error="invalid file, please check your files"))

			out_path = request.form.get('file_name')
			phread_score = int(request.form.get('phred_score'))
			distance_5_prime = int(request.form.get('distance_5_prime'))
			distance_3_prime = int(request.form.get('nucleotide_match'))
			


			variable_sites_number = int(request.form.get('variable_sites_number'))
		
			out_path = FVS.main(fasta_file.filename,fastq_file.filename,out_path,True,True,phread_score,distance_5_prime,distance_3_prime,variable_sites_number)
			try:
				return send_file(out_path, as_attachment=True)
			except:
				return redirect(url_for('Error' , error=out_path.split(":")[1]))
		except:
			return redirect(url_for('Error' , error="EEE"))
		
@app.route('/Error/<error>') 
def Error(error): 
	return f'<html><body> <script> function change_page(){'{window.location.href = "/";}'} </script><h1>An Error has been detected:</h1><h2>{error}\n\n</h2></body></html> <input type="button" value="main page" onclick="change_page()"/>'



"""
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
"""


#public_url = ngrok.connect(5000)
#print(f"public_url:\t{public_url}")

if __name__ == '__main__': 
	app.run(port=5000,debug=True)