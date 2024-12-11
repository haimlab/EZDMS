from fileinput import filename 
import find_variable_sites as FVS

import Re_scaling as RSC

import os

from flask import *
from pyngrok import ngrok
from flask_cors import CORS
app = Flask(__name__) 
CORS(app)

@app.route('/',methods=['GET', 'POST']) 
def main():
	return f'<html><body> <script> function change_page(){'{window.location.href = "/find_variable_sites";}'}  </script> <script> function change_page_pref(){'{window.location.href = "/preference";}'}  </script> <h1>Main page:</h1><h2>\n\n</h2></body></html> <input type="button" value="main page" onclick="change_page()"/>  <input type="button" value="preference" onclick="change_page_pref()"/>'
@app.route('/find_variable_sites',methods=['GET', 'POST']) 
def find_variable_sites(): 
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
			print(variable_sites_number)
		
			out_path = FVS.main(fasta_file.filename,fastq_file.filename,out_path,True,True,phread_score,distance_5_prime,distance_3_prime,variable_sites_number)
			try:
				return send_file(out_path, as_attachment=True)
			except:
				return redirect(url_for('Error' , error=out_path.split(":")[1]))
		except:
			return redirect(url_for('Error' , error="EEE"))

@app.route('/preference',methods=['GET', 'POST']) 
def preference():
	if request.method == 'GET':
		return render_template("preference.html")
	if request.method == 'POST': 
		try:
			try:
				file_fasta_1 = request.files['file_fasta_1']
				file_fastq_1 = request.files['file_fastq_1']
				file_fastq_2 = request.files['file_fastq_2'] 
				file_fasta_1.save(file_fasta_1.filename) 
				file_fastq_1.save(file_fastq_1.filename) 
				file_fastq_2.save(file_fastq_2.filename) 
			except:
				return redirect(url_for('Error' , error="invalid file, please check your files"))
			
			out_path = request.form.get('file_name')
			phread_score = int(request.form.get('phred_score'))


			distance_5_prime = int(request.form.get('distance_5_prime'))
			distance_3_prime = int(request.form.get('nucleotide_match'))

			print(file_fasta_1.filename)

			variable_sites = int(request.form.get('variable_sites_number'))
		
			pre_amino_dict = FVS.main(file_fasta_1.filename,file_fastq_1.filename,"",True,True,phread_score,distance_5_prime,distance_3_prime,variable_sites)
			post_amino_dict = FVS.main(file_fasta_1.filename,file_fastq_2.filename,"",True,True,phread_score,distance_5_prime,distance_3_prime,variable_sites)
			out_path = RSC.main(pre_amino_dict,post_amino_dict,out_path)
			try:
				return send_file(out_path, as_attachment=True)
			except:
				return redirect(url_for('Error' , error=out_path.split(":")[1]))
		except:
			return redirect(url_for('Error' , error=pre_amino_dict))

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