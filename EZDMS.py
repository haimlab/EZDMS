from fileinput import filename 

import webbrowser

import threading

import os

import time

import sys

from pathlib import Path


#imports

import webview

from flask import *

from flask_cors import CORS

import requests


sys.path.insert(0, '/modules')

import modules.find_variable_sites as FVS
import modules.BuildPrimer as BP
import modules.Re_scaling as RSC
import modules.customErrors as CE

app = Flask(__name__, template_folder='./templates') 

# Function to open the link in the default system browser
class JSAPI:
	def open_link(self, url):
		webbrowser.open(url)  # Opens the URL in the default browser

def start_flask():
	
	CORS(app)

	# Determine where to save files (either current directory or temp directory in executable)
	if getattr(sys, 'frozen', False):
		# If running from a PyInstaller packaged executable
		base_path = sys._MEIPASS  # Temporary directory where PyInstaller unpacks resources
	else:
		# If running in development mode
		base_path = os.path.dirname(__file__)  # Current working directory

	upload_folder = os.path.join(base_path, 'uploads')
	if not os.path.exists(upload_folder):
		os.makedirs(upload_folder)  # Ensure the 'uploads' folder exists

	app.config['UPLOAD_FOLDER'] = upload_folder

	@app.route('/',methods=['GET', 'POST']) 
	def home():
		files = os.listdir(upload_folder)
		for file in files:
			file_path = os.path.join(upload_folder, file)
			if os.path.isfile(file_path):
				os.remove(file_path)
		return render_template("home_page.html")
	
	@app.route('/css/<path:filename>')
	def send_css(filename):
		return send_from_directory('templates/css', filename)
	

	@app.route('/EGG',methods=['GET', 'POST']) 
	def test():
		return render_template("EGG.html")

	@app.route('/find_variable_sites',methods=['GET', 'POST']) 
	def find_variable_sites(): 
		out_path = None
		print(out_path)
		if request.method == 'GET':
			return render_template("FQ_VSearch.html")
		if request.method == 'POST': 
				try:
					fasta_file = request.files['file1']
					fastq_file = request.files['file2'] 
					fasta_file.save(os.path.join(app.config['UPLOAD_FOLDER'],fasta_file.filename))
					fastq_file.save(os.path.join(app.config['UPLOAD_FOLDER'],fastq_file.filename))

					out_path_name = os.path.join(app.config['UPLOAD_FOLDER'],request.form.get('file_name'))
					in_fastq_file = os.path.join(app.config['UPLOAD_FOLDER'],fastq_file.filename)
					in_fasta_file = os.path.join(app.config['UPLOAD_FOLDER'],fasta_file.filename)
					phread_score = int(request.form.get('phred_score'))
					distance_5_prime = int(request.form.get('distance_5_prime'))
					distance_3_prime = int(request.form.get('nucleotide_match'))
					barcode = request.form.get('barcode').split(",")
					variable_sites_number = int(request.form.get('variable_sites_number'))
				except:
					raise CE.FileNotFound("invalid file, please check your files")

				out_path = FVS.main(in_fasta_file,in_fastq_file,out_path_name,True,True,phread_score,distance_5_prime,distance_3_prime,int(variable_sites_number),barcode)
				print(out_path)

				if os.path.exists(out_path):
					return send_file(out_path, as_attachment=True)
				else:
					raise CE.UnknownError("The program has failed to return a valid file path")

			
	@app.errorhandler(Exception)
	def handle_invalid_input_error(error):
		# Render a custom error page for invalid input
		return render_template('error.html', message=str(error)), 400

	@app.route('/acknowledgement')
	def acknowledgement():
		# This can be a page where users can download files
		return render_template('acknowledgement_page.html')
	
	@app.route('/contact')
	def contact():
		# This can be a page where users can download files
		return render_template('Contacts.html')
	
	@app.route('/help')
	def help():
		# This can be a page where users can download files
		return render_template('Help_page.html')
	
	
	@app.route('/startpage')
	def downloads():
		# This can be a page where users can download files
		return render_template('startpage.html')
	
	@app.route('/build_primer',methods=['GET', 'POST'])
	def build_primer():
		if request.method == 'GET':
			return render_template("build_primer.html")
		if request.method == 'POST': 

			try:
				fasta_file = request.files['file1']
				fasta_file.save(os.path.join(app.config['UPLOAD_FOLDER'],fasta_file.filename))
				out_path_name = os.path.join(app.config['UPLOAD_FOLDER'],request.form.get('file_name'))
			except:
				raise CE.FileNotFound("invalid file, please check your files")

			out_path = BP.main(os.path.join(app.config['UPLOAD_FOLDER'],fasta_file.filename),out_path_name)

			if os.path.exists(out_path):
				return send_file(out_path, as_attachment=True)
			else:
				raise CE.UnknownError("The program has failed to return a valid file path")


	

	UPLOAD_FOLDER = '/Volumes/rdss_hhaim/LAB PROJECTS/Sam/FQ_VSearch/uploads/'

	@app.route('/preference',methods=['GET', 'POST']) 
	def preference():
		if request.method == 'GET':
			return render_template("preference.html")
		if request.method == 'POST': 
			try:
				file_fasta_1 = request.files['file_fasta_1']
				file_fastq_1 = request.files['file_fastq_1']
				file_fastq_2 = request.files['file_fastq_2'] 
				file_fasta_1.save(os.path.join(app.config['UPLOAD_FOLDER'],file_fasta_1.filename)) 
				file_fastq_1.save(os.path.join(app.config['UPLOAD_FOLDER'],file_fastq_1.filename)) 
				file_fastq_2.save(os.path.join(app.config['UPLOAD_FOLDER'],file_fastq_2.filename))
			except:
				raise CE.FileNotFound("invalid file, please check your files")
			

			barcode_1 = request.form.get('barcode_1').split(",")
			barcode_2 = request.form.get('barcode_2').split(",")
			out_path = os.path.join(app.config['UPLOAD_FOLDER'],request.form.get('file_name'))
			phread_score = int(request.form.get('phred_score'))
			distance_5_prime = int(request.form.get('distance_5_prime'))
			distance_3_prime = int(request.form.get('nucleotide_match'))

			variable_sites_number = int(request.form.get('variable_sites_number0'))
			print(variable_sites_number)

			wild_type_amino_acid = request.form.get('wild_type_amino_acid_1')
			wild_type_amino_acid_2 = request.form.get('wild_type_amino_acid_2')

			if variable_sites_number > 1:
				wild_type_amino_acid = f"{wild_type_amino_acid},{wild_type_amino_acid_2}"

			print("wild_type_amino_acid",wild_type_amino_acid,wild_type_amino_acid_2)

			file_fasta_1_filename = os.path.join(app.config['UPLOAD_FOLDER'],file_fasta_1.filename)
			file_fastq_1_filename = os.path.join(app.config['UPLOAD_FOLDER'],file_fastq_1.filename)
			file_fastq_2_filename = os.path.join(app.config['UPLOAD_FOLDER'],file_fastq_2.filename)
		
			pre_amino_dict = FVS.main(file_fasta_1_filename,file_fastq_1_filename,"",True,True,phread_score,distance_5_prime,distance_3_prime,int(variable_sites_number),barcode_1)

			post_amino_dict = FVS.main(file_fasta_1_filename,file_fastq_2_filename,"",True,True,phread_score,distance_5_prime,distance_3_prime,int(variable_sites_number),barcode_2)

			print(out_path)

			out_path = RSC.main(wild_type_amino_acid ,pre_amino_dict,post_amino_dict,out_path)
			
			if os.path.exists(out_path):
				return send_file(out_path, as_attachment=True)
			else:
				raise CE.UnknownError("The program has failed to return a valid file path")

	app.run(port=5000, use_reloader=False, host='127.0.0.1')  # Turn off reloader to avoid multiple instances of the app

def start_gui():
    # Start the webview window after Flask is started.
    webview.create_window("EZDMS", 'http://127.0.0.1:5000/startpage', js_api=JSAPI(), width=200, height=200)
    webview.start()

if __name__ == '__main__':
	# Start Flask in a separate thread
	flask_thread = threading.Thread(target=start_flask)
	flask_thread.daemon = True  # Make the thread daemon so it exits when the main program exits
	flask_thread.start()
	time.sleep(2)
	# Start the webview GUI
	start_gui()

    # When the webview window is closed, the program will exit gracefully
	sys.exit()  # Ensure that the main thread exits properly after the window is closed