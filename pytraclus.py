#! /usr/bin/python3
# -*- coding: UTF-8 -*-

# Thanks due to Pierre Poulain
# Olivier July 2012
import string
import re
import os
import sys
import glob
import cgi
import time
import zipfile
import random
import subprocess
import datetime
import HTML #import HTML.py (in same folder) this script takes care of most of
			#communication with the index.html sheet, does it ?

# data
#=======================================================================

# genetic code
genetic_code = {"GCT":"A","GCC":"A","CAC":"H","GCA":"A","GCG":"A","TGT":"C","TGC":"C","GAT":"D","GAC":"D","GAA":"E","GAG":"E","TTT":"F","TTC":"F","GGT":"G","GGC":"G","GGA":"G","GGG":"G","CAT":"H","ATT":"I","ATC":"I","ATA":"I","AAA":"K","AAG":"K","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L","ATG":"M","AAT":"N","AAC":"N","CCT":"P","CCC":"P","CCA":"P","CCG":"P","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R","AGG":"R","CAA":"Q","CAG":"Q","TCT":"S","TCC":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GTT":"V","GTC":"V","GTA":"V","GTG":"V","TGG":"W","TAT":"Y","TAC":"Y","TAA":"*","TGA":"*","TAG":"*"} 

toggle= ""





#=======================================================================
# FUNCTIONS
#=======================================================================

leet = str.maketrans('ATGCN', 'TACGN') #ROFL THX Doug !

def isNucleicSequence(seq):
	# return True if the input sequence is a nucleic sequence
	#if len(seq) == 0:
	#	return False
	for i in seq:
		if i not in ["A", "T", "C", "G", "N"]:
			return False
	return True



#======================================================================#
# MAIN PROGRAM                                                         #
#======================================================================#


# debug
#=======================================================================
import cgitb
cgitb.enable()
sys.stderr = sys.stdout
# end of debugging

# print HTML header
#=======================================================================
HTML.Header()

# get rid of previous results
#=======================================================================
files =glob.glob("pzu/*.*")
for f in files: 
	os.remove (f)


# Recover zip 
#=======================================================================

print ("<h3>Reading data from HTML form</h3>")

form = cgi.FieldStorage()
seq_file = form["seq_0"]
print
uploaded_name = os.path.basename(seq_file.filename)
print   ("""uploaded_name here : %s """ % uploaded_name)

# write file on server disk
#=======================================================================
toggle = ""
try:
	zip_file = zipfile.ZipFile(seq_file.file)
	#print("""<br>PROOOOOT SUCCESS""")
	print ("""<br />file <span class="seq">%s</span> successfully uploaded <br />""" %(uploaded_name))

except:
	print ("Uploaded file is not a zip file; hopefully a text file with Fasta sequences !!")
	toggle="txt"


# treatment of zip file
#=======================================================================

if toggle != "txt" :
	zip_name = "./pzu" + "/seq.zip"
	print   ("""<br>zip_name here : %s """ % zip_name)
	# write a binary file
	f_out = open(zip_name, "wb")
	# be sure to be at the start of the file
	seq_file.file.seek(0)
	f_out.write(seq_file.file.read())
	f_out.close()
	


# inflate zip 
#=======================================================================
command_line = "unzip %s -d %s" %("./pzu/seq.zip", "./pzu")

proc = subprocess.Popen(command_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	
if proc.wait() != 0:
	print ("<hr />")
	print (proc.stderr.read())
	print ("<br />")
	print ("ERROR while decompressing zip file")
	sys.exit(1)
	
# list sequence files
#===================================================================

harvest = glob.glob("pzu/*seq")
# print(harvest)

command_line = 'find pzu -name "*.*"' 
proc = subprocess.Popen(command_line, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
if proc.wait() != 0:
	print ("<hr />")
	print (proc.stderr.read())
	print ("<br />")
	print ("ERROR while listing inflated file")
	sys.exit(1)

print ("""<br>%d files found in zip archive as .seq files... <br />""" % (len(harvest)))

# recover 5' and 3' signatures and used checkboxes
#=======================================================================

signatures = {} #why a dictionnary ???


if 'ch_b0' in form : #ch_b0 is for reverse comlementation of 5' signature
	ch_b0 = form["ch_b0"].value
else :
	ch_b0 = ""

try :
	'5_signature' in form    
	signature_5 = form["5_signature"].value.upper()
	if signature_5 != "":
		print("""<span class="warning">:-0 </span> 5 signature %s""" % str(signature_5))
	else: 
		print("""<span class="warning">:-0 </span> 5 signature box empty """ )
		print("""<br><span class="warning">If you are not interested in clustering but want only see nucleotide sequences and translations you might go for HULK """)
		problem_5 = "on"
	print("""   </br><br />""")
	# ci dessous inutile tout se fait sur l'index html
	if isNucleicSequence(signature_5) == False:
		print("""<br><span class="warning">5' signature contains weird characters, go back with your browser and correct the entry""")
		sys.exit(0)
		

except :
	print ("""<span class="warning">ERROR ERROR ERROR </span>: 5' signature  is not a nucleic sequence </br><br />""")
	#il faut check la vérité des signatures ZOB inside ?
	#sys.exit(1)
	problem="on"

	 
if ch_b0=="rev_comp_signature":
	signature_5c=str(signature_5).translate(leet)
	signature_5rc=signature_5c[::-1]
	signatures["5' signature"]=signature_5rc
	print ("""Entered 5' signature was reverse complemented, as requested: <span class="seq">%s</span><br />""") %signature_5rc
else :
	pass


'3_signature' in form    
signature_3 = form["3_signature"].value.upper()
if signature_3 != "":
	print("""<span class="warning">:-0 </span> 3 signature %s""" % str(signature_3))
else :
	print("""<span class="warning">:-0 </span> 3 signature box empty""")
	problem_3 = "on"
	
print("""   </br><br />""")

# read sequence files, store sequences dispaly them on screen
#======================================================================= 
sequences = {}
if toggle != "txt":      
	nuc_name = "pzu" + "/a" # "+/a" or whatever is mandatory * is all right, . not :-0
	# pzu is name of directory where zip has been inflated
	# nuc_name = full path of file containing sequences
	print("""<span class="warning">:-0 </span> nuc_name here : %s""" % nuc_name)  
	
	fo_ut=open(nuc_name, "a")
	for f_name in harvest:
		print("""<br><br><span class="warning">:-0 </span> f_name : %s""" % f_name)
		f_in = open(f_name, "r")
		lines = f_in.readlines()
		print("""<br><span class="warning">:-0 </span> lines one after the other  : %s""" % lines)
		# we get a list : first item is first line i.e. the stuff given by GATC Biotech e.g. >26974237.seq-ID:PCR2 MDCJ Lu Triple... etc. Second item is corresponding nucleotide sequence.
		f_in.close()
		# create sequence name from file name
		f_name_NOT_short = f_name.split("/")
		print("""<br><span class="warning">:-0 </span> f_name_NOT_short after split : a list !! : %s""" % f_name_NOT_short)
		f_name_short = f_name.split("/")[-1]
		print("""<br><span class="warning">:-0 </span> f_name_short  : %s""" % f_name_short)
		seq_name = "_".join(f_name_short.split(".")[0:-1])
		print("""<br><span class="warning">:-0 </span> seq_name  : %s""" % seq_name)
		# check for correct fasta format
		if ">" != lines[0][0]:
			print ("""<span class="warning">WARNING: no > found at the beginning of %s [discarded]</span> <br />""") %(f_name_short)
			continue
		# read sequence
		seq = ""
		for line in lines[1:]:
			if ">" in line:
				print ("""<span class="warning">WARNING: more than one sequence in %s? [discarded]</span> <br />""" %(f_name_short))
				break
			seq += line.strip().upper()
		if seq == "":
			print ("""<span class="warning">WARNING: sequence appears empty in %s [discarded]</span> <br />""" %(f_name_short))
			continue

		# finally accept sequence
		
		
		sequences[seq_name] = seq
		# print (seq_name)
		# print (sequences[seq_name])
		
		# print ("<hr />")
		# sequences is a dictionary seq_name is key, seq is value
	fo_ut.close()
	fi_n = open(nuc_name, "r")
	fi_n.seek(0, 0)
	nuc_content = fi_n.read()
	print ("<br/>Below meaningful nucleotides sequences extracted from archive :")
	
	zogzog= ""
	current_date = str(datetime.date.today())
	for key, value in sequences.items() :
		zogdata=[">",key,"-",current_date,"\n",value,"\n"]
		datazog="".join(zogdata)
		#fo_ut.write(datazog)
		zogzog = zogzog + datazog
	
	print ( """<textarea  readonly cols =112  rows =10 style = "background-color: #BEF781;">%s</textarea >"""%zogzog)
	fi_n.close()
	
	print ("""<br /><span class="comment"style ="color: #204ED8">Note that areas with results are read-only, text can be copied and pasted to a text editor, motives may be searched using the search funtion of your browser (ctrl-F), whole area (which is scrollable) may be enlarged in most browsers by dragging handle (lower right corner).<br /><br /></span>""")
	
	print ("<h3>Looking for signature in sequences, and trimming sequence ends</h3>")

if 'ch_b' in form : #ch_b0 is for reverse comlementation of 5' signature
	ch_b = form["ch_b"].value
	
for name, seq in sequences.items():
	print ("""<br><span class="warning">name is  %s </span> <br />""" %name)
	print ("""<br><span class="warning">seq is  %s </span> <br />""" %seq)
	y = seq.find(str(signature_5).upper())
	print( """here y in seq: %d	<br>""" % y)
	seq_rc = str(seq).translate(leet)
	seq_rc = seq_rc[::-1]
	x = seq_rc.find(signature_5)
	print( """here x in seq_rc: %d	<br>""" % x)
		
	
