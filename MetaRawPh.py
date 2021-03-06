#!/usr/bin/env python


##### ##### ##### ##### ##### ##### ##### #####
#                                             #
#             MetaRawPh.py v.1.0              #
#                                             #
#        This is a script for doing the       #
#       gene centric analyses in Joel's       #
#               Honours project.              #
#                                             #
##### ##### ##### ##### ##### ##### ##### #####

##### Messsages #####
intro = '''
## MetaRawPh.py  v.0.00000001        ##
## Searches raw sequences for genes  ##
## Joel Boyd, Honours 2014           ##
'''

import argparse
import re
try:
    from Bio import SeqIO
except ImportError:
    print "Whoops! MetaRawPhy.py needs a whole ton of stuff installed! pynast biopython emboss fxtract hmmer seqmagick pplacer/2.6.32-Erick-Hack"
    exit(1)
import math
import glob
import subprocess
import os.path
from datetime import datetime

##### Classes #####

class Cl(object):
  pass

##### Input Files #####

parser = argparse.ArgumentParser(description='''MetaRawPh.py  v.0.1'''
                          , epilog='Joel Boyd - Honours 2014.')
parser.add_argument('-f', metavar = 'forward read or fasta file', type = str, nargs = 1 , help = 'Froward raw sequence file in .fa, or .fq.gz format.', required = True)
parser.add_argument('-r', metavar = 'reverse read', type = str, nargs = 1, help = 'Optional reverse raw sequence file in .fa, or .fq.gz format.', default = argparse.SUPPRESS)
parser.add_argument('-t', metavar = 'prot or dna', type = str, nargs = 1, help = 'dna (like 16S) or prot (like mcrA)', choices = ['prot','dna'], required = True)
parser.add_argument('-m', metavar = 'hmm_file', type = str, nargs = 1, help = 'HMM file', required = True)
parser.add_argument('-c', metavar = 'reference_package', type = str, nargs = 1, help = 'Reference package of gene family', required = True)
parser.add_argument('-g', metavar = 'gg_database', type = str, nargs = 1, help = 'Aligned gg database for pynast alignment (dna sequences only)', default = argparse.SUPPRESS)
parser.add_argument('-v', action = 'version', version = 'MetaRawPhy.py v.0.1')

cl = Cl()
args=parser.parse_args(namespace = cl)

##### Functions #####

# splits fasta file for hmmaligning
def fasta_parser(to_split):
  num_seqs = 1
  all_fasta = list(SeqIO.parse(to_split, 'fasta'))
  num_files = int(math.ceil(len(all_fasta)/float(num_seqs)))
  i = 1
  while i <= num_files:
    start = int(i-1) * int(num_seqs)
    end = int(i) * int(num_seqs)
    filename = to_split.split('.')[:-1][0] + '_' + str(i) + '.g.fasta'
    SeqIO.write(all_fasta[start:end], filename, 'fasta')
    i += 1

# run nhmmer
def nhmmer(hmm, raw_seq_file):
  fasta_match = re.search('\.fa$', raw_seq_file)
  fastqgz_match = re.search('\.fq\.gz$', raw_seq_file)
  path = './'+replace_name(raw_seq_file, '')+'.1.hmmout.csv'
  if os.path.isfile(path) == True:
    hmmout_table_title = replace_name(raw_seq_file, '.2.hmmout.csv')
    if fasta_match:
      subprocess.check_call(["/bin/bash", "-c", " nhmmer --tblout " +hmmout_table_title+ " " +hmm+ " <(sed 's/:/_/g' " +raw_seq_file+ ") 2>&1 > /dev/null"])
    if fastqgz_match:
      subprocess.check_call(["/bin/bash", "-c", " nhmmer --tblout " +hmmout_table_title+ " " +hmm+ " <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat " +raw_seq_file+ " | sed 's/:/_/g')) 2>&1 > /dev/null"])
  elif os.path.isfile(path) == False:
    hmmout_table_title = replace_name(raw_seq_file, '.1.hmmout.csv')
    if fasta_match:
      print 'hello'
      subprocess.check_call(["/bin/bash", "-c", " nhmmer --tblout " +hmmout_table_title+ " " +hmm+ " <(sed 's/:/_/g' " +raw_seq_file+ ") 2>&1 > /dev/null"])
    if fastqgz_match:
      subprocess.check_call(["/bin/bash", "-c", " nhmmer --tblout " +hmmout_table_title+ " " +hmm+ " <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat " +raw_seq_file+ " | sed 's/:/_/g')) 2>&1 > /dev/null"])
  else:
    print 'ERROR: Suffix on %s not recegnised. Please submit an .fq.gz or .fa file\n' % (raw_seq_file)
    exit(1)


# run hmmsearch
def hmmsearch(hmm, raw_seq_file):
  fasta_match = re.search('\.fa$', raw_seq_file)
  fastqgz_match = re.search('\.fq\.gz$', raw_seq_file)
  path = './'+replace_name(raw_seq_file, '')+'.1.hmmout.csv'
  if os.path.isfile(path) == True:
    hmmout_table_title = replace_name(raw_seq_file, '.2.hmmout.csv')
    if fasta_match:
      subprocess.check_call(["/bin/bash", "-c", " hmmsearch --domtblout " +hmmout_table_title+ " " +hmm+ " <(getorf -sequence  <(sed 's/:/_/g' " +raw_seq_file+ ") -outseq >(cat) -minsize 98 2>/dev/null) 2>&1 > /dev/null"])
      return
    if fastqgz_match:
      subprocess.check_call(["/bin/bash", "-c", " hmmsearch --domtblout " +hmmout_table_title+ " " +hmm+ " <(getorf -sequence <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat " +raw_seq_file+ " | sed 's/:/_/g')) \
      -outseq >(cat) -minsize 98 2>/dev/null) 2>&1 > /dev/null"])
      return
    else:
      print 'ERROR: Suffix on %s not recegnised. Please submit an .fq.gz or .fa file\n' % (raw_seq_file)
      exit(1)
  else:
    hmmout_table_title = replace_name(raw_seq_file, '.1.hmmout.csv')
    if fasta_match:
      subprocess.check_call(["/bin/bash", "-c", " hmmsearch --domtblout " +hmmout_table_title+ " " +hmm+ " <(getorf -sequence  <(sed 's/:/_/g' " +raw_seq_file+ ") -outseq >(cat) -minsize 98 2>/dev/null) 2>&1 > /dev/null"])
      return
    if fastqgz_match:
      subprocess.check_call(["/bin/bash", "-c", " hmmsearch --domtblout " +hmmout_table_title+ " " +hmm+ " <(getorf -sequence <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat " +raw_seq_file+ " | sed 's/:/_/g')) \
      -outseq >(cat) -minsize 98 2>/dev/null) 2>&1 > /dev/null"])
      return
    else:
      print 'ERROR: Suffix on %s not recegnised. Please submit an .fq.gz or .fa file\n' % (raw_seq_file)
      exit(1)

# run pynast
def pynast(ts_file, gg_db_path):
  subprocess.check_call('pynast -l 0 -i ' +replace_name(ts_file, '.ts.fasta')+ ' -t ' +gg_db_path+ ' -a ' +replace_name(ts_file, '.aln.fasta'), shell=True)

# run pplacer
def pplacer(refpkg, ts_file):
  subprocess.check_call('pplacer --verbosity 0 -c ' +refpkg+ ' ' +ts_file, shell=True)

# run guppy classify
def guppy_class(rpkg, jplace_file):
  subprocess.check_call('guppy classify -c ' +rpkg+ ' '+jplace_file+' > ' +replace_name(jplace_file, '.guppy'), shell=True)

# Removes unneeded files
def rubbish(unneeded_file):
  subprocess.check_call("rm " +str(unneeded_file)+ ' ', shell=True)

# HMM aligning
def hmmalign(hmm, sequencefile):
  subprocess.check_call("mkdir "+replace_name(cl.f[0], ''), shell=True)
  base_name = replace_name(sequencefile, '')
  subprocess.check_call("mv " +base_name+".*  "+replace_name(cl.f[0], ''), shell = True)
  if cl.t[0] == 'prot':
    fasta_parser((replace_name(cl.f[0], '')+'/'+replace_name(sequencefile, '.orf.fasta')))
  elif cl.t[0] == 'dna':
    fasta_parser((replace_name(cl.f[0], '')+'/'+replace_name(sequencefile, '.ts.fasta')))
  path = './'+replace_name(cl.f[0], '')+'/*.g.fasta'
  files = glob.glob(path)
  for file in files:
    subprocess.check_call('hmmalign --trim -o '+file+'.sto ' +hmm+ ' ' +file, shell = True)
  path = './'+replace_name(cl.f[0], '')+'/*.sto'
  files = glob.glob(path)
  for file in files:
    subprocess.check_call('seqmagick convert ' +file+ ' ' +file+'conv_.fa', shell=True)
  path = './'+replace_name(cl.f[0], '')+'/*conv_.fa'
  files = glob.glob(path)
  for file in files:
    subprocess.check_call("sed 's/[a-z]//g' " +file+ " > " +file+'.ready', shell=True)
  subprocess.check_call("rm "+replace_name(cl.f[0], '')+"/*.sto "+replace_name(cl.f[0], '')+"/*conv_.fa "+replace_name(cl.f[0], '')+"/*.fasta", shell=True)
  path = './'+replace_name(cl.f[0], '')+'/*.ready'
  files = glob.glob(path)
  out_aligned_name = replace_name(sequencefile, '.aln.fasta')
  title_file_open = open(out_aligned_name, 'w')
  for file in files:
    for line in open(file):
      title_file_open.write(line)
  title_file_open.close()
  subprocess.check_call('rm '+replace_name(cl.f[0], '')+'/*.ready', shell=True)

# process hmmsearch/nhmmer results into 'ts_file's
def csv_to_titles_for(hmmout, type_, raw_seq_file):
  raw_titles = []
  line_counter = 0
  for line in open(hmmout):
    match = re.search('^\#', line)
    if match:
      continue
    else:
      line_counter += 1
      split = line.split(' ', 1)
      new_title = str(split[0].replace('_',':'))
      if type_ == 'dna':
        raw_titles.append(new_title)
      if type_ == 'prot':
        if 'FCC' in line:
          raw_titles.append(new_title[:-4]+'/1')
        else:
          raw_titles.append(new_title[:-2])
  if line_counter == 0:
    time = datetime.now().strftime('%H:%M:%S')
    print '['+time+']: 0 Reads found!'
    exit(1)
  else:
    time = datetime.now().strftime('%H:%M:%S')
    print '['+time+']: Found %s read(s) in %s.' % (str(line_counter), replace_name(raw_seq_file, ''))
  file_name = replace_name(hmmout, '.titles.1.txt')
  title_file_open = open(file_name, 'w')
  for title in raw_titles:
    title = str(title) + '\n'
    title_file_open.write(title)
  title_file_open.close()

def csv_to_titles_rev(hmmout, type_, raw_seq_file):
  raw_titles = []
  line_counter = 0
  for line in open(hmmout):
    match = re.search('^\#', line)
    if match:
      continue
    else:
      line_counter += 1
      split = line.split(' ', 1)
      new_title = str(split[0].replace('_',':'))
      if type_ == 'dna':
        raw_titles.append(new_title)
      if type_ == 'prot':
        if 'FCC' in line:
          raw_titles.append(new_title[:-4]+'/1')
        else:
          raw_titles.append(new_title[:-2])
  if line_counter == 0:
    time = datetime.now().strftime('%H:%M:%S')
    print '['+time+']: 0 Reads found!'
    exit(1)
  else:
    time = datetime.now().strftime('%H:%M:%S')
    print '['+time+']: Found %s read(s) in %s.' % (str(line_counter), replace_name(raw_seq_file, ''))
  file_name = replace_name(hmmout, '.titles.2.txt')
  title_file_open = open(file_name, 'w')
  for title in raw_titles:
    title = str(title) + '\n'
    title_file_open.write(title)
  title_file_open.close()

def extract_from_raw_reads(raw_seq_file, hmm):
  all_titles = []
  for_titles = []
  file_name = replace_name(raw_seq_file, '.titles.txt')
  x = open('./'+replace_name(cl.f[0], '')+'.titles.1.txt', 'r')
  for line in x:
    for_titles.append(line)
  x.close()
  try:
    x = open('./'+replace_name(cl.f[0], '')+'.titles.2.txt', 'r')
    for line in x:
      if line in for_titles:
        all_titles.append(line)
    x.close()
  except IOError:
    all_titles = for_titles
  title_file_open = open(file_name, 'w')
  for x in all_titles:
    title_file_open.write(x)
  title_file_open.close()
  file_name = replace_name(raw_seq_file, '.ts')
  match = re.search('\.fa$', raw_seq_file)
  if match:
    time = datetime.now().strftime('%H:%M:%S')
    print '['+time+']: Extracting reads.'
    subprocess.check_call("fxtract -H -f " +replace_name(raw_seq_file, '.titles.txt')+ " " +raw_seq_file+ " > " +file_name+ " ", shell=True)
  elif '.fq.gz' in raw_seq_file:
    time = datetime.now().strftime('%H:%M:%S')
    print '['+time+']: Extracting reads.'
    subprocess.check_call("fxtract -H -f " +replace_name(raw_seq_file, '.titles.txt')+ " " +raw_seq_file+ " | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > " +file_name+ " ", shell=True)
  file_name2 = replace_name(raw_seq_file, '.ts.fasta')
  fasta_file_open = open(file_name2, 'w')
  for line in open(file_name):
    if '>' in line:
      line = line.replace(':', '_')
      fasta_file_open.write(line)
    else:
      line = line
      fasta_file_open.write(line)
  fasta_file_open.close()
  file_name = replace_name(raw_seq_file, '.ts')
  file_name2 = replace_name(raw_seq_file, '.ts.fasta')
  fasta_file_open = open(file_name2, 'w')
  for line in open(file_name):
    if '>' in line:
      line = line.replace(':', '_')
      fasta_file_open.write(line)
    else:
      fasta_file_open.write(line)
  fasta_file_open.close()
  if cl.t[0] == 'dna':
    return
  raw_orf_title = replace_name(raw_seq_file, '.orf')
  subprocess.check_call('getorf -sequence ' +file_name2+ ' -outseq ' +raw_orf_title+ ' -minsize 98 2>/dev/null', shell=True)
  hmm_out_title = replace_name(raw_seq_file, '.tmp.csv')
  subprocess.check_call("hmmsearch --tblout "+hmm_out_title+" "+hmm+" "+raw_orf_title+" 2>&1 > /dev/null", shell=True)
  raw_titles = []
  for line in open(hmm_out_title):
    hmm_table_parser = re.search('^#', line)
    if hmm_table_parser:
      continue
    else:
      split = line.split(' ', 1)
      raw_titles.append(split[0])
  file_name = replace_name(raw_seq_file, '.orf.titles')
  title_file_open = open(file_name, 'w')
  for title in raw_titles:
    title = str(title) + '\n'
    title_file_open.write(title)
  title_file_open.close()
  hmm_out_title = replace_name(raw_seq_file, '.orf.fasta')
  subprocess.check_call('fxtract -H -f '+file_name+' '+raw_orf_title+' > '+hmm_out_title, shell=True)

def replace_name(old_title,new_suffix):
  if '/' in old_title:
    split = old_title.split('/')
    title2 = split[len(split) - 1]
    if '.' in title2:
      split = title2.split('.')
      return split[0] + new_suffix
    else:
      return 'Suffix on '+old_title+' not recognized.'
  else:
    if '.' in old_title:
      split = old_title.split('.')
      return split[0] + new_suffix
    else:
      return old_title


 ##### Running Script #####

print intro,
seq_file_title = replace_name(cl.f[0], '')
hmm_file_title = replace_name(cl.m[0], '.hmm')



if cl.t[0] == 'prot':
  date = datetime.now().strftime('%Y/%m/%d')
  time = datetime.now().strftime('%H:%M:%S')
  print '\n            '+date
  print '['+time+']: Searching %s using %s.' % (seq_file_title, hmm_file_title)
  hmmsearch(cl.m[0], cl.f[0])
  csv_to_titles_for(replace_name(cl.f[0], '.1.hmmout.csv'), cl.t[0], cl.f[0])
  try:
    hmmsearch(cl.m[0], cl.r[0])
    csv_to_titles_rev(replace_name(cl.r[0], '.2.hmmout.csv'), cl.t[0], cl.f[0])
  except AttributeError:
    time = datetime.now().strftime('%H:%M:%S')
    print '['+time+']: No reverse read found. Continuing with forward read only.'
  extract_from_raw_reads(cl.f[0], cl.m[0])
  time = datetime.now().strftime('%H:%M:%S')
  print '['+time+']: Aligning to hmm.'
  hmmalign(cl.m[0], replace_name(cl.f[0], '.ts.fasta'))
  time = datetime.now().strftime('%H:%M:%S')
  print '['+time+']: Placing in refpkg tree.'
  pplacer(cl.c[0], replace_name(cl.f[0], '.aln.fasta'))
  time = datetime.now().strftime('%H:%M:%S')
  print '['+time+']: Creating Guppy file.'
  guppy_class(cl.c[0], replace_name(cl.f[0],'.aln.jplace'))
  time = datetime.now().strftime('%H:%M:%S')
  print '['+time+']: Finished %s!' % (seq_file_title)
  exit(1)

if cl.t[0] == 'dna':
  date = datetime.now().strftime('%Y-%m-%d')
  time = datetime.now().strftime('%H:%M:%S')
  print '\n            '+date
  print '['+time+']: Searching %s using %s.' % (seq_file_title, hmm_file_title)
  nhmmer(cl.m[0], cl.f[0])
  csv_to_titles_for(replace_name(cl.f[0], '.1.hmmout.csv'), cl.t[0], cl.f[0])
  try:
    nhmmer(cl.m[0], cl.r[0])
    csv_to_titles_rev(replace_name(cl.r[0], '.2.hmmout.csv'), cl.t[0], cl.f[0])
  except AttributeError:
    time = datetime.now().strftime('%H:%M:%S')
    print '['+time+']: No reverse read found. Continuing with forward read only.'
  extract_from_raw_reads(cl.f[0], cl.m[0])
  time = datetime.now().strftime('%H:%M:%S')
  print '['+time+']: Aligning to hmm.'
  pynast(cl.f[0], cl.g[0])
  time = datetime.now().strftime('%H:%M:%S')
  print '['+time+']: Placing in refpkg tree.'
  pplacer(cl.c[0], replace_name(cl.f[0],'.aln.fasta'))
  time = datetime.now().strftime('%H:%M:%S')
  print '['+time+']: Creating Guppy file.'
  guppy_class(cl.c[0], replace_name(cl.f[0], '.aln.jplace'))
  time = datetime.now().strftime('%H:%M:%S')
  print '['+time+']: Finished %s!' % (seq_file_title)
  exit(1)


