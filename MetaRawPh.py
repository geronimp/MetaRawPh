#!/usr/bin/env python


##### ##### ##### ##### ##### ##### ##### #####
#                                             #
#             MetaRawPh.py v.1.0              #
#                                             #
#    This is a master script for doing the    #
#       gene centric analyses in Joel's       #
#               Honours project.              #
#                                             #
##### ##### ##### ##### ##### ##### ##### #####

##### Messsages #####
intro = '''
## MetaRawPh.py  v.1.0               ##
## Searches raw sequences for genes  ##
## Joel Boyd, Honours 2014           ##
'''

import argparse
import re
try:
    from Bio import SeqIO
except ImportError:
    print "Whoops! MetaRawPhy.py needs a whole ton of stuff installed! biopython emboss fxtract hmmer pv seqmagick pplacer/2.6.32-Erick-Hack"
    exit(1)
import math
import glob
import subprocess

##### Classes #####

class Cl(object):
  pass

##### Input Files #####

parser = argparse.ArgumentParser(description='''MetaRawPh.py  v.1.0
                          Searches raw sequences for genes'''
                          , epilog='Joel Boyd - Honours 2014.')
parser.add_argument('-s', metavar = 'raw_sequence_file', type=str, help='Raw sequence file in fasta, or fq.gz format.', required=True)
parser.add_argument('-t', metavar = 'prot or dna', type=str, nargs = 1, help = 'dna or prot', choices=['prot','dna'], required=True)
parser.add_argument('-m', metavar = 'hmm_file', type=str, nargs = 1, help = 'hmm file built from database', required=True)
parser.add_argument('-c', metavar = 'reference_package', type=str, nargs = 1, help = 'reference package of same database used to build hmm', required=True)
parser.add_argument('-g', metavar = 'gg_database', type=str, nargs = '?', help = 'Aligned gg database for pynast alignment', default=argparse.SUPPRESS)
parser.add_argument('-v', action='version', version='MetaRawPhy.py v.0.1')

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
  out_title = newname(raw_seq_file, '.hmmout.csv')
  subprocess.check_call(["/bin/bash", "-c", " nhmmer --tblout " +out_title+ " " +hmm+ " \
  <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(pv --eta -p " +raw_seq_file+ " | zcat | sed 's/:/_/g')) 2>&1 > /dev/null"])

# run hmmsearch
def hmmsearch(hmm, raw_seq_file):
  out_title = newname(raw_seq_file, '.hmmout.csv')
  subprocess.check_call(["/bin/bash", "-c", " hmmsearch --domtblout " +out_title+ " " +hmm+ " \
  <(getorf -sequence <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(pv --eta --timer -p " +raw_seq_file+ " | zcat | sed 's/:/_/g')) \
  -outseq >(cat) -minsize 99) 2>&1 > /dev/null"])

# run pynast
def pynast(ts_file, gg_db_path):
  subprocess.check_call('pynast -l 0 -i ' +newname(ts_file, '.ts.fasta')+ ' -t ' +gg_db_path+ ' -a ' +newname(ts_file, '.aln.fasta'), shell=True)

# run pplacer
def pplacer(refpkg, ts_file):
  subprocess.check_call('pplacer --verbosity 0 -c ' +refpkg+ ' ' +ts_file, shell=True)

# run guppy classify
def guppy_class(rpkg, jplace_file):
  subprocess.check_call('guppy classify -c ' +rpkg+ ' '+jplace_file+' > ' +newname(jplace_file, '.guppy'), shell=True)

# Removes unneeded files
def rubbish(unneeded_file):
  subprocess.check_call("rm " +str(unneeded_file)+ ' ', shell=True)

# HMM aligning
def hmmalign(hmm, sequencefile):
  subprocess.check_call("mkdir "+newname(cl.s, ''), shell=True)
  base_name = newname(sequencefile, 'keep')[:-8]
  subprocess.check_call("mv " +base_name+"*  "+newname(cl.s, ''), shell=True)
  if cl.t[0] == 'prot':
    fasta_parser((newname(cl.s, '')+'/'+newname(sequencefile, '.orf.fasta')))
  elif cl.t[0] == 'dna':
    fasta_parser((newname(cl.s, '')+'/'+newname(sequencefile, '.ts.fasta')))
  path = './'+newname(cl.s, '')+'/*.g.fasta'
  files = glob.glob(path)
  for file in files:
    subprocess.check_call('hmmalign --trim -o '+file+'.sto ' +hmm+ ' ' +file, shell=True)
  path = './'+newname(cl.s, '')+'/*.sto'
  files = glob.glob(path)
  for file in files:
    subprocess.check_call('seqmagick convert ' +file+ ' ' +file+'.fa', shell=True)
  path = './'+newname(cl.s, '')+'/*.fa'
  files = glob.glob(path)
  for file in files:
    subprocess.check_call("sed 's/[a-z]//g' " +file+ " > " +file+'.ready', shell=True)
  subprocess.check_call("rm "+newname(cl.s, '')+"/*.sto "+newname(cl.s, '')+"/*.fa "+newname(cl.s, '')+"/*.fasta", shell=True)
  out_aligned_name = newname(sequencefile, '.aln.fasta')
  path = './'+newname(cl.s, '')+'/*.ready'
  files = glob.glob(path)
  title_file_open = open(out_aligned_name, 'w')
  for file in files:
    for line in open(file):
      title_file_open.write(line)
  title_file_open.close()
  subprocess.check_call('rm '+newname(cl.s, '')+'/*.ready', shell=True)


# process hmmsearch/nhmmer results into 'ts_file's
def csv_to_titles(hmmout, raw_seq_file, hmm, type_):
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
        raw_titles.append(new_title[:-2])
  if line_counter == 0:
    print '''
0 Reads found!
    '''
    exit(1)
  else:
    print '''
Found %s read(s) in %s
Extracting reads
    ''' % (str(line_counter), newname(raw_seq_file, 'keep'))
  file_name = newname(hmmout, '.titles.txt')
  title_file_open = open(file_name, 'w')
  for title in raw_titles:
    title = str(title) + '\n'
    title_file_open.write(title)
  title_file_open.close()

  file_name = newname(raw_seq_file, '.ts')
  subprocess.check_call("fxtract -H -f " +newname(hmmout, '.titles.txt')+ " " +raw_seq_file+ " | awk '{print \">\" substr($0,2);getline;print;getline;getline}' > " +file_name+ " ", shell = True)

  file_name2 = newname(raw_seq_file, '.ts.fasta')
  fasta_file_open = open(file_name2, 'w')
  for line in open(file_name):
    if '>' in line:
      line = line.replace(':', '_')
      fasta_file_open.write(line)
    else:
      line = line
      fasta_file_open.write(line)
  fasta_file_open.close()

  file_name = newname(raw_seq_file, '.ts')
  file_name2 = newname(raw_seq_file, '.ts.fasta')
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
  raw_orf_title = newname(raw_seq_file, '.orf')
  subprocess.check_call('getorf -sequence ' +file_name2+ ' -outseq ' +raw_orf_title+ ' -minsize 99 ', shell=True)
  hmm_out_title = newname(raw_seq_file, '.tmp.csv')
  subprocess.check_call("hmmsearch --tblout "+hmm_out_title+" "+hmm+" "+raw_orf_title+" 2>&1 > /dev/null", shell=True)
  raw_titles = []
  for line in open(hmm_out_title):
    if '#' in line:
      continue
    else:
      split = line.split(' ', 1)
      raw_titles.append(split[0])
  file_name = newname(raw_seq_file, '.orf.titles')
  title_file_open = open(file_name, 'w')
  for title in raw_titles:
    title = str(title) + '\n'
    title_file_open.write(title)
  title_file_open.close()
  hmm_out_title = newname(raw_seq_file, '.orf.fasta')
  subprocess.check_call('fxtract -H -f '+file_name+' '+raw_orf_title+' > '+hmm_out_title, shell=True)


def newname(old_title,new_suffix):
  if '/' in old_title:
    split = old_title.split('/')
    title = split[len(split) - 1]
    if new_suffix == 'keep':
      return title
    elif '.' in title:
      split = title.split('.')
      return split[0] + new_suffix
    else:
      return 'Suffix on '+title+' not recognized.'
  else:
    if new_suffix == 'keep':
      return old_title
    elif '.' in old_title:
      split = old_title.split('.')
      return split[0] + new_suffix
    else:
      return 'Suffix on '+title+' not recognized.'


 ##### Running Script #####

print intro
seq_file_title = newname(cl.s, 'keep')
hmm_file_title = newname(cl.m, 'keep')



if cl.t[0] == 'prot':
  print '''
Searching %s using %s
  ''' % (seq_file_title, hmm_file_title)
  hmmsearch(cl.m, cl.s)
  csv_to_titles(newname(cl.s, '.hmmout.csv'), cl.s, cl.m, cl.t[0])
  print '''
Aligning...
  '''
  hmmalign(cl.m, newname(seq_file_title, '.ts.fasta'))
  print '''
pplacing...
  '''
  pplacer(cl.c[0], newname(seq_file_title, '.aln.fasta'))
  print '''
Creating Guppy file...
  '''
  guppy_class(cl.c[0], newname(seq_file_title, '.aln.jplace'))
  print '''
Finished %s!
  ''' % (seq_file_title)
  exit(1)

if cl.t[0] == 'dna':
  print '''
Searching %s using %s
  ''' % (seq_file_title, hmm_file_title)
  nhmmer(cl.m, cl.s)
  csv_to_titles(newname(cl.s, '.hmmout.csv'), cl.s, cl.m, cl.t[0])
  print '''
Aligning...
  '''
  pynast(cl.s, cl.g)
  print '''
pplacing...
  '''
  pplacer(cl.c[0], newname(seq_file_title, '.aln.fasta'))
  print '''
Creating Guppy file...
  '''
  guppy_class(cl.c[0], newname(seq_file_title, '.aln.jplace'))
  print '''
Finished %s!
  ''' % (seq_file_title),
  exit(1)
