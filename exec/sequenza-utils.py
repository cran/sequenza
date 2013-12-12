#!/usr/bin/env python2.7


import argparse, os, sys, gzip, math, multiprocessing, time, logging, json, gc
if sys.version_info[0] < 3:
   from itertools import izip_longest
else:
   from itertools import zip_longest as izip_longest
   from codecs import getreader, getwriter
from multiprocessing.pool import ThreadPool
from functools import partial
from multiprocessing.queues import SimpleQueue
try:
   import rpy2.robjects as robjects
   from rpy2.robjects.packages import importr
   from rpy2.rinterface import RRuntimeError
   try:
      sequenza = importr("sequenza")
      RPY2 = True
   except RRuntimeError:
      print "The R sequenza package is not installed in your R installation (" + os.environ['R_HOME'] + ")."
      RPY2 = False
except ImportError:
   RPY2 = False

VERSION = "1.0.0"
DATE    = "28 November 2013"
AUTHOR  = "Favero Francesco"
MAIL    = "favero@cbs.dtu.dk"

def check_dir(directory):
   '''
   Check the given directory and return the whole path
   or set the default to working directory.
   '''
   if directory == '':
      directory = None
   if directory == None:
      directory = os.getcwd()
   if os.path.isdir(directory):
      directory = os.path.abspath(directory)
   else:
      sys.exit("Not existing directory, please create the directory first")
   return directory

def walklevel(directory, level=1):
   '''

   '''
   directory = directory.rstrip(os.path.sep)
   assert os.path.isdir(directory)
   num_sep = directory.count(os.path.sep)
   for root, dirs, files in os.walk(directory):
      yield root, dirs, files
      num_sep_this = root.count(os.path.sep)
      if num_sep + level <= num_sep_this:
         del dirs[:]

def xopen(filename, mode='r'):
   '''
   Replacement for the "open" function that can also open
   files that have been compressed with gzip. If the filename ends with .gz,
   the file is opened with gzip.open(). If it doesn't, the regular open()
   is used. If the filename is '-', standard output (mode 'w') or input
   (mode 'r') is returned.
   '''
   if sys.version_info[0] < 3:
      assert isinstance(filename, basestring)
   else:
      assert isinstance(filename, str)
   if filename == '-':
      return sys.stdin if 'rb' in mode else sys.stdout
   if filename.endswith('.gz'):
      if sys.version_info[0] < 3:
         return gzip.open(filename, mode)
      else:
         if 'rb' in mode:
            return getreader('ascii')(gzip.open(filename, mode))
         else:
            return getwriter('ascii')(gzip.open(filename, mode))
   else:
      return open(filename, mode)

class IterableQueue(SimpleQueue):
   _sentinel = object()
   def next(self):
      item = self.get()
      if item == "EOF":
         self.close()
         raise StopIteration
      return item
   def __iter__(self):
      return (iter(self.next, self._sentinel))
   def close(self):
      self.put(self._sentinel)


def grouper(n, iterable, padvalue=None):
   """grouper(3, 'abcdefg', 'x') -->
   ('a','b','c'), ('d','e','f'), ('g','x','x')"""
   return izip_longest(*[iter(iterable)]*n, fillvalue=padvalue)

def pileup_partial_split(pileup_line):
   '''
   split pileup line in chromosome, position and rest of line,
   without even strip
   '''
   chromosome, position, line = pileup_line.split('\t', 2)
   return [chromosome, int(position), line]

class multiPileups:
   '''
   Definig an iterable object merging two pilups
   '''
   def __init__(self, p1, p2):
      self.p1 = p1
      self.p2 = p2
      self._last_chromosome = None
   _sentinel        = object()
   def next(self):
      chromosome1, position1, line1 = pileup_partial_split(self.p1.next())
      chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
      if chromosome1 == chromosome2:
         self._last_chromosome = chromosome1
      going_on = True
      while going_on:
         '''
         If one of the files ever finish it would rise StopIteration
         and it would stop the iteration cleanly
         '''
         try:
            if chromosome1 == chromosome2:
               self._last_chromosome = chromosome1
               if position1 == position2:
                  ref1, line1 = line1.split('\t', 1)
                  ref2, line2 = line2.split('\t', 1)
                  return [chromosome1, position1, ref1, line1.strip(), line2.strip()]
                  going_on = False
               elif position1 > position2:
                  chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
               elif position1 < position2:
                  chromosome1, position1, line1 = pileup_partial_split(self.p1.next())
            else:
               if chromosome1 != self._last_chromosome:
                  chromosome2, position2, line2 = pileup_partial_split(self.p2.next())
               else:
                  chromosome1, position1, line1 = pileup_partial_split(self.p1.next())
         except StopIteration:
            going_on = False
            raise StopIteration
   '''
   def close(self):
      self.p1.close()
      self.p2.close()
   '''
   def __iter__(self):
      return (iter(self.next, self._sentinel))

def next_gcfile(self):
   gc_line = self.gc.next()
   try:
      pos, self._last_GC    = gc_line.split("\t", 1)
      self._last_window_s   = int(pos)
      self._last_window_e   = self._last_window_s + self._last_span
      return self
   except ValueError:
      line, chrom, span     = gc_line.strip().split()
      self._chromosome = chrom.split('chrom=')[1]
      self._last_span       = int(span.split('span=')[1])
      pos, self._last_GC    = self.gc.next().strip().split("\t", 1)
      self._last_window_s   = int(pos)
      self._last_window_e   = self._last_window_s + self._last_span
      return self
   except ValueError:
      if self._chromosome:
         raise StopIteration
      else:
         sys.exit("Error! The GC-content file is supported as the goldenpath hg19.gc5Base.txt.gz format.")

class GCmultiPileups:
   '''
   Add GC-content to the pilup, and still make an iterable object
   '''
   def __init__(self,multiPileups,gc):
      self.gc          = gc
      self             = next_gcfile(self)
      self.mpileup     = multiPileups
      self._last_chromosome = None
   _sentinel        = object()
   def next(self):
      self.pile_line   = self.mpileup.next()
      going_on = True
      while going_on:
         if self._chromosome == self.pile_line[0]:
            self._last_chromosome = self._chromosome
            if self.pile_line[1] >= self._last_window_s and self.pile_line[1] < self._last_window_e:
               self.pile_line.append(self._last_GC.strip())
               return self.pile_line
               going_on = False
            elif self.pile_line[1] < self._last_window_s:
               self.pile_line = self.mpileup.next()
            elif self.pile_line[1] >= self._last_window_e:
               self = next_gcfile(self)
         else:
            if self._last_chromosome != self._chromosome and self._last_chromosome != None:
               self.pile_line = self.mpileup.next()
            else:
               self = next_gcfile(self)
   '''
   def close(self):
      self.mpileup.close()
      self.gc.close()
   '''
   def __iter__(self):
      return (iter(self.next, self._sentinel))

def parse_pileup(ref_base, line, qlimit=20, qformat='sanger'):
   '''
   Parse the pileup format
   '''
   depth, mut_list, mut_qual = line
   ref_base = ref_base.upper()
   depth    = int(depth)
   if ref_base != "N":
      freq_dict = parse_pileup_seq(mut_list, mut_qual, depth, ref_base, qlimit, qformat)
      return [ref_base, depth, freq_dict['A'], freq_dict['C'], freq_dict['G'], freq_dict['T']]
   else:
      return [ref_base, depth, 0, 0, 0, 0]

def parse_pileup_str(line, min_depth, qlimit=20, qformat='sanger'):
   '''
   Parse the pileup format
   '''
   if line.strip():
      line   = line.strip()
      chromosome, n_base, ref_base, depth, mut_list, mut_qual = line.split()
      rd     = int(depth)
      rebase = ref_base.upper()
      if rd >= min_depth and rebase != "N":
         freq_dict = parse_pileup_seq(mut_list, mut_qual, rd, rebase, qlimit, qformat)
         line = [chromosome, n_base, rebase, depth, freq_dict['A'], freq_dict['C'], freq_dict['G'], freq_dict['T']]
         return '\t'.join(map(str,line))
      else:
         pass
   else:
      pass

def parse_pileup_seq(seq, quality, depth, reference, qlimit=20, qformat='sanger', noend=False, nostart=False):
   '''
   Parse the piled up sequence and return the frequence of mutation.
   '''
   nucleot_dict = {'A':0, 'C':0, 'G':0, 'T':0}
   reference    = reference.strip().upper()
   H            = 0
   n            = 0
   if qformat == 'sanger':
      qlimit = qlimit + 33
   elif qformat == 'illumina':
      qlimit = qlimit + 64
   quality      = quality.strip()
   block = {'seq':'','length':0}
   start     = False
   del_ins   = False
   l_del_ins = ''
   last_base = ''
   ins_del_length = 0
   for base in seq.strip():
      if block['length'] == 0:
         if base == '$':
            if noend:
               nucleot_dict[last_base] -= 1
         elif base == '^':
            start            = True
            block['length'] += 1
            block['seq']     = base
         elif base == '+' or base == '-':
            del_ins          = True
            block['length'] += 1
            block['seq']     = base
         elif base == '.' or base == ',':
            if ord(quality[n]) >= qlimit:
               nucleot_dict[reference] += 1
            last_base        = reference
            n += 1
         elif base.strip().upper() in nucleot_dict:
            if ord(quality[n]) >= qlimit:
               nucleot_dict[base.strip().upper()] += 1
            last_base        = base.strip().upper()
            n +=1
         else:
            n += 1
      else:
         if start:
            block['length'] += 1
            block['seq']    += base
            if block['length'] == 3:
               if not nostart:
                  if base == '.' or base == ',':
                     if ord(quality[n]) >= qlimit:
                        nucleot_dict[reference] += 1
                  elif base.strip().upper() in nucleot_dict:
                     if ord(quality[n]) >= qlimit:
                        nucleot_dict[base.strip().upper()] += 1
               block['length'] = 0
               block['seq']    = ''
               start           = False
               n += 1
         elif del_ins:
            if base.isdigit():
               l_del_ins       += base
               block['seq']    += base
               block['length'] += 1
            else:
               ins_del_length = int(l_del_ins)+ 1 + len(l_del_ins)
               block['seq']    += base
               block['length'] += 1
               if block['length'] == ins_del_length:
                  block['length'] = 0
                  block['seq']    = ''
                  l_del_ins       = ''
                  ins_del         = False
                  ins_del_length  = 0
      # Debug line
      #print str(n) + " " + base + " " + seq
   if n == depth:
      return nucleot_dict
   else:
      #print seq + " " + str(depth)+" " + str(n)
      sys.exit('Something went wrong on the block counting!!')


def line_worker(line, depth_sum, qlimit=20, qformat='sanger', hom_t=0.85, het_t=0.35):
   '''
   After the 3 files are syncronized we need to transform the pileup in a
   readable format, find the alleles, and compute the allele frequency in tumor
   '''
   if line:
      bases_list = ['A', 'C', 'G', 'T']
      chromosome, position, ref, p1_str, p2_str, gc = line
      p1_list = p1_str.split()
      p2_list = p2_str.split()
      if int(p1_list[0]) + int(p2_list[0]) >= depth_sum and len(p1_list) == len(p2_list):
         p1_mu  = parse_pileup(ref, p1_list, qlimit, qformat)
         sum_p1 = float(sum(p1_mu[2:6]))
         if sum_p1 > 0:
            p1_freq = [x/sum_p1 for x in p1_mu[2:6]]
            sort_freq_p1 = sorted(p1_freq, reverse = True)
            if sort_freq_p1[0] >= hom_t:
               # Homozygous positions here
               i        = p1_freq.index(sort_freq_p1[0])
               p2_mu    = parse_pileup(ref, p2_list, qlimit, qformat)
               sum_p2   = float(sum(p2_mu[2:6]))
               if sum_p2 > 0:
                  p2_freq  = [x/sum_p2 for x in p2_mu[2:6]]
                  homoz_p2 = p2_freq[i]
                  no_zero_idx   = [val for val in range(len(p2_freq)) if p2_freq[val] > 0]
                  no_zero_bases = [str(bases_list[ll])+str(round(p2_freq[ll], 3)) for ll in no_zero_idx if ll != i]
                  if no_zero_bases == []:
                     no_zero_bases = '.'
                  else:
                     no_zero_bases = ":".join(map(str,no_zero_bases))
                  #homoz_p2 = p2_mu[2 + i]/sum_p2
                  # chromosome, n_base, base_ref, depth_normal, depth_sample, depth.ratio, Af, Bf, ref.zygosity, GC-content, percentage reads above quality, AB.ref, AB.tum
                  line_out = [chromosome, position, p1_mu[0], p1_mu[1], p2_mu[1], round(p2_mu[1]/float(p1_mu[1]), 3), round(homoz_p2, 3), 0, 'hom', gc, int(sum_p2), bases_list[i], no_zero_bases]
                  return line_out
               else:
                  pass
            else:
               if sort_freq_p1[1] >= het_t:
                  # Heterozygous position here
                  allele = list()
                  for b in xrange(4):
                     if p1_freq[b] >= het_t:
                        allele.append(bases_list[b])
                  if len(allele) == 2:
                     p2_mu    = parse_pileup(ref, p2_list, qlimit, qformat)
                     sum_p2   = float(sum(p2_mu[2:6]))
                     if sum_p2 > 0:
                        i        = bases_list.index(allele[0])
                        ii       = bases_list.index(allele[1])
                        het_a_p2 = p2_mu[2 + i]/sum_p2
                        het_b_p2 = p2_mu[2 + ii]/sum_p2
                        if  het_a_p2 >= het_b_p2:
                           line_out = [chromosome, position, p1_mu[0], p1_mu[1], p2_mu[1], round(p2_mu[1]/float(p1_mu[1]), 3), round(het_a_p2 , 3), round(het_b_p2 , 3) , 'het', gc, int(sum_p2), bases_list[i]+bases_list[ii], "."]
                           return line_out
                        elif het_a_p2 < het_b_p2:
                           line_out = [chromosome, position, p1_mu[0], p1_mu[1], p2_mu[1], round(p2_mu[1]/float(p1_mu[1]), 3), round(het_b_p2 , 3), round(het_a_p2 , 3) , 'het', gc, int(sum_p2), bases_list[i]+bases_list[ii], "."]
                           return line_out
                  else:
                     pass
         else:
            pass
      else:
         pass
   else:
      pass

def stream_fasta(fa_file, window_size, out_queue):
   """
   Read and stream a fasta file in a Pipe
   """
   tmp_line   = ''
   for line in fa_file:
      line = line.strip()
      if line.startswith('>'):
         chromosome = line
         if tmp_line == '':
            #print chromosome
            out_queue.put(chromosome)
         else:
            groups = grouper(window_size, tmp_line.upper())
            for group in groups:
               #print group
               out_queue.put(group)
            tmp_line = ''
            #print chromosome
            out_queue.put(chromosome)
      else:
         tmp_line = tmp_line + line
         if len(tmp_line)%window_size == 0:
            groups = grouper(window_size, tmp_line.upper())
            for group in groups:
               #print group
               out_queue.put(group)
            tmp_line = ''
         else:
            #print tmp_line
            pass
   if tmp_line == '':
      pass
   else:
      groups = grouper(window_size, tmp_line.upper())
      for group in groups:
         #print group
         out_queue.put(group)
   out_queue.put('EOF')

def seq_map(seq_list):
   """
   ...
   """
   counter = {"A":0,"C":0,"G":0,"N":0,"T":0}
   for nucleotide in seq_list:
      counter[nucleotide] +=1
   return counter

def process_gc_from_pipe(in_queue, window_size):
   """
   whatever...
   """
   base_counter  = 1
   window_size   = float(window_size)
   while True:
      seq_list   = in_queue.get()
      #print seq_list
      if seq_list != 'EOF':
         if type(seq_list) == str:
            if seq_list.startswith('>'):
               chromosome = seq_list.strip('>')
               print "variableStep chrom=" + chromosome + " span=" + str(int(window_size))
               base_counter = 1
         else:
            seq_list = filter(None,seq_list)
            stats    = seq_map(seq_list)
            seq_len  = len(seq_list)
            if seq_len == window_size:
               if stats['N']/window_size >= 0.75:
                  base_counter = base_counter + window_size
                  pass
               else:
                  gc_percent = 100 * (stats['G'] + stats['C'])/window_size
                  print str(int(base_counter)) + "\t" + str(gc_percent)
                  base_counter = base_counter + window_size
            else:
               if stats['N'] / seq_len >= 0.75:
                  base_counter = base_counter + seq_len
                  pass
               else:
                  gc_percent   = 100 * (stats['G'] + stats['C'])/seq_len
                  print str(int(base_counter)) + "\t" + str(gc_percent)
                  base_counter = base_counter + seq_len
      else:
         break

class SubcommandHelpFormatter(argparse.RawDescriptionHelpFormatter):
    def _format_action(self, action):
        parts = super(argparse.RawDescriptionHelpFormatter, self)._format_action(action)
        if action.nargs == argparse.PARSER:
            parts = "\n".join(parts.split("\n")[1:])
        return parts

class DefaultHelpParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

def RPy2sqeezeABfreq(abfreq, loop, tag, out, kmin, gamma, mufreq_treshold):
   """
    Process an abfreq file and store the relvant (small and easy to re-acces) information
   """
   is_gz = False
   if os.path.split(abfreq)[1][-2:] == "gz":
      is_gz = True
   if loop:
      print("loading GC-content into memory")
      gc_stats = sequenza.gc_sample_stats(abfreq)
      chr_vect = robjects.r('as.character')(gc_stats.rx2('file.metrics').rx2('chr'))
      gc_vect  = robjects.r.setNames(gc_stats.rx2('raw.mean'), gc_stats.rx2('gc.values'))
   else:
      chr_vect = list()
      chr_vect.append(1)
   windows_baf   = robjects.ListVector({})
   windows_ratio = robjects.ListVector({})
   mutation_list = robjects.ListVector({})
   segments_list = robjects.ListVector({})
   for chr in chr_vect:
      if loop:
         file_lines = gc_stats.rx2('file.metrics').rx(chr, True)
         abf_data = sequenza.read_abfreq(abfreq, gz = is_gz, n_lines = robjects.IntVector((file_lines.rx('start')[0][0],  file_lines.rx('end')[0][0])))
      else:
         print("loading all the file in memory")
         abf_data = sequenza.read_abfreq(abfreq, fast = True, gz = is_gz)
         gc_stats = sequenza.gc_norm(x = abf_data.rx(True, 'depth.ratio'), gc = abf_data.rx(True, 'GC.percent'))
         gc_vect  = robjects.r.setNames(gc_stats.rx2('raw').rx(True, '50%'), gc_stats.rx2('gc.values'))
         chr_vect = robjects.r.unique(abf_data.rx(True, 'chromosome'))
         chr_vect = robjects.r('as.character')(chr_vect)
      abf_data = robjects.r.cbind(abf_data, adjusted_ratio = robjects.r.round(abf_data.rx(True, 'depth.ratio').ro / gc_vect.rx(robjects.r('as.character')(abf_data.rx(True, 'GC.percent'))), 3))
      abf_data.names[-1] = 'adjusted.ratio'
      abf_hom  = abf_data.rx(True, 'ref.zygosity').ro == 'hom'
      abf_het  = abf_data.rx(abf_hom.ro != True, True)
      abf_r_win = sequenza.windowValues(x = abf_data.rx(True, 'adjusted.ratio'), positions = abf_data.rx(True, 'n.base'), chromosomes = abf_data.rx(True, 'chromosome'), window = 1e6, overlap = 1, weight = abf_data.rx(True, 'depth.normal'))
      abf_b_win = sequenza.windowValues(x = abf_het.rx(True, 'Bf'), positions = abf_het.rx(True, 'n.base'), chromosomes = abf_het.rx(True, 'chromosome'), window = 1e6, overlap = 1, weight = robjects.r.round(abf_het.rx(True, 'good.s.reads'), 0))
      breaks = sequenza.find_breaks(abf_het, gamma = gamma, kmin = kmin, baf_thres = robjects.FloatVector((0, 0.5)))
      seg_s1 = sequenza.segment_breaks(abf_data, breaks = breaks)
      mut_tab = sequenza.mutation_table(abf_data, mufreq_treshold = mufreq_treshold, min_reads = 40, max_mut_types = 1, min_type_freq = 0.9, segments = seg_s1)
      gc.collect()
      if loop:
         windows_baf.rx2[chr]   = abf_b_win.rx2(1)
         windows_ratio.rx2[chr] = abf_r_win.rx2(1)
         mutation_list.rx2[chr] = mut_tab
         segments_list.rx2[chr] = seg_s1
      else:
         windows_baf   = abf_b_win
         windows_ratio = abf_r_win
         for chr in chr_vect:
            mutation_list.rx2[chr] = mut_tab.rx(mut_tab.rx(True, 'chromosome').ro == chr, True)
            segments_list.rx2[chr] = seg_s1.rx(seg_s1.rx(True, 'chromosome').ro == chr, True)
   windows_baf.names   = chr_vect
   windows_ratio.names = chr_vect
   mutation_list.names = chr_vect
   segments_list.names = chr_vect
   subdir = out +'/' + tag
   if not os.path.exists(subdir):
      os.makedirs(subdir)
   res = robjects.ListVector({'BAF' : windows_baf, 'ratio' : windows_ratio, 'mutations' : mutation_list,
               'segments' : segments_list, 'chromosomes' : chr_vect, 'gc' : gc_stats})
   robjects.r('assign')(x = tag + '_sequenza_extract', value = res)
   robjects.r('save')(list = tag + '_sequenza_extract', file = subdir + '/' + tag + '_sequenza_extract.Rdata')
   # robjects.r('save')(list = tag + '_windows_ratio', file = subdir + '/' + tag + '_windows_ratio.Rdata')
   # robjects.r('write.table')(x = gc_stats.rx2('raw'), file = subdir +'/' + tag + '_raw_GC.txt', col_names = True, row_names = False, sep = "\t")
   # robjects.r('write.table')(x = gc_stats.rx2('adj'), file = subdir +'/' + tag + '_adj_GC.txt', col_names = True, row_names = False, sep = "\t")
   # robjects.r('assign')(x = tag + '_windows_Bf', value = windows_baf)
   # robjects.r('save')(list = tag + '_windows_Bf', file = subdir + '/' + tag + '_windows_Bf.Rdata')
   # robjects.r('assign')(x = tag + '_windows_ratio', value = windows_ratio)
   # robjects.r('save')(list = tag + '_windows_ratio', file = subdir + '/' + tag + '_windows_ratio.Rdata')
   # robjects.r('assign')(x = tag + '_mutation_list', value = mutation_list)
   # robjects.r('save')(list = tag + '_mutation_list', file = subdir + '/' + tag + '_mutation_list.Rdata')
   # robjects.r('assign')(x = tag + '_segments_list', value = segments_list)
   # robjects.r('save')(list = tag + '_segments_list', file = subdir + '/' + tag + '_segments_list.Rdata')

def RPy2doAllSequenza(data_dir, is_male = True, tag = None, X = "X", Y = "Y", ncores = 4, ratio_priority = False, segment_filter = 10e6, priors = {'CN':[1,2,3,4], 'value' : [1,2,1,1]}):
   '''
   Load the information stored from the functions above and infer cellularity/ploidy
   and save plots and table of mutations and CNV calls
   '''
   xy = {'X':X, 'Y' : Y}
   if tag == None:
      tag = os.path.split(data_dir)[-1]
   robjects.r.load(data_dir +'/' + tag + '_sequenza_extract.Rdata')
   extract = robjects.r.eval(robjects.r('as.name')(tag + '_sequenza_extract'))
   windows_Bf    = extract.rx2('BAF')
   windows_ratio = extract.rx2('ratio')
   mutation_list = extract.rx2('mutations')
   segments_list = extract.rx2('segments')
   chr_vect      = extract.rx2('chromosomes')
   # obj_list = robjects.StrVector(('adj_GC.txt', 'raw_GC.txt', 'windows_Bf.Rdata', 'mutation_list.Rdata', 'windows_ratio.Rdata', 'segments_list.Rdata'))
   # obj_list = robjects.r.paste(tag, obj_list, sep = "_")
   # gc_tab   = robjects.r('read.table')(data_dir +'/' + obj_list[0], header = True, sep = '\t')
   # for i in range(2,6):
   #    robjects.r.load(data_dir +'/' + obj_list[i])
   # windows_Bf    = robjects.r.eval(robjects.r('as.name')(tag + '_windows_Bf'))
   # windows_ratio = robjects.r.eval(robjects.r('as.name')(tag + '_windows_ratio'))
   # mutation_list = robjects.r.eval(robjects.r('as.name')(tag + '_mutation_list'))
   # segments_list = robjects.r.eval(robjects.r('as.name')(tag + '_segments_list'))
   # chr_vect      = windows_Bf.names
   avg_depth_ratio = robjects.r.mean(extract.rx2('gc').rx2('adj').rx(True, 2))
   segs_all      = robjects.r('do.call')('rbind', segments_list)
   mut_all       = robjects.r('do.call')('rbind', mutation_list)
   mut_all       = robjects.r('na.exclude')(mut_all)
   segs_len      = segs_all.rx(True, 'end.pos').ro - segs_all.rx(True, 'start.pos')
   segs_filt     = segs_len.ro >= segment_filter
   if is_male:
      segs_is_xy = segs_all.rx(True, 'chromosome').ro == xy["X"] or segs_all.rx(True, 'chromosome').ro == xy["Y"]
      mut_is_xy  = mut_all.rx(True, 'chromosome').ro == xy["X"] or mut_all.rx(True, 'chromosome').ro == xy["Y"]
   else:
      segs_is_xy = segs_all.rx(True, 'chromosome').ro == xy["Y"]
   filt_test  = segs_is_xy.ro == False
   filt_test  = segs_filt.ro & filt_test
   seg_test   = segs_all.rx(filt_test, True)
   weights_seg = robjects.r.round(segs_len.ro / 1e6, 0).ro + 150
   robjects.r('''
   wrapBafBayes <- function (Bf, depth_ratio , weight_ratio, weight_Bf,
                             avg_depth_ratio,  cellularity, priors_table,
                             ploidy, mc_cores, ratio_priority, ...) {
                    baf.model.fit(Bf = Bf, depth.ratio = depth_ratio,
                    weight.ratio = weight_ratio, weight.Bf = weight_Bf,
                    avg.depth.ratio = avg_depth_ratio, cellularity = cellularity,
                    ploidy = ploidy, priors.table = priors_table,
                    mc.cores = mc_cores, ratio.priority = ratio_priority)
   }
   ''')
   CP  = robjects.r('wrapBafBayes')(Bf = seg_test.rx(True, 'Bf'), depth_ratio = seg_test.rx(True, 'depth.ratio'),
                    weight_ratio = weights_seg.rx(filt_test).ro * 2,
                    weight_Bf = weights_seg.rx(filt_test), avg_depth_ratio = avg_depth_ratio,
                    cellularity = robjects.r.seq(0.1, 1, 0.01) , priors_table = robjects.DataFrame({'CN':robjects.IntVector(priors["CN"]), 'value' : robjects.IntVector(priors["value"])}),
                    ploidy = robjects.r.seq(1, 7, 0.1), mc_cores = ncores, ratio_priority = ratio_priority)

   cint = sequenza.get_ci(CP)

   robjects.r.pdf(data_dir +'/'+ tag + "_CP_ci.pdf")
   robjects.r.par(mfrow = robjects.IntVector((2, 2)))
   sequenza.cp_plot(CP)
   robjects.r.plot(cint.rx2('values.y'), ylab = "Cellularity", xlab = "Likelihood", type = "l")
   robjects.r.abline(h = cint.rx2('confint.y'), lty = 2, lwd = 0.5, col = "red")
   robjects.r.plot(cint.rx2('values.x'), ylab = "Ploidy", xlab = "Likelihood", type = "l")
   robjects.r.abline(v = cint.rx2('confint.x'), lty = 2, lwd = 0.5, col = "red")
   robjects.r('dev.off()')
   robjects.r.pdf(data_dir +'/'+ tag + "_CP_contours.pdf")
   sequenza.cp_plot(CP)
   sequenza.cp_plot_contours(CP, add = True, likThresh = 0.95)
   robjects.r('dev.off()')

   if is_male:
      seg_res  = sequenza.baf_bayes(Bf = segs_all.rx(True, 'Bf').rx(segs_is_xy.ro == False),
                         depth_ratio = segs_all.rx(True, 'depth.ratio').rx(segs_is_xy.ro == False),
                         avg_depth_ratio = avg_depth_ratio,
                         weight_ratio = 2*200, ratio_priority = ratio_priority,
                         weight_Bf = 200, CNt_max = 20,
                         cellularity = cint.rx2('max.y'),
                         ploidy = cint.rx2('max.x'), CNn = 2)
      mut_alleles = sequenza.mufreq_bayes(mufreq = mut_all.rx(True, 'F').rx(mut_is_xy.ro == False),
                                          depth_ratio = mut_all.rx(True, 'adjusted.ratio').rx(mut_is_xy.ro == False),
                                          cellularity = cint.rx2('max.y'), ploidy = cint.rx2('max.x'),
                                          avg_depth_ratio = avg_depth_ratio, CNn = 2, CNt_max = 20)
      seg_res_xy  = sequenza.baf_bayes(Bf = segs_all.rx(True, 'Bf').rx(segs_is_xy),
                         depth_ratio = segs_all.rx(True, 'depth.ratio').rx(segs_is_xy),
                         avg_depth_ratio = avg_depth_ratio,
                         weight_ratio = 2*200, ratio_priority = ratio_priority,
                         weight_Bf = 200, CNt_max = 20,
                         cellularity = cint.rx2('max.y'),
                         ploidy = cint.rx2('max.x'), CNn = 1)
      mut_alleles_xy = sequenza.mufreq_bayes(mufreq = mut_all.rx(True, 'F').rx(mut_is_xy),
                                             depth_ratio = mut_all.rx(True, 'adjusted.ratio').rx(mut_is_xy),
                                             cellularity = cint.rx2('max.y'), ploidy = cint.rx2('max.x'),
                                             avg_depth_ratio = avg_depth_ratio, CNn = 1, CNt_max = 20)
      seg_res = robjects.r.cbind(robjects.r.rbind(segs_all.rx(segs_is_xy.ro == False,True),segs_all.rx(segs_is_xy, True)), robjects.r.rbind(seg_res, seg_res_xy))
      mut_res = robjects.r.cbind(robjects.r.rbind(mut_all.rx(mut_is_xy.ro == False,True), mut_all.rx(mut_is_xy, True)), robjects.r.rbind(mut_alleles, mut_alleles_xy))
   else:
      seg_res  = sequenza.baf_bayes(Bf = segs_all.rx(True, 'Bf'),
                         depth_ratio = segs_all.rx(True, 'depth.ratio'),
                         avg_depth_ratio = avg_depth_ratio,
                         weight_ratio = 2*200, ratio_priority = ratio_priority,
                         weight_Bf = 200, CNt_max = 20,
                         cellularity = cint.rx2('max.y'),
                         ploidy = cint.rx2('max.x'), CNn = 2)
      mut_alleles = sequenza.mufreq_bayes(mufreq = mut_all.rx(True, 'F'),
                                          depth_ratio = mut_all.rx(True, 'adjusted.ratio'),
                                          cellularity = cint.rx2('max.y'), ploidy = cint.rx2('max.x'),
                                          avg_depth_ratio = avg_depth_ratio, CNn = 2, CNt_max = 20)
      seg_res = robjects.r.cbind(segs_all, seg_res)
      mut_res = robjects.r.cbind(mut_all, mut_alleles)
   robjects.r('write.table')(seg_res, data_dir +'/'+ tag + "_segments.txt", col_names = True, row_names = False, sep = "\t")
   robjects.r('write.table')(mut_res, data_dir +'/'+ tag + "_mutations.txt", col_names = True, row_names = False, sep = "\t")


   robjects.r.pdf(data_dir +'/'+ tag + "_chromosome_view.pdf")
   for chrom in chr_vect:
      if is_male:
         if chrom == xy["X"] or chrom == xy["Y"]:
            CNn = 1
         else:
            CNn = 2
      else:
         CNn = 2
      sequenza.chromosome_view(baf_windows = windows_Bf.rx2(chrom),
                      ratio_windows = windows_ratio.rx2(chrom), min_N_ratio = 1,
                      cellularity = cint.rx('max.y')[0][0], ploidy = cint.rx('max.x')[0][0],
                      segments = seg_res.rx(seg_res.rx(True, 'chromosome').ro == chrom, True), mut_tab = mutation_list.rx2(chrom),
                      main = chrom, avg_depth_ratio = avg_depth_ratio, CNn = CNn, BAF_style = "lines")
   robjects.r('dev.off()')
   res_seg_xy = seg_res.rx(True, 'chromosome').ro == xy["Y"]


   barscn = robjects.DataFrame({'size' : (seg_res.rx(True, 'end.pos').rx(res_seg_xy.ro == False).ro - seg_res.rx(True, 'start.pos').rx(res_seg_xy.ro == False)),
                                'CNt' : seg_res.rx(True, 'CNt').rx(res_seg_xy.ro == False)})
   cn_sizes = robjects.r.split(barscn.rx(True,'size'),barscn.rx(True,'CNt'))
   cn_sizes = robjects.r.sapply(cn_sizes, 'sum')
   robjects.r.pdf(data_dir +'/'+ tag + "_CN_bars.pdf")
   robjects.r.barplot(robjects.r.round((cn_sizes.ro/robjects.r.sum(cn_sizes)).ro * 100, 0), names = cn_sizes.names, las = 1,
           ylab = "Percentage (%)", xlab = "copy number")
   robjects.r('dev.off()')
   res_ci_tab = robjects.DataFrame({'cellularity' : robjects.FloatVector((cint.rx2('confint.y')[0], cint.rx2('max.y')[0], cint.rx2('confint.y')[1])),
                                     'ploidy' : robjects.FloatVector((cint.rx2('confint.x')[0], cint.rx2('max.x')[0], cint.rx2('confint.x')[1])),
                                     'ploidy.mean'      : robjects.r('weighted.mean')(x=robjects.r('as.integer')(cn_sizes.names), w = cn_sizes)})
   robjects.r('write.table')(res_ci_tab, data_dir +'/'+ tag + "_confints_CP.txt", col_names = True, row_names = False, sep = "\t")

def RPy2SequenzaOverride(data_dir, is_male = True, tag = None, X = "X", Y = "Y", ratio_priority = False, cellularity = None, ploidy = None):
   '''
   Fit the data with arbitrary values of cellularity and ploidy.
   Useful when cellularity is so low that can't be estimate correctly.
   '''
   xy = {'X':X, 'Y' : Y}
   if tag == None:
      tag = os.path.split(data_dir)[-1]
   robjects.r.load(data_dir +'/' + tag + '_sequenza_extract.Rdata')
   extract = robjects.r.eval(robjects.r('as.name')(tag + '_sequenza_extract'))
   windows_Bf    = extract.rx2('BAF')
   windows_ratio = extract.rx2('ratio')
   mutation_list = extract.rx2('mutations')
   segments_list = extract.rx2('segments')
   chr_vect      = extract.rx2('chromosomes')
   #obj_list = robjects.StrVector(('adj_GC.txt', 'raw_GC.txt', 'windows_Bf.Rdata', 'mutation_list.Rdata', 'windows_ratio.Rdata', 'segments_list.Rdata'))
   #obj_list = robjects.r.paste(tag, obj_list, sep = "_")
   #gc_tab   = robjects.r('read.table')(data_dir +'/' + obj_list[0], header = True, sep = '\t')
   #avg_depth_ratio = robjects.r.mean(gc_tab.rx(True, 2))
   #for i in range(2,6):
   #   robjects.r.load(data_dir +'/' + obj_list[i])
   #windows_Bf    = robjects.r.eval(robjects.r('as.name')(tag + '_windows_Bf'))
   #windows_ratio = robjects.r.eval(robjects.r('as.name')(tag + '_windows_ratio'))
   #mutation_list = robjects.r.eval(robjects.r('as.name')(tag + '_mutation_list'))
   #segments_list = robjects.r.eval(robjects.r('as.name')(tag + '_segments_list'))
   #chr_vect      = windows_Bf.names
   avg_depth_ratio = robjects.r.mean(extract.rx2('gc').rx2('adj').rx(True, 2))
   segs_all      = robjects.r('do.call')('rbind', segments_list)
   mut_all       = robjects.r('do.call')('rbind', mutation_list)
   mut_all       = robjects.r('na.exclude')(mut_all)
   if is_male:
      segs_is_xy = segs_all.rx(True, 'chromosome').ro == xy["X"] or segs_all.rx(True, 'chromosome').ro == xy["Y"]
      mut_is_xy  = mut_all.rx(True, 'chromosome').ro == xy["X"] or mut_all.rx(True, 'chromosome').ro == xy["Y"]
   else:
      segs_is_xy = segs_all.rx(True, 'chromosome').ro == xy["Y"]
   if is_male:
      seg_res  = sequenza.baf_bayes(Bf = segs_all.rx(True, 'Bf').rx(segs_is_xy.ro == False),
                         depth_ratio = segs_all.rx(True, 'depth.ratio').rx(segs_is_xy.ro == False),
                         avg_depth_ratio = avg_depth_ratio,
                         weight_ratio = 2*200, ratio_priority = ratio_priority,
                         weight_Bf = 200, CNt_max = 20,
                         cellularity = cellularity,
                         ploidy = ploidy, CNn = 2)
      mut_alleles = sequenza.mufreq_bayes(mufreq = mut_all.rx(True, 'F').rx(mut_is_xy.ro == False),
                                          depth_ratio = mut_all.rx(True, 'adjusted.ratio').rx(mut_is_xy.ro == False),
                                          cellularity = cellularity, ploidy = ploidy,
                                          avg_depth_ratio = avg_depth_ratio, CNn = 2, CNt_max = 20)
      seg_res_xy  = sequenza.baf_bayes(Bf = segs_all.rx(True, 'Bf').rx(segs_is_xy),
                         depth_ratio = segs_all.rx(True, 'depth.ratio').rx(segs_is_xy),
                         avg_depth_ratio = avg_depth_ratio,
                         weight_ratio = 2*200, ratio_priority = ratio_priority,
                         weight_Bf = 200, CNt_max = 20,
                         cellularity = cellularity,
                         ploidy = ploidy, CNn = 1)
      mut_alleles_xy = sequenza.mufreq_bayes(mufreq = mut_all.rx(True, 'F').rx(mut_is_xy),
                                             depth_ratio = mut_all.rx(True, 'adjusted.ratio').rx(mut_is_xy),
                                             cellularity = cellularity, ploidy = ploidy,
                                             avg_depth_ratio = avg_depth_ratio, CNn = 1, CNt_max = 20)
      seg_res = robjects.r.cbind(robjects.r.rbind(segs_all.rx(segs_is_xy.ro == False,True),segs_all.rx(segs_is_xy, True)), robjects.r.rbind(seg_res, seg_res_xy))
      mut_res = robjects.r.cbind(robjects.r.rbind(mut_all.rx(mut_is_xy.ro == False,True), mut_all.rx(mut_is_xy, True)), robjects.r.rbind(mut_alleles, mut_alleles_xy))
   else:
      seg_res  = sequenza.baf_bayes(Bf = segs_all.rx(True, 'Bf'),
                         depth_ratio = segs_all.rx(True, 'depth.ratio'),
                         avg_depth_ratio = avg_depth_ratio,
                         weight_ratio = 2*200, ratio_priority = ratio_priority,
                         weight_Bf = 200, CNt_max = 20,
                         cellularity = cellularity,
                         ploidy = ploidy, CNn = 2)
      mut_alleles = sequenza.mufreq_bayes(mufreq = mut_all.rx(True, 'F'),
                                          depth_ratio = mut_all.rx(True, 'adjusted.ratio'),
                                          cellularity = cellularity, ploidy = ploidy,
                                          avg_depth_ratio = avg_depth_ratio, CNn = 2, CNt_max = 20)
      seg_res = robjects.r.cbind(segs_all, seg_res)
      mut_res = robjects.r.cbind(mut_all, mut_alleles)
   override_prefix = "_C"+ str(cellularity)+ "_P" + str(ploidy)
   robjects.r('write.table')(seg_res, data_dir +'/'+ tag + override_prefix + "_segments.txt", col_names = True, row_names = False, sep = "\t")
   robjects.r('write.table')(mut_res, data_dir +'/'+ tag + override_prefix + "_mutations.txt", col_names = True, row_names = False, sep = "\t")
   robjects.r.pdf(data_dir +'/'+ tag + override_prefix + "_chromosome_view.pdf")
   for chrom in chr_vect:
      if is_male:
         if chrom == xy["X"] or chrom == xy["Y"]:
            CNn = 1
         else:
            CNn = 2
      else:
         CNn = 2
      sequenza.chromosome_view(baf_windows = windows_Bf.rx2(chrom),
                      ratio_windows = windows_ratio.rx2(chrom), min_N_ratio = 1,
                      cellularity = cellularity, ploidy = ploidy,
                      segments = seg_res.rx(seg_res.rx(True, 'chromosome').ro == chrom, True), mut_tab = mutation_list.rx2(chrom),
                      main = chrom, avg_depth_ratio = avg_depth_ratio, CNn = CNn, BAF_style = "lines")
   robjects.r('dev.off()')
   res_seg_xy = seg_res.rx(True, 'chromosome').ro == xy["Y"]

   barscn = robjects.DataFrame({'size' : (seg_res.rx(True, 'end.pos').rx(res_seg_xy.ro == False).ro - seg_res.rx(True, 'start.pos').rx(res_seg_xy.ro == False)),
                                'CNt' : seg_res.rx(True, 'CNt').rx(res_seg_xy.ro == False)})
   cn_sizes = robjects.r.split(barscn.rx(True,'size'),barscn.rx(True,'CNt'))
   cn_sizes = robjects.r.sapply(cn_sizes, 'sum')
   robjects.r.pdf(data_dir +'/'+ tag + override_prefix + "_CN_bars.pdf")
   robjects.r.barplot(robjects.r.round((cn_sizes.ro/robjects.r.sum(cn_sizes)).ro * 100, 0), names = cn_sizes.names, las = 1,
           ylab = "Percentage (%)", xlab = "copy number")
   robjects.r('dev.off()')


def pileup2acgt(parser, subparser):
   subparser.add_argument('pileup', metavar='pileup',
                   help='pileup (SAMtools) file in input, if a gzipped file will be selected it will be opened in gzip mode, if file name is - it would be loaded from STDIN.')
   subparser.add_argument('-n', dest='n', type=int,
                   help='The minimum read depth on the base to consider the mutation on it.')
   parser_pup2muoutput      = subparser.add_argument_group(title='Output', description='Argument that involve the output destination.')
   parser_pup2muoutput.add_argument('-o', '--output', dest='output', type = str, default = '-',
                       help='Destination of the output file. To use gzip compression name the file ending by .gz. (default STDOUT).')
   parser_pup2muoutput.add_argument('--quiet', dest='quiet', action="store_true",
                       help='Do not output additional debugging.')
   parser_pup2muperformance = subparser.add_argument_group(title='Performance', description='Argument that can effect the performance.')
   parser_pup2muperformance.add_argument('-p', '--processes', dest='nproc', default="0", type=int,
                   help='Set the number of processes to split the parsing job. If it is set to 0 (default), the job will occurs with no forking to other processes. If it is bigger then 0 it is more efficient with an adequate chunk size, otherwise with smaller chuncks (eg.: < 1000) it will loose performance.')
   parser_pup2muperformance.add_argument('-c', '--chunk', dest='chunk', default="0", type=int,
                   help='Set the size (number of lines) of the portion of the file to assign to each process. If is set to 0 (default) wil set to default also the --processes parameter (-p 0). An adequate chunk size defends on the number of processes and on the file size (chunk size bigger then total number of line is not good). However a chunk size ~ 1000 leads to better performance.')
   parser_pup2muqualitysets = subparser.add_argument_group(title='Quality and Format', description='Argument that change the quality threshold or the quality format.')
   parser_pup2muqualitysets.add_argument('-q', '--qlimit', dest='qlimit', default=20,type=int,
                   help='Minimum nucleotide quality score for consider in the counts.')
   parser_pup2muqualitysets.add_argument('-f', '--qformat', dest='qformat', default="sanger",
                   help='Quality format, options are sanger or illumina, it will add an offset of 33 or 64 respectively to the qlimit value.')
   parser_pup2muinfo        = subparser.add_argument_group(title='Info',description='Other Information.')
   parser_pup2muinfo.add_argument('-v', '--version', action="version", version='%(prog)s version: ' + VERSION + ". " + DATE,
                   help='Just display the version information and exit.')
   return parser.parse_args()

def pileup2abfreq(parser, subparser):
   parser_ABinput    = subparser.add_argument_group(title='Input Files',description='Required input files.')
   parser_ABgenotype    = subparser.add_argument_group(title='Genotyper',description='Options regarding the genotyping.')
   parser_ABperformance = subparser.add_argument_group(title='Performance', description='Options affecting the performance.')
   parser_ABqualitysets = subparser.add_argument_group(title='Quality and Format', description='Options that change the quality threshold and format.')
   parser_ABinput.add_argument('-r', '--reference', dest = 'reference', required = True,
                   help='The pileup of the reference/normal sample')
   parser_ABinput.add_argument('-s', '--sample', dest = 'sample', required = True,
                   help='The pileup of the tumor sample')
   parser_ABinput.add_argument('-gc', dest = 'gc', metavar = 'gc', required = True,
                   help='The GC-content file coming from UCSC genome browser, or generated in the same UCSC format')
   parser_ABqualitysets.add_argument('-q', '--qlimit', dest = 'qlimit', default = 20, type = int,
                   help='Minimum nucleotide quality score for consider in the counts. Default 20.')
   parser_ABqualitysets.add_argument('-f', '--qformat', dest = 'qformat', default = "sanger",
                   help='Quality format, options are sanger or illumina, it will add an offset of 33 or 64 respectively to the qlimit value. Default "sanger".')
   parser_ABqualitysets.add_argument('-n', dest = 'n', type = int, default = 20,
                   help='Threshold to filter positions by the sum of read depth of the two samples. Default 20.')
   parser_ABgenotype.add_argument('--hom', dest = 'hom', type = float, default = 0.9,
                   help='Threshold to select homozygous positions. Default 0.9.')
   parser_ABgenotype.add_argument('--het', dest = 'het', type = float, default = 0.25,
                   help='Threshold to select heterozygous positions. Default 0.25.')
   parser_ABperformance.add_argument('-p', '--processes', dest='nproc', default="0", type=int,
                   help='Set the number of processes to split the genotyping. If it is set to 0 (default), the job will occurs with no forking to other processes. If it is bigger then 0 it is more efficient with an adequate chunk size, otherwise with smaller chunks (eg.: < 1000) it will loose performance. Default 0')
   parser_ABperformance.add_argument('-c', '--chunk', dest='chunk', default="1", type=int,
                   help='Set the number of lines to assign to each process. If is set to 1 (default) will set to default also the --processes parameter (-p 0). An adequate chunk size defends on the number of processes and on the file size (chunk size bigger then total number of line is not good). However a chunk size ~ 1000 leads to better performance. Default 1.')
   return parser.parse_args()

def GC_windows(parser, subparser):
   subparser.add_argument('fasta', metavar = 'fasta',
                   help='the fasta file. It can be a file name or \"-\" to take the input from STDIN')
   subparser.add_argument('-w', dest = 'window', type = int, default = 50,
                   help='The window size to calculate the GC-content percentage')
   return parser.parse_args()

def merge_pileups(parser, subparser):
   subparser.add_argument('-1', '--pileup1', dest = 'p1', required = True,
                   help='The first pileup')
   subparser.add_argument('-2', '--pileup2', dest = 'p2', required = True,
                   help='The second pileup, will show as the last columns set')
   return parser.parse_args()

def sequenzaExtract(parser, subparser):
   parser_io      = subparser.add_argument_group(title='Input and output',description='Input ABfreq files and output options.')
   parser_segment = subparser.add_argument_group(title='Segmentation',description='Option to control the segmentation.')
   parser_mut     = subparser.add_argument_group(title='Mutations',description='Option to filter the mutations by variant allele frequency.')
   parser_misc    = subparser.add_argument_group(title='Misc',description='Miscellaneous options.')
   parser_io.add_argument('--abfreq', dest = 'abfreq', required = True,
                   help='An existing abfreq file')
   parser_io.add_argument('-o', '--out', dest = 'dir', default = "./",
                   help='Path where to make the directory containing the processed data')
   parser_io.add_argument('-t', '--tag', dest = 'tag', required = True,
                   help='Tag to name the directory and the prefix of the generated files')
   parser_segment.add_argument('-k', '--kmin', dest = 'kmin', type = int, default = 10,
                   help='minimum number of position per segment. default 10 (WGS is suggested to set to 500 or so)')
   parser_segment.add_argument('-g', '--gamma', dest = 'gamma', type = int, default = 80,
                   help='gamma parameter for the segmentation, higher is less sensible smaller is more. default 80')
   parser_mut.add_argument('-f', '--mut-threshold', dest = 'mufreq', type = float, default = 0.1,
                   help='Threshold on the variant allele frequency to filter out the mutations. Default 0.1.')
   parser_misc.add_argument('--no-loop', dest = 'loop', action='store_false', default = True,
                   help='Boolean flag indicating if to loop over chromosomes one by one (default), or load all the file in memory')
   return parser.parse_args()

def sequenzaFit(parser, subparser):
   parser_io      = subparser.add_argument_group(title='Input',description='Input options.')
   parser_gender  = subparser.add_argument_group(title='Gender',description='Option to control the gender and X/Y chromosome.')
   parser_model   = subparser.add_argument_group(title='Model',description='Options to control the Bayesian inference.')
   parser_perf    = subparser.add_argument_group(title='Performance',description='Options to control performance.')
   parser_io.add_argument('--dir', dest = 'dir', required = True,
                   help='The directory where the data to load are stored')
   parser_io.add_argument('-t', '--tag', dest = 'tag', default = None,
                   help='Tag indicating the prefix of the data, if not specified it is assumed as the name of the specified directory')
   parser_gender.add_argument('--is-male', dest = 'isMale', action='store_true', default = False,
                   help='Boolean flag indicating if the sequencing data are from a male or female, and consequently properly handle chromosome X and Y')
   parser_perf.add_argument('-p', dest = 'ncpu', type = int, default = 4,
                   help='The number of core to use when performing the Bayesian inference. Default 4.')
   parser_gender.add_argument('-X', "--chrX", dest = 'X', type = str, default = "X",
                   help='Character defining chromosome X. Default X.')
   parser_gender.add_argument('-Y', "--chrY", dest = 'Y', type = str, default = "Y",
                   help='Character defining chromosome Y. Default Y.')
   parser_model.add_argument('-r', "--only-ratio", dest = 'onlyratio', action='store_true', default = False,
                   help='Do not take into account the BAF in the Bayesian inference, but only the depth ratio.')
   parser_model.add_argument('-f', "--segment-filter", dest = 'segfilt', type = float, default = 10e6,
                   help='Size in base-pair, to filter the segments to use in the Bayesian inference. Default 10e6.')
   parser_model.add_argument('-l', "--priors", dest = 'priors', type = str, default = '{"CN" :[1, 2, 3, 4], "value" : [1, 2, 1, 1]}',
                   help='Set the priors on the copy-number. Default 2 on CN = 2, 1 for all the other CN state \'{"CN" :[1, 2, 3, 4], "value" : [1, 2, 1, 1]}\'.')
   return parser.parse_args()

def sequenzaOverride(parser, subparser):
   parser_io      = subparser.add_argument_group(title='Input',description='Input options.')
   parser_param   = subparser.add_argument_group(title='Parameters',description='Model Parameters.')
   parser_gender  = subparser.add_argument_group(title='Gender',description='Option to control the gender and X/Y chromosome.')
   parser_io.add_argument('--dir', dest = 'dir', required = True,
                   help='The directory where the data to load are stored')
   parser_io.add_argument('-t', '--tag', dest = 'tag', default = None,
                   help='Tag indicating the prefix of the data, if not specified it is assumed as the name of the contanitor directory')
   parser_param .add_argument('-c', '--cellularity', required = True, type = float,
                   help = 'Cellularity parameter to pass to the model.')
   parser_param .add_argument('-p', '--ploidy', required = True, type = float,
                   help = 'Ploidy parameter to pass to the model.')
   parser_gender.add_argument('--is-male', dest = 'isMale', action='store_true', default = False,
                   help='Boolen flag indicating if the sequencing data are from a male or female, and consequently properly handle chromosome X and Y')
   parser_gender.add_argument('-X', "--chrX", dest = 'X', type = str, default = "X",
                   help='Character defining chromosome X. Default X.')
   parser_gender.add_argument('-Y', "--chrY", dest = 'Y', type = str, default = "Y",
                   help='Character defining chromosome Y. Default Y.')
   parser_param.add_argument('-r', "--only-ratio", dest = 'onlyratio', action='store_true', default = False,
                   help='Do not take into account the BAF in the Bayesian inference, but only the depth ratio.')
   return parser.parse_args()

def main():
   '''
   Execute the function with args
   '''
   parser = DefaultHelpParser(prog = __file__, formatter_class=lambda prog: SubcommandHelpFormatter(prog, max_help_position=20, width=75),
                              description='Sequenza Utils is an ensemble of tools capable of perform various tasks, primarily aimed to convert bam/pileup files to a format usable by the sequenza R package',
                              usage= '%(prog)s module [options]', epilog = 'This is version {0} - Francesco Favero - 2013'.format(VERSION))
   subparsers = parser.add_subparsers(dest='module')
   subparsers.metavar = None
   parser_pileup2abfreq  = subparsers.add_parser('pileup2abfreq', help = ' given a paired set of pileup (normal and matching tumor), and GC-content genome-wide information returns the common positions with A and B alleles frequencies',formatter_class=lambda prog: SubcommandHelpFormatter(prog,max_help_position=39, width=90))
   parser_pup2mu = subparsers.add_parser('pileup2acgt', help = 'convert pileup format to ACGT format',formatter_class=lambda prog: SubcommandHelpFormatter(prog,max_help_position=30, width=90))
   parser_gc_window  = subparsers.add_parser('GC-windows', help = 'Given a fasta file and a window size it computes the GC percentage across the sequences, and returns a file in the same format as gc5Base from UCSC')
   parser_merge_pileups = subparsers.add_parser('merge-pileups', help = 'Merging two pileups, it finds the common positions and return an mpileup file adding the second pilep as last 3 columns.')
   if RPY2 == True:
      parser_squeezeAB  = subparsers.add_parser('sequenzaExtract', help = 'Uses the R sequenza package to extract and save meaningful information from an ABfreq file')
      parser_doAllSequenza = subparsers.add_parser('sequenzaFit', help = 'Uses the R sequenza package to infer purity, ploidy and annotate mutation and CNV.')
      parser_overrideSequenza = subparsers.add_parser('sequenzaOverride', help = 'Uses the R sequenza package to calculate the CNV given the values for cellularity and ploidy.')
   try:
      used_module =  sys.argv[1]
      if used_module == "pileup2acgt":
         args = pileup2acgt(parser, parser_pup2mu)
         if args.chunk == 0:
            args.chunk = 1
            args.nproc  = 0
         if args.nproc >= 1:
            p = multiprocessing.Pool(processes=args.nproc)
         if not args.quiet:
            logging.basicConfig(format='%(message)s')
            start = time.clock()
            if args.pileup != '-':
               file_size = (os.stat(args.pileup).st_size/(1024*1024))
               logging.warning("Converting " + args.pileup + " -- size = %0.1f MB --" % file_size + " to ACGT..." )
            else:
               logging.warning("Converting " + args.pileup + " from STDIN to ACGT..." )
            logging.warning("Using chunks of " + str(args.chunk) + " line(s), and splitting the job in " + str(args.nproc+1)  + " process(es).")
         with xopen(args.output, "wb") as fileout:
            with xopen(args.pileup, "rb") as f:
               fileout.write('chr' + "\t" + 'n_base' + "\t" + 'ref_base' + "\t" +  'read.depth' + "\t" + 'A' + "\t" + 'C' + "\t" + 'G' + "\t" + 'T' + '\n')
               parse_pileup_partial = partial(parse_pileup_str, min_depth=args.n, qlimit=args.qlimit, qformat=args.qformat)
               counter = 0
               for chunk in grouper(args.chunk, f):
                  if args.nproc >= 1:
                     try:
                        results = p.map_async(parse_pileup_partial, chunk).get(99)
                     except AttributeError:
                        pass
                  else:
                     try:
                        results = map(parse_pileup_partial, chunk)
                     except AttributeError:
                        pass
                  for r in results:
                     counter = counter + 1
                     if r:
                        fileout.write(r + '\n')
            if not args.quiet:
               end = time.clock()
               seconds =  end-start
               logging.warning("Pileup to Mufreq: processed " + str(counter) + " lines in " + str(seconds) + " seconds")

      elif used_module == "pileup2abfreq":
         args = pileup2abfreq(parser, parser_pileup2abfreq)
         with xopen('-', "wb") as fileout:
            out_header = ["chromosome", "n.base", "base.ref", "depth.normal", "depth.sample", "depth.ratio", "Af", "Bf", "ref.zygosity", "GC.percent", "good.s.reads", "AB.germline", "AB.sample"]
            p1 = args.reference
            p2 = args.sample
            gc = args.gc
            stream_mpileup = IterableQueue()
            line_worker_partial = partial(line_worker, depth_sum=args.n, qlimit=args.qlimit, qformat=args.qformat, hom_t=args.hom, het_t=args.het)
            with xopen(p1, 'rb') as normal, xopen(p2, 'rb') as tumor, xopen(gc, 'rb') as gc_file:
               pup = multiPileups(normal,tumor)
               pup = GCmultiPileups(pup, gc_file)
               fileout.write("\t".join(out_header) + '\n')
               if args.chunk > 1 or args.nproc > 0:
                  #p = ThreadPool(processes=args.nproc)
                  p = multiprocessing.Pool(processes=args.nproc)
                  for res in p.imap(line_worker_partial, pup,chunksize=args.chunk):
                     #for res in results.get(99):
                     if res:
                        fileout.write('\t'.join(map(str,res))+'\n')
               else:
                  for line in pup:
                     res = line_worker_partial(line)
                     if res:
                        fileout.write('\t'.join(map(str,res))+'\n')


      elif used_module == "GC-windows":
         args = GC_windows(parser, parser_gc_window)
         gc_queue   = SimpleQueue()
         gc_process = multiprocessing.Process(target = process_gc_from_pipe, args = (gc_queue, args.window))
         gc_process.deamon = True
         gc_process.start()
         with xopen(args.fasta, 'rb') as fa_file:
            stream_fasta(fa_file, args.window, gc_queue)
         #gc_process.terminate()
         #gc_process.join()
      elif used_module == "merge-pileups":
         args = merge_pileups(parser, parser_merge_pileups)
         with xopen(args.p1, 'rb') as pileup1, xopen(args.p2, 'rb') as pileup2:
            pup = multiPileups(pileup1,pileup2)
            for line in pup:
               print('\t'.join(map(str, line)) + '\n')

      elif RPY2 == True and used_module == "sequenzaExtract":
         args = sequenzaExtract(parser, parser_squeezeAB)
         RPy2sqeezeABfreq(args.abfreq, args.loop, args.tag, check_dir(args.dir), args.kmin, args.gamma, args.mufreq)
      elif RPY2 == True and used_module == "sequenzaFit":
         args = sequenzaFit(parser, parser_doAllSequenza)
         priors_dict = json.loads(args.priors)
         RPy2doAllSequenza(check_dir(args.dir), args.isMale, args.tag, args.X, args.Y, args.ncpu, args.onlyratio, args.segfilt, priors_dict)
      elif RPY2 == True and used_module == "sequenzaOverride":
         args = sequenzaOverride(parser, parser_overrideSequenza)
         RPy2SequenzaOverride(check_dir(args.dir), args.isMale, args.tag, args.X, args.Y, args.onlyratio, args.cellularity, args.ploidy)
      else:
         return parser.parse_args()

   except IndexError:
      args = parser.parse_args()

if __name__ == "__main__":
   main()
