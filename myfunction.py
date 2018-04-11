#!/usr/bin/python
# -*- coding: utf-8 -*-
import commands
import string
import types
import sys, os
import platform
# import pysam

class MyError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def trans_code_by_os(s):
    if platform.system() == 'Windows':
        s = s.decode('gbk').encode('utf8')
    return s

def read_in_file_list(file_path):
    f = open(file_path,'rU')
    content = f.read()
    flist = content.split('\n')
    # 删去空白行
    while '' in flist:
        flist.remove('')
    return flist

def read_in_file_list_by_fid(f):
    content = f.read()
    flist = content.split('\n')
    # 删去空白行
    while '' in flist:
        flist.remove('')
    return flist

def call_cmd(cmd_str, verbose = True):
    if verbose:
        print cmd_str
    (status, output) = commands.getstatusoutput(cmd_str)
    if status:
        raise MyError(output)
    return output

def call_cmd_list(cmd_list):
    if type(cmd_list) is types.ListType or type(cmd_list) is types.TupleType:
        for cmd_str in cmd_list:
            call_cmd(cmd_str)
        return 0
    else:
        raise MyError("Unexpected data type of cmd_list! Should be ListType or TupleType.")
        

def call_cmd_with_retry(cmd_str):
    try:
        return call_cmd(cmd_str)
    except:
        return call_cmd(cmd_str)

def call_cmd_list_with_retry(cmd_list):
    try:
        return call_cmd_list(cmd_list)
    except:
        return call_cmd_list(cmd_list)

def count_pattern_per_read(fin, fout, option, expstr):
    total_num = 0
    line_id = 0

    def process_read(read_lines):
        if not len(read_lines) == 2 or not read_lines[0].startswith(">"):
            raise MyError('Format error occurs in line %d - %d.\n%s' %(line_id + 1, line_id + len(read_lines),lines))

        title =  read_lines[0]
        seq = read_lines[1]
        if option == '-s':
            num = string.count(seq, expstr)
        elif option == '-e':
            num = len(re.findall(expstr, seq))
        else:
            print "Undefined input of option: %s" %option
            sys.exit(1)
        title = title + '\t' + option + ':' + expstr + "=" + str(num)
        print >> fout, title + '\n' + seq
        return num

    # 读下一个以'>'开头的行之前的行，并用process_read函数进行处理
    curr_read = []
    for line in fin:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if curr_read:
                total_num += process_read(curr_read)
                line_id += len(curr_read)
                if line_id % 1000000 == 0:
                    print line_id
                curr_read = []
            else:
                line_id = 0
        curr_read.append(line)
    total_num += process_read(curr_read)
    print >> fout, 'total_num = %d' %total_num
    return total_num

def check_filepath(filepath):
    if not os.path.isfile(filepath):
        filepath = trans_code_by_os(filepath)
        raise MyError('Error! The file: %s does not exist.' %filepath)


def get_chrom_len(genome_length_file):
    fin = open(genome_length_file, 'rU')
    chrlen_map = {}
    for line in fin:
        line = line.strip()
        eles = line.split()
        if len(eles) != 2:
            raise MyError('Format error for genome length file: %s\n%s\n' %(genome_length_file, line))
        chrm = eles[0]
        chrmlen = int(eles[1])
        if (chrm in chrlen_map) and chrlen_map[chrm] != chrmlen:
            raise MyError('Duplicated chromosome with different length in genome length file: %s\nchrm:%s' 
                %(genome_length_file, chrm))
        chrlen_map[chrm] = chrmlen
    fin.close()
    return chrlen_map

def read_and_process_fasta_seq(genome_fasta_file, chroms = []):
    def get_length_and_seq(curr_seq):
        if len(curr_seq) <2 or not curr_seq[0].startswith('>'):
            raise MyError('Format error occurs around line %d\n%s.' %(line_id + 1, repr(curr_seq)))
        title_eles = curr_seq[0].lstrip('>').split()
        title = title_eles[0]
        seq = ''.join(curr_seq[1:])
        return (title, seq)

    def fill_in_dict():
        if title in seq_len_hash:
            raise MyError('Duplicated title in fasta file : %s' %genome_fasta_file)
        seq_len_hash[title] = seq.upper()

    fin = open(genome_fasta_file, 'rU')
    seq_len_hash = {}
    line_id = 0
    curr_seq = []
    for line in fin:
        line = line.strip()
        if line.startswith('>'):
            if curr_seq:
                line_id += len(curr_seq)
                (title, seq) = get_length_and_seq(curr_seq)
                if not chroms or title in chroms:
                    fill_in_dict()
                curr_seq = []
        curr_seq.append(line)
    (title, seq) = get_length_and_seq(curr_seq)
    fill_in_dict()
    fin.close()
    return seq_len_hash


# def open_sam_or_bam_file_for_read(filename, check_header=True, check_sq=True):
#     extension = filename.split('.')
#     if len(extension) < 2 or extension[-1] not in ['bam', 'sam']:
#         raise MyError('Unknown file type: %s' %filename)
#     elif extension[-1] == 'sam': 
#         mode = 'r'
#     else:
#         mode = 'rb'
#     if check_header == False: check_sq = False
#     return pysam.AlignmentFile(filename, mode, check_header=check_header, check_sq=check_sq)


# def open_sam_or_bam_file_for_write(filename, template=None, 
#     reference_names=None, reference_lengths=None, text=None, header=None, add_sq_text=False):
#     extension = filename.split('.')
#     if len(extension) < 2 or extension[-1] not in ['bam', 'sam']:
#         raise MyError('Unknown file type: %s' %filename)
#     elif extension[-1] == 'sam':
#         mode = 'wh' 
#     else:
#         mode = 'wb'
#     return pysam.AlignmentFile(filename, mode, template=template, reference_names=reference_names, 
#         reference_lengths=reference_lengths, text=text, header=header, add_sq_text=add_sq_text)


def complementary_and_reverse_seq(seq):
    trans_table = string.maketrans('aAtTcCgG','tTaAgGcC')
    compl_seq = string.translate(seq, trans_table)
    compl_reverse_seq = compl_seq[::-1]
    return compl_reverse_seq

