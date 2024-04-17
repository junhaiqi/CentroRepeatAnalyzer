import os
import sys
import subprocess
import argparse
from tqdm import tqdm
import numpy as np
import pyfastx
from statistics import mode

cdhit_max_seq_name_len = 19  # The maximum length of sequence names cd-hit displays is 19.
class AnalyzerPipeline: 
    def __init__(self, fa_file, out_dir_path, rep_redio_cutoff = 0.5, rep_count_cutoff = 10, \
                 len_cutoff = 5000, k_mer_size = 11, hpc = False, clustering_idn_cutoff = 0.9,
                 rep_idn = 0.85, cluster_size_cutoff = 5, thread_num = 4, mono_tem_file = None):
        self.fa_file = fa_file
        self.rep_redio_cutoff = rep_redio_cutoff
        self.out_dir_path = out_dir_path
        self.len_cutoff = len_cutoff
        self.rep_count_cutoff = rep_count_cutoff
        self.k_mer_size = k_mer_size
        self.hpc = hpc
        self.clustering_idn_cutoff = clustering_idn_cutoff
        self.rep_idn = rep_idn
        self.cluster_size_cutoff = cluster_size_cutoff
        self.thread_num = thread_num
        self.mono_tem_file = mono_tem_file

    def scan_potential_repeat_reads(self):
        
        if os.path.exists(self.out_dir_path) == False:
            os.makedirs(self.out_dir_path)

        cmd = f'./bin/centroAnno {self.fa_file} \
                -L {self.len_cutoff} \
                -r {self.rep_redio_cutoff} \
                -S true -o {self.out_dir_path} \
                > {self.out_dir_path}/potential_repeat_reads.fa'

        print(f"Scanning potential repeat reads from {self.fa_file}...")
        os.system(cmd)
        print("Scan is over!")

    def anno_repeat_reads(self, potential_rep_fa_file):
        if os.path.exists(self.out_dir_path) == False:
            os.makedirs(self.out_dir_path)

        if self.hpc == True:
            cmd = f'./bin/centroAnno {potential_rep_fa_file} \
                    -L {self.len_cutoff} \
                    -r {self.rep_redio_cutoff} \
                    -k {self.k_mer_size} \
                    -c {self.hpc} \
                    -o {self.out_dir_path} \
                    > /dev/null'
        else:
            cmd = f'./bin/centroAnno {potential_rep_fa_file} \
                    -L {self.len_cutoff} \
                    -r {self.rep_redio_cutoff} \
                    -k {self.k_mer_size} \
                    -o {self.out_dir_path} \
                    > /dev/null'

        print(f"Annotating {potential_rep_fa_file}...")
        os.system(cmd)
        print("Annotation is over!")

    def scan_decompose_res_get_prefix_list(self, _dir):
        dir_list = os.listdir(path = _dir)
        prefix_set = set()
        for item in dir_list:
            idx = 0
            for i in range( len(item) ):
                if item[i] == '_':
                    idx = i
            prefix_set.add( item[0:idx] )
        return list( prefix_set )

    def fa_dict(self, fa_file, seq_name_len = cdhit_max_seq_name_len):
        fa = pyfastx.Fastx(fa_file)
        seq_dict = {}
        for name,seq in fa:
            seq_dict[ name[0:seq_name_len] ] = seq
        return seq_dict

    def get_centro_block_indentity_list(self, decompose_csv_file):
        inden_list = []
        with open(decompose_csv_file) as f:
            lines = f.readlines()
            for line in lines[1:]:
                idx1 = 0
                idx2 = 0
                t = 0
                for i in range( len(line) - 1, 0, -1):
                    if line[i] == ',':
                        if t == 0:
                            idx1 = i
                            t += 1
                        elif t == 1:
                            idx2 = i
                            break
                        else:
                            t += 1
                inden_list.append( float(line[idx2+1:idx1]) )

        return inden_list

    def get_monos_of_centroAnno_out(self, decompose_csv_file):
        mono_list = []
        with open(decompose_csv_file) as f:
            lines = f.readlines()
            for line in lines[1:]:
                idx1 = 0
                idx2 = 0
                t = 0
                for i in range( len(line) - 1, 0, -1):
                    if line[i] == ',':
                        if t == 3:
                            idx1 = i
                            t += 1
                        elif t == 4:
                            idx2 = i
                            break
                        else:
                            t += 1
                mono_list.append( line[idx2+1:idx1] )

        return mono_list

    def get_potential_repeats(self, _dir, repeat_fa):
        prefix_name_list = self.scan_decompose_res_get_prefix_list(_dir)
        r_f = open(repeat_fa, 'w')
        if self.mono_tem_file != None:
            cen_file = open(f'{self.out_dir_path}/monomers-related_reads.fa', 'w')
            monos = self.fa_dict( self.mono_tem_file )
            mono_lens = [ len(monos[key]) for key in monos]
            potential_repeat_reads = self.fa_dict(f'{self.out_dir_path}/potential_repeat_reads.fa', 1000)  # Don't include spaces in $sequence name (fastX file)! 
        print("Scanning for potential repeats and monomers-related reads...")
        for item in tqdm(prefix_name_list):
            _file = _dir + '/' + f'{item}_decomposedResult.csv'
            try:
                score = self.get_centro_block_indentity_list( _file )
                if np.mean( score ) > self.rep_idn:
                    mono_list = self.get_monos_of_centroAnno_out(_file)
                    mode_mono = mode(mono_list)
                    if mode_mono[-1] != "'" and mono_list.count(mode_mono) > self.rep_count_cutoff:
                        fa_file = _dir + '/' + f'{item}_monomerTemplates.fa'
                        seq_dict = self.fa_dict(fa_file)
                        r_f.write(f'>{item},repeat length = { len ( seq_dict[mode_mono] ) }\n')
                        r_f.write(f'{seq_dict[mode_mono]}\n')
                        if self.mono_tem_file != None:
                            is_cen = self.is_centromeric_read(mono_lens, len ( seq_dict[mode_mono] ))
                            if is_cen:
                                cen_file.write(f'>{item}\n')
                                cen_file.write(f'{ potential_repeat_reads[item] }\n')
            except:
                pass
        r_f.close()
        if self.mono_tem_file != None:
            cen_file.close()
        print("Scan is over.")

    def clustering_for_repeats(self, repeat_fa):
        print("Clustering to find repeats with high quality...")
        cmd = f'./bin/cd-hit -i {repeat_fa} \
                -c {self.clustering_idn_cutoff} \
                -o {self.out_dir_path}/rep_clusters \
                -M 0 \
                -T {self.thread_num} \
                 > /dev/null'
        os.system(cmd)
        print("Clustering ends.")

    def read_cdhit_res(self, cluster_file):  # get the clusters which size larger $self.rep_count_cutoff
        cluster_dict = {}
        with open(cluster_file) as f:
            lines = f.readlines()
            cluster_name = ''
            for line in lines:
                if line.startswith('>Cluster'):
                    cluster_name = line.strip('\n')
                    cluster_dict[ cluster_name ] = []
                else:
                    if line.strip('\n')[-1] != '*':
                        cluster_info = line.strip('\n').split('>')[-1].split('...')[0]
                    else:
                        cluster_info = line.strip('\n').split('>')[-1]
                    cluster_dict[ cluster_name ].append( cluster_info )

        if cluster_dict == {}:
            print('No repeats are inferred!')
            exit(-1)

        high_quality_clusters = {}
        for cluster in cluster_dict:
            if len( cluster_dict[ cluster ] ) > self.cluster_size_cutoff:
                high_quality_clusters[ cluster ] = cluster_dict[ cluster ]
        
        return high_quality_clusters

    def get_seq_len_mode(self, seq_dict, name_list):
        len_list = []
        for key in name_list:
            if key[-1] != '*':
                len_list.append( len(seq_dict[ key ]) )
            else:
                len_list.append(-1)
        mode_len = mode(len_list)
        mode_idx = 0
        for i in range( len(len_list) ):
            if len_list[i] == mode_len:
                mode_idx = i
                break
        return name_list[mode_idx], seq_dict[ name_list[mode_idx][0:cdhit_max_seq_name_len] ]

    def get_representative_repeats(self, cluster_file, repeat_fa):
        clustering_res = self.read_cdhit_res(cluster_file)
        if clustering_res == {}:
            print('No repeats are inferred!')
            return 0
        seq_dict = self.fa_dict(repeat_fa)
        representative_repeats_file = f'{self.out_dir_path}/representative_repeats.fa'
        _file = open(representative_repeats_file, 'w')
        for cluster in clustering_res:
            try:
                info = self.get_seq_len_mode(seq_dict, clustering_res[cluster])
                rep_name = info[0]
                rep_seq = info[1]
                _file.write(f'>{rep_name},length = {len(rep_seq)},read count = {len(clustering_res[cluster])}\n')
                _file.write(f'{rep_seq}\n')
            except:
                pass
        _file.close()

    def get_centromeric_loss(self, true_len = 171, now_len = 345):
        if now_len / true_len > 0.9:
            p1 = (now_len / true_len) - int(now_len / true_len)
            loss = min( abs(1-p1), p1 )
        else:
            loss = 1
        return loss

    def is_centromeric_read(self, mono_lens, this_len, loss_cutoff = 0.05):
        is_cen = False
        for _len in mono_lens:
            loss = self.get_centromeric_loss(_len, this_len)
            if loss < loss_cutoff:
                is_cen = True
                break
        return is_cen

    def refine_anno_of_mono_related_reads(self):
        if self.mono_tem_file != None:
            print('The annotation results of monomers-related these reads are being refined...')
            if self.hpc == True:
                cmd = f'./bin/centroAnno {self.out_dir_path}/monomers-related_reads.fa \
                        -L {self.len_cutoff} \
                        -r {self.rep_redio_cutoff} \
                        -k {self.k_mer_size} \
                        -c {self.hpc} \
                        -o {self.out_dir_path}/monomers-related_refine_annotations \
                        -m {self.mono_tem_file} \
                        > /dev/null'
            else:
                cmd = f'./bin/centroAnno {self.out_dir_path}/monomers-related_reads.fa \
                        -L {self.len_cutoff} \
                        -r {self.rep_redio_cutoff} \
                        -k {self.k_mer_size} \
                        -o {self.out_dir_path}/monomers-related_refine_annotations \
                        -m {self.mono_tem_file} \
                        > /dev/null'
            os.system(cmd)
            print('Refine ends.')

def get_parameters():
    parser = argparse.ArgumentParser()

    parser.add_argument('fastX_file', type = str,
                        help = 'It is the input filename in fasta format, required, can be in .gz format.')

    parser.add_argument('output_folder_path', type = str,
                        help = 'It is the path to the output folder, required.')

    parser.add_argument('-k', '--kmer-size', type = int, required = False, default = 9,
                        help = 'It is the size of k-mer, which is used to filter repeat sequences and infer repeats (default = 9).')

    parser.add_argument('-l', '--len-cutoff', type = int, required = False, default = 5000,
                        help = 'It is the length cutoff, reads longer than this value are considered repeat sequences (default = 5000).')

    parser.add_argument('-c', '--count-cutoff', type = int, required = False, default = 10,
                        help = 'It is the count cutoff, repeats that appear more than this value are considered high-quality repeats (default = 10).')

    parser.add_argument('-r', '--rep-redio-cutoff', type = float, required = False, default = 0.5,
                        help = 'It is the redio cutoff of non-unique k-mer, sequences greater than this value are considered potential repeat sequences (default = 0.5).')

    parser.add_argument('-p', '--close-hpc', action = "store_true", 
                        help = 'It is the bool flag, when it appears, homopolymer compression technology will not be applied (default = False).')

    parser.add_argument('-s', '--cluster-size-cutoff', type = int, required = False, default = 5,
                        help = 'It is the size cutoff, the clusters which size larger this values are used to inferred repeats (default = 5).')

    parser.add_argument('-d', '--clus-identity', type = float, required = False, default = 0.9,
                        help = 'It is the sequence identity threshold for cdhit (default = 0.9).')

    parser.add_argument('-e', '--rep-identity', type = float, required = False, default = 0.85,
                        help = 'It is the identity threshold to filter reads with high-quality annotation results (default = 0.85).')

    parser.add_argument('-m', '--monomer-template-file', type = str, required = False, default = None,
                        help = 'It is the fasta file which includes monomer templates (default = None), If it is given, the program will scan out all reads potentially related to these monomers and complete more accurate annotations based on these monomers.')

    parser.add_argument('-y', '--only-infer', action = "store_true",
                        help = 'It is the bool flag, when it appears, it means that there is already annotation results, and the annotation results are stored in $output_folder_path (default = False).')

    parser.add_argument('-t', '--thread-num', type=int, required=False, default = 4,
                        help = 'It is the number of threads used to clustering (default = 4).')

    args = parser.parse_args()

    return args

def main():

    args = get_parameters()

    pipeline = AnalyzerPipeline(fa_file = args.fastX_file, out_dir_path = args.output_folder_path, \
                                rep_redio_cutoff = args.rep_redio_cutoff, rep_count_cutoff = args.count_cutoff, \
                                len_cutoff = args.len_cutoff, k_mer_size = args.kmer_size, hpc = args.close_hpc, 
                                clustering_idn_cutoff = args.clus_identity, cluster_size_cutoff = args.cluster_size_cutoff, \
                                rep_idn = args.rep_identity, thread_num = args.thread_num, mono_tem_file = args.monomer_template_file)

    if args.only_infer == False:
        pipeline.scan_potential_repeat_reads()

        pipeline.anno_repeat_reads(potential_rep_fa_file = f'{args.output_folder_path}/potential_repeat_reads.fa')
        rep_fa_path = f'{args.output_folder_path}/potential_repeats.fa'

        pipeline.get_potential_repeats(_dir = args.output_folder_path, repeat_fa = rep_fa_path)

        pipeline.clustering_for_repeats(repeat_fa = rep_fa_path)
        cluster_file_path = f'{args.output_folder_path}/rep_clusters.clstr'

        pipeline.get_representative_repeats(cluster_file = cluster_file_path, repeat_fa = rep_fa_path)

        pipeline.refine_anno_of_mono_related_reads()
    else:
        rep_fa_path = f'{args.output_folder_path}/potential_repeats.fa'
        pipeline.get_potential_repeats(_dir = args.output_folder_path, repeat_fa = rep_fa_path)
        
        pipeline.clustering_for_repeats(repeat_fa = rep_fa_path)
        cluster_file_path = f'{args.output_folder_path}/rep_clusters.clstr'
        pipeline.get_representative_repeats(cluster_file = cluster_file_path, repeat_fa = rep_fa_path)

        pipeline.refine_anno_of_mono_related_reads()


if __name__ == "__main__":
    main()
