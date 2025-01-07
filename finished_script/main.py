import argparse
import os
import psutil
import pandas as pd
from dbCAN_database import DBDownloader
from input_process import InputProcessor
from diamond import DiamondProcessor
from pyhmmer_search import PyHMMERProcessor
from process_dbcan_sub import DBCANProcessor
from OverviewGenerator import OverviewGenerator
from generate_cgc_gff import GFFProcessor
from CGCFinder_pandas_class import CGCFinder


def run_dbCAN_database(dbCAN_database_config):
    downloader = DBDownloader(dbCAN_database_config)
    downloader.download_file()
    downloader.extract_tar_file()

def run_dbCAN_input_process(dbCAN_input_process_config):
    input_processor = InputProcessor(dbCAN_input_process_config)
    input_processor.process_input()

def run_dbCAN_CAZyme_annotation(method, diamond_config, hmm_config, dbcan_sub_config, generate_overview_config):
    if 'diamond' in method:
        diamond_processor = DiamondProcessor(diamond_config)
        diamond_processor.run_diamond()
        diamond_processor.format_results()
    if 'hmm' in method:
        hmm_processor = PyHMMERProcessor(hmm_config)
        hmm_processor.hmmsearch()
    if 'dbCANsub' in method:
        dbcan_sub_processor = PyHMMERProcessor(dbcan_sub_config)
        dbcan_sub_processor.hmmsearch()
        parser_dbcan_sub = DBCANProcessor(dbcan_sub_config)
        parser_dbcan_sub.process_dbcan_sub()
    overview_generator = OverviewGenerator(generate_overview_config)
    overview_generator.run()



def run_CGC_annotation_preprocess(CGC_info_config,tc_config,tf_config,stp_config):
    cgc_infoProcessor = GFFProcessor(CGC_info_config)
    cgc_infoProcessor.generate_non_cazyme_faa()
#    tc_hmm_processor   = PyHMMERProcessor(tc_config)
#    tc_hmm_processor.hmmsearch()   #used to process tcDoms but now use diamond

    tc_diamond_processor = DiamondProcessor(tc_config)
    tf_hmm_processor   = PyHMMERProcessor(tf_config)
    stp_hmm_processor  = PyHMMERProcessor(stp_config)
    tc_diamond_processor.run_tcdb_diamond()
    tc_diamond_processor.format_results_tcdb()
    tf_hmm_processor.hmmsearch()
    stp_hmm_processor.hmmsearch()

#    tc_df = pd.read_csv(tc_config['output_file'], sep='\t')

    columns = ['Annotate Name', 'Annotate Length', 'Target Name', 'Target Length', 'i-Evalue', 'Annotate From', 'Annotate To', 'Target From', 'Target To', 'Coverage', 'Annotate File Name']
    tc_df = pd.read_csv(tc_config['output_file'], names=columns, header=0, sep='\t')
    tf_df = pd.read_csv(tf_config['output_file'], names=columns, header=0, sep='\t')
    stp_df = pd.read_csv(stp_config['output_file'], names=columns, header=0, sep='\t')
    total_function_annotation_df = pd.concat([tc_df, tf_df, stp_df], ignore_index=True)
    tf_hmm_processor.filter_overlaps(total_function_annotation_df).to_csv(os.path.join(tf_config["output_dir"], 'total_cgc_info.tsv'), index=False, sep='\t')
    cgc_infoProcessor.process_gff()

def run_CGC_annotation(cgc_config):
    cgc_finder = CGCFinder(cgc_config)
    cgc_finder.read_gff()
    cgc_finder.mark_signature_genes()
    clusters = cgc_finder.find_cgc_clusters()   
    cgc_finder.output_clusters(clusters)


def main():
    parser = argparse.ArgumentParser(description='CAZyme analysis workflow management.')
    parser.add_argument('command', choices=['database', 'input_process', 'CAZyme_annotation', 'CGC_info', 'CGC_annotation'],
                        help='The function to run: \n'
'database - Downloads and prepares the dbCAN database files.\n'
'input_process - Processes the input fasta files for annotation.\n'
'CAZyme_annotation - Runs the CAZyme annotation using specified methods.\n'
'CGC_info - Prepares the input files for CGC annotation.\n'
'CGC_annotation - Runs the CGC annotation using specified methods.')
    parser.add_argument('--db_dir', default='./dbCAN_databases', help='Specify the directory to store the dbCAN database files.')
    parser.add_argument('--input_raw_data', help='Specify the input fasta file for data preprocessing.')
    parser.add_argument('--input_gff', help='Specify the input gff file for CGC annotation.')
    parser.add_argument('--mode', help='Check the mode of the input file.')
    parser.add_argument('--output_dir', default='./output', help='Output folder.')
    parser.add_argument('--input_format', default='NCBI', choices=['NCBI', 'JGI'], help='Specify the input format for protein sequences.')
    parser.add_argument('--input_gff_format', default='NCBI', choices=['NCBI_euk', 'JGI', 'NCBI_prok', 'prodigal'], help='Specify the input format for protein sequences.')

    parser.add_argument('--methods', nargs='+', choices=['diamond', 'hmm', 'dbCANsub'],
                        default=['diamond', 'hmm', 'dbCANsub'],help='Specify the annotation methods to use. Options are diamond, hmm, and dbCANsub.')

    parser.add_argument('--diamond_evalue', type=float, default=1e-102, help='E-value threshold for diamond annotation.')
    parser.add_argument('--diamond_threads', type=int, default=psutil.cpu_count(), help='Number of threads for diamond annotation.')
    parser.add_argument('--diamond_verbose_option', action='store_true', default=False,  help='Enable verbose output for diamond.')

    parser.add_argument('--dbcan_hmm_evalue', type=float, default=1e-15, help='E-value threshold for HMM annotation.')
    parser.add_argument('--dbcan_hmm_coverage', type=float, default=0.35, help='Coverage threshold for HMM annotation.')
    parser.add_argument('--hmmer_cpu',type=int,default=psutil.cpu_count(), help='Number of threads for HMM annotation.')

    parser.add_argument('--dbcansub_hmm_evalue', type=float, default=1e-15, help='E-value threshold for dbCANsub annotation.')
    parser.add_argument('--dbcansub_hmm_coverage', type=float, default=0.35, help='Coverage threshold for dbCANsub annotation.')


    parser.add_argument('--tc_evalue', type=float, default=1e-15, help='E-value threshold for TC annotation.')
    parser.add_argument('--tc_coverage', type=float, default=0.35, help='Coverage threshold for TC annotation.')
    parser.add_argument('--tf_evalue', type=float, default=1e-15, help='E-value threshold for TF annotation.')
    parser.add_argument('--tf_coverage', type=float, default=0.35, help='Coverage threshold for TF annotation.')
    parser.add_argument('--stp_evalue', type=float, default=1e-15, help='E-value threshold for STP annotation.')
    parser.add_argument('--stp_coverage', type=float, default=0.35, help='Coverage threshold for STP annotation.')


    parser.add_argument('--additional_genes', nargs='+', default=["TC"], help='Specify additional gene types for CGC annotation, including TC, TF, and STP')
    parser.add_argument('--num_null_gene', type=int, default=2, help='Maximum number of null genes allowed between signature genes.')
    parser.add_argument('--base_pair_distance', type=int, default=15000, help='Base pair distance of sig genes for CGC annotation.')
    parser.add_argument('--use_null_genes', action='store_true', default=True, help='Use null genes in CGC annotation.')
    parser.add_argument('--use_distance', action='store_true', default=False, help='Use base pair distance in CGC annotation.')


    args = parser.parse_args()

    if args.command == 'CAZyme_annotation' and not args.methods:
        print("Error: At least one method must be specified for CAZyme annotation.")
        parser.print_help()
        return

    if args.command == 'database' and not args.db_dir:
        print("Error: Database directory is required for database preparation.")
    elif args.command == 'database' and args.db_dir:
        dbCAN_database_config = {'db_dir': args.db_dir}
        run_dbCAN_database(dbCAN_database_config)

    
    if args.command == 'input_process':
        if args.input_raw_data and args.mode:
            dbCAN_input_process_config = {
                'input_raw_data': args.input_raw_data, 
                'mode': args.mode, 
                'output_dir': args.output_dir,
                'input_format': args.input_format}
            run_dbCAN_input_process(dbCAN_input_process_config)
        else:
            print("Error: Input file is required for input processing.")
            parser.print_help()
        
    if args.command == 'CAZyme_annotation':
        diamond_config = {
            'diamond_db': os.path.join(args.db_dir, 'CAZy.dmnd'),
            'input_faa': os.path.join(args.output_dir, 'uniInput.faa'),
            'output_file': os.path.join(args.output_dir, 'diamond_results.tsv'),
            'e_value_threshold': args.diamond_evalue,
            'threads': args.diamond_threads,
            'verbose_option': args.diamond_verbose_option
        }
        hmm_config = {
            'hmm_file': os.path.join(args.db_dir, 'dbCAN.hmm'),
            'input_faa': os.path.join(args.output_dir, 'uniInput.faa'),
            'output_file': os.path.join(args.output_dir, 'dbCAN_hmm_results.tsv'),
            'e_value_threshold': args.dbcan_hmm_evalue,
            'coverage_threshold': args.dbcan_hmm_coverage,
            'hmmer_cpu': args.hmmer_cpu
        }
        dbcan_sub_config = {
            'hmm_file': os.path.join(args.db_dir, 'dbCAN_sub.hmm'),
            'input_faa': os.path.join(args.output_dir, 'uniInput.faa'),
            'e_value_threshold': args.dbcansub_hmm_evalue,
            'coverage_threshold': args.dbcansub_hmm_coverage,
            'hmmer_cpu': args.hmmer_cpu,
            'output_file': os.path.join(args.output_dir, 'dbCAN-sub.substrate.tsv'),
            'mapping_file': os.path.join(args.db_dir, 'fam-substrate-mapping.tsv')
        }
        generate_overview_config = {'output_dir': args.output_dir}

        run_dbCAN_CAZyme_annotation(args.methods, diamond_config, hmm_config, dbcan_sub_config, generate_overview_config)

    if args.command == 'CGC_info':

        CGC_info_config = { 
            'input_total_faa': os.path.join(args.output_dir, 'uniInput.faa'),
            'output_dir': args.output_dir,
            'cazyme_overview': os.path.join(args.output_dir, 'overview.tsv'),
            'cgc_sig_file': os.path.join(args.output_dir, 'total_cgc_info.tsv'),
            'input_gff': args.input_gff,
            'output_gff': os.path.join(args.output_dir, 'cgc.gff'),
            'gff_type': args.input_gff_format
        }

        tc_config = {
            'diamond_db': os.path.join(args.db_dir, 'tcdb.dmnd'),
            'input_faa': os.path.join(args.output_dir, 'non_CAZyme.faa'),
            'output_file': os.path.join(args.output_dir, 'TC_results.tsv'),
            'e_value_threshold': args.tc_evalue,
            'coverage_threshold_tc': args.tc_coverage,
            'threads': args.diamond_threads,
            'verbose_option': args.diamond_verbose_option
        }

        tf_config = {
            'hmm_file': os.path.join(args.db_dir, 'TF.hmm'),
            'input_faa': os.path.join(args.output_dir, 'non_CAZyme.faa'),
            'output_file': os.path.join(args.output_dir, 'TF_results.tsv'),
            'e_value_threshold': args.tf_evalue,
            'coverage_threshold': args.tf_coverage,
            'hmmer_cpu': args.hmmer_cpu,
            'output_dir': args.output_dir
        }

        stp_config = {
            'hmm_file': os.path.join(args.db_dir, 'STP.hmm'),
            'input_faa': os.path.join(args.output_dir, 'non_CAZyme.faa'),
            'output_file': os.path.join(args.output_dir, 'STP_results.tsv'),
            'e_value_threshold': args.stp_evalue,
            'coverage_threshold': args.stp_coverage,
            'hmmer_cpu': args.hmmer_cpu,
            'output_dir': args.output_dir
        }

        run_CGC_annotation_preprocess(CGC_info_config,tc_config,tf_config,stp_config)


    if args.command == 'CGC_annotation':
        cgc_config = {
            'output_dir': args.output_dir,
            'filename': os.path.join(args.output_dir, 'cgc.gff'),
            'num_null_gene': args.num_null_gene,
            'base_pair_distance': args.base_pair_distance,
            'use_null_genes': args.use_null_genes,
            'use_distance': args.use_distance,
            'additional_genes': args.additional_genes
        }
        run_CGC_annotation(cgc_config)



if __name__ == '__main__':
    main()