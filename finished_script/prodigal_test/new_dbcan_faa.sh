python ../main.py input_process --input_raw EscheriaColiK12MG1655.faa --mode protein --output_dir new_output_faa --input_format NCBI

python ../main.py CAZyme_annotation --output_dir new_output_faa --db_dir ../../part1_script/dbCAN_db/

python ../main.py CGC_info --output_dir new_output_faa --db_dir ../../part1_script/dbCAN_db/ --input_gff EscheriaColiK12MG1655.gff --input_gff_format NCBI_prok

python ../main.py CGC_annotation --output_dir new_output_faa

python ../main.py CGC_substrate_prediction -i new_output_faa --db_dir ../../part1_script/dbCAN_db/ --output_dir new_output_faa

python ../main.py CGC_substrate_plot -i new_output_faa/substrate.out -b new_output_faa/PUL_blast.out --cgc new_output_faa/cgc_standard_out.tsv --db_dir ../../part1_script/dbCAN_db --output_dir new_output_faa
