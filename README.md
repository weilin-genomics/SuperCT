# SuperCT
This project is for the development of a supervised classifier of cell types based on the single-cell RNAseq digital expression profiles.

You may run the standalone command as follow:

python pred_types.py -i <input_dge.csv> -o <cell_type_output.csv> -p -m <model.h5> -s <human or mouse> -t <id2type.csv>
