# SuperCT
This project is for the development of a supervised classifier of cell types based on the single-cell RNAseq digital expression profiles.

You may run the standalone command as follow:

python pred_types.py -i <input_dge.csv> -o <cell_type_output.csv> -p -m <model.h5> -s <human/mouse> -t <id2type.csv>

To prepare input_dge.csv from 10xGenomics output files for the above prediction command, you may run the command as follow:

 python transform_10x_to_dge.py -m <matrix.mtx> -g <genes.tsv> -b <barcodes.tsv> -o <transformed_dge.csv>
