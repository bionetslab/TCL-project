from melc import DA_randomly_assign_cellTypeAnnotations
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    DA_randomly_assign_cellTypeAnnotations.da_randomly_assign_cell_type_annotations(args.adata_pickle_path)