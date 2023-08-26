from melc import generate_protein_importance_for_clustering
import argparse


def _get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('adata_pickle_path', type=str, help='Please specify path to anndata pickle data!')
    parser.add_argument('essential_proteins_per_patient_pickle_path', type=str, help='Please specify path to essential proteins dictionary generated from the clustering step!')
    return parser


if __name__ == '__main__':
    args = _get_parser().parse_args()
    generate_protein_importance_for_clustering.generate_protein_importance_for_clustering(args.adata_pickle_path, args.essential_proteins_per_patient_pickle_path)
    