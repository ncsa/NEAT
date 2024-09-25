from .utils import make_qual_score_list, apply_markov_chain, save_file
from ..models import QualityScoreModel

__all__ = [
    "quality_score_model_runner"
]


def quality_score_model_runner(bam_file_path, csv_file_path, pickle_file_path):
    quality_df = make_qual_score_list(bam_file_path)
    markov_preds_df = apply_markov_chain(quality_df)
    save_file(markov_preds_df, csv_file_path, pickle_file_path)

    # ?
    # final_model = QualityScoreModel(markov_preds_df)
    # save_file(final_model, csv_file_path, pickle_file_path)
