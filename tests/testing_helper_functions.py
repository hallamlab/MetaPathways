""" Functions that run MetaPathways specific operations and return outputs that can be tested """
import os
from .testing_utils import get_test_data


def run_func_search(tmp_dir: str, db_name: str, test_sample_name: str) -> (str, str):
    from metapathways import MetaPathways_func_search

    faa_input = get_test_data(
        os.path.join(
            "output", test_sample_name, "orf_prediction", test_sample_name + ".qced.faa"
        )
    )

    ref_blast_db = get_test_data(
        os.path.join("ref_data", "functional", "formatted", db_name)
    )

    os.makedirs(
        os.path.join(tmp_dir, test_sample_name, "orf_prediction"), exist_ok=True
    )
    os.makedirs(os.path.join(tmp_dir, test_sample_name, "blast_results"), exist_ok=True)

    output_blast_results = os.path.join(
        tmp_dir,
        test_sample_name,
        "blast_results",
        test_sample_name + "." + db_name + ".BLASTout",
    )

    expect_output_blast_results = get_test_data(
        os.path.join(
            "output",
            test_sample_name,
            "blast_results",
            test_sample_name + "." + db_name + ".BLASTout",
        )
    )

    args = [
        "--algorithm",
        "BLAST",
        "--blast_executable",
        "blastp",
        "--num_threads",
        "4",
        "--blast_max_target_seqs",
        "5",
        "--blast_outfmt",
        "6",
        "--blast_query",
        faa_input,
        "--blast_evalue",
        "0.000001",
        "--blast_db",
        ref_blast_db,
        "--blast_out",
        output_blast_results,
    ]

    MetaPathways_func_search.main(args)

    return output_blast_results, expect_output_blast_results


def run_parse_blast(tmp_dir: str, db_name: str, sample_name: str) -> (str, str):
    from metapathways import MetaPathways_parse_blast

    ref_blast_db_annots = get_test_data(
        os.path.join("ref_data", "functional", "formatted", db_name + "-names.txt")
    )

    input_blast_results = get_test_data(
        os.path.join(
            "output",
            sample_name,
            "blast_results",
            sample_name + "." + db_name + ".BLASTout",
        )
    )

    input_blast_refscores = get_test_data(
        os.path.join(
            "output", sample_name, "blast_results", sample_name + ".refscores.BLAST"
        )
    )

    os.makedirs(os.path.join(tmp_dir, sample_name, "blast_results"), exist_ok=True)

    output_parsed_blast_results = os.path.join(
        tmp_dir,
        sample_name,
        "blast_results",
        sample_name + "." + db_name + ".BLASTout.parsed.txt",
    )

    expected_parsed_blast_results = get_test_data(
        os.path.join(
            "output",
            sample_name,
            "blast_results",
            sample_name + "." + db_name + ".BLASTout.parsed.txt",
        )
    )
    args = [
        "-b",
        input_blast_results,
        "-r",
        input_blast_refscores,
        "-m",
        ref_blast_db_annots,
        "-d",
        db_name,
        "-o",
        output_parsed_blast_results,
        "--min_bsr",
        "0.4",
        "--min_score",
        "20",
        "--min_length",
        "45",
        "--max_evalue",
        "0.000001",
        "--algorithm",
        "BLAST",
    ]

    MetaPathways_parse_blast.main(args)
    return output_parsed_blast_results, expected_parsed_blast_results


def run_rrna_stats_calculator(
    out_folder_name: str, test_sample: str, rna_db_id: str
) -> (str, str):
    from metapathways import MetaPathways_rRNA_stats_calculator

    ref_blast_db = get_test_data(os.path.join("ref_data", "taxonomic", rna_db_id))

    os.makedirs(
        os.path.join(out_folder_name, test_sample, "blast_results"), exist_ok=True
    )
    os.makedirs(
        os.path.join(out_folder_name, test_sample, "results", "rRNA"), exist_ok=True
    )

    input_blast_results = get_test_data(
        os.path.join(
            "output",
            test_sample,
            "blast_results",
            test_sample + ".rRNA." + rna_db_id + ".BLASTout",
        )
    )

    output_rrna_stats = get_test_data(
        os.path.join(
            "output",
            test_sample,
            "results",
            "rRNA",
            test_sample + "." + rna_db_id + ".rRNA.stats.txt",
        )
    )
    expected_rrna_stats = get_test_data(
        os.path.join(
            "output",
            test_sample,
            "results",
            "rRNA",
            test_sample + "." + rna_db_id + ".rRNA.stats.txt",
        )
    )
    args = [
        "-o",
        output_rrna_stats,
        "-b",
        "50",
        "-e",
        "0.000001",
        "-s",
        "20",
        "-i",
        input_blast_results,
        "-d",
        ref_blast_db,
    ]

    MetaPathways_rRNA_stats_calculator.main(args)
    return output_rrna_stats, expected_rrna_stats
