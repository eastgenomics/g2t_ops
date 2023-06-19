import argparse
from pathlib import Path

from g2t_ops import convert_symbols, utils, transcript_assigner


def main(args):
    output_path = utils.create_output_folder(args["output_folder"])
    hgnc_dump = utils.parse_tsv(args["hgnc_dump"])

    if args["command"] == "convert_symbols":
        if args["symbols"]:
            symbols = args["symbols"]
        else:
            symbols = utils.parse_file(args["symbol_file"])

        symbols = [symbol.upper() for symbol in symbols]

        hgnc_ids = convert_symbols.convert_symbols(symbols, hgnc_dump)
        convert_symbols.write_hgnc_ids(hgnc_ids, output_path)

    elif args["command"] == "assign_transcripts":
        db, user, pwd = args["hgmd_credentials"]
        session, meta = utils.connect_to_local_database(user, pwd, db)

        if args["hgnc_ids"]:
            for arg in args["hgnc_ids"]:
                assert arg.startswith("HGNC"), f"{arg} does not start with HGNC"

            hgnc_ids = args["hgnc_ids"]
        else:
            hgnc_ids = utils.parse_file(args["hgnc_file"])

        exon_file_data = utils.parse_exon_file(args["exon_file"])
        g2t_data = transcript_assigner.get_transcripts(
            hgnc_ids, exon_file_data
        )
        mane_data = transcript_assigner.parse_mane_file(
            args["mane_select"], hgnc_dump
        )
        clinical_tx_data = transcript_assigner.assign_transcripts(
            session, meta, mane_data, g2t_data
        )

        transcript_assigner.write_g2t(clinical_tx_data, output_path)
        transcript_assigner.write_sql_queries(clinical_tx_data, output_path)
        transcript_assigner.write_transcript_status(
            clinical_tx_data, output_path
        )
        print(f"Created output in {output_path}")
    else:
        print((
            f"{args['command']} is not an accepted command. Please use "
            "'convert_symbols' or 'assign_transcripts'"
        ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-o", "--output_folder", default=f"./{utils.get_date()}_results",
        help="Output folder in which to create the output files"
    )
    parser.add_argument(
        "-hgnc_dump", "--hgnc_dump", required=True,
        help="Path to HGNC dump downloaded from genenames.org"
    )

    subparser = parser.add_subparsers(dest="command")

    genes = subparser.add_parser(
        "convert_symbols", help="Convert list of gene symbols to HGNC ids"
    )
    genes.add_argument(
        "-symbols", "--symbols", nargs="+", help="Symbols to convert"
    )
    genes.add_argument(
        "-symbol_file", "--symbol_file", type=Path, help=(
            "File containing the symbols to convert"
        )
    )

    assigner = subparser.add_parser(
        "assign_transcripts", help=(
            "Assign clinical transcript to list of genes"
        )
    )
    assigner.add_argument(
        "-hgnc_ids", "--hgnc_ids", nargs="+", help="HGNC ids to convert"
    )
    assigner.add_argument(
        "-hgnc_file", "--hgnc_file", type=Path, help=(
            "File containing the HGNC ids to convert"
        )
    )
    assigner.add_argument(
        "-mane", "--mane_select", required=True,
        help="MANE Select file downloaded from http://tark.ensembl.org/web/mane_GRCh37_list/"
    )
    assigner.add_argument(
        "-hgmd", "--hgmd_credentials", required=True, nargs="+",
        help="Name of the MySQL HGMD database as well a read user and its password"
    )
    assigner.add_argument(
        "-exon_file", "--exon_file",
        default="project-Fkb6Gkj433GVVvj73J7x8KbV:file-GF611Z8433Gk7gZ47gypK7ZZ",
        help="Exon file id in DNAnexus"
    )

    args = vars(parser.parse_args())
    main(args)
