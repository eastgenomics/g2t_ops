import argparse
from collections import defaultdict
import datetime
import gzip
import logging
import os
from pathlib import Path
import regex

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.schema import MetaData
import xlrd

from hgnc_queries import get_id as hq_get_id


def create_log(output_folder):
    logging.basicConfig(
        filename=f"{output_folder}/g2t.log", level=logging.DEBUG
    )
    logger = logging.getLogger(__name__)

    return logger


def get_date():
    """ Return today's date in YYMMDD format

    Returns:
        str: Date
    """

    return str(datetime.date.today())[2:].replace("-", "")


def connect_to_db(user: str, passwd: str, host: str, database: str):
    """ Return cursor of panel_database

    Args:
        user (str): Username for the database
        passwd (str): Password for the user
        host (str): Host for the database
        database (str): Name of the database to connect to

    Returns:
        tuple: SQLAlchemy session obj, SQLAlchemy meta obj
    """

    try:
        db = create_engine(
            f"mysql://{user}:{passwd}@{host}/{database}"
        )
    except Exception as e:
        raise e
    else:
        meta = MetaData()
        meta.reflect(bind=db)
        Session = sessionmaker(bind=db)
        session = Session()
        return session, meta


def parse_hgnc_dump(hgnc_file: str):
    """ Parse the hgnc dump and return a dict of the data in the dump

    Args:
        hgnc_file (str): Path to the hgnc file

    Returns:
        dict: Dict of hgnc data, symbol data, alias data, previous symbol data
    """

    data = {}
    symbol_dict = {}
    alias_dict = {}
    prev_dict = {}

    with open(hgnc_file) as f:
        for i, line in enumerate(f):
            # first line is headers
            if i == 0:
                reformatted_headers = []
                headers = line.strip().split("\t")

                for header in headers:
                    # need transform the header name to the table attribute name
                    if "supplied" in header:
                        # external links provided always have: "(supplied by ...)"
                        # split on the (, get the first element, strip it
                        # (there's spaces sometimes), lower the characters and
                        # replace spaces by underscores
                        header = header.split("(")[0].strip().lower().replace(" ", "_")
                        # they also need ext_ because they're external links
                        header = f"ext_{header}"
                    else:
                        header = header.lower().replace(" ", "_")

                    reformatted_headers.append(header)

                # gather positions for the following headers
                symbol_index = get_header_index(
                    "approved_symbol", reformatted_headers
                )
                alias_index = get_header_index(
                    "alias_symbols", reformatted_headers
                )
                prev_index = get_header_index(
                    "previous_symbols", reformatted_headers
                )
                ensg_index = get_header_index(
                    "ext_ensembl_id", reformatted_headers
                )
                hgnc_index = get_header_index(
                    "hgnc_id", reformatted_headers
                )

            else:
                line = line.strip("\n").split("\t")

                for j, ele in enumerate(line):
                    if j == 0:
                        hgnc_id = ele
                        data.setdefault(hgnc_id, {})
                    else:
                        # we have the index of the line so we can automatically
                        # get the header and use it has a subkey in the dict
                        data[hgnc_id][reformatted_headers[j]] = ele

                alias_symbols = line[alias_index].split(",")
                prev_symbols = line[prev_index].split(",")

                symbol_dict[line[symbol_index]] = {
                    "hgnc_id": line[hgnc_index],
                    "ensg_id": line[ensg_index]
                }

                for symbol in alias_symbols:
                    alias_dict[symbol.strip()] = {
                        "hgnc_id": line[hgnc_index],
                        "ensg_id": line[ensg_index]
                    }

                for symbol in prev_symbols:
                    prev_dict[symbol.strip()] = {
                        "hgnc_id": line[hgnc_index],
                        "ensg_id": line[ensg_index]
                    }

    return data, symbol_dict, alias_dict, prev_dict


def find_id_using_ensg(ensg_id: str, hgnc_data: dict):
    """ Find HGNC id using the ENSG id
    Args:
        ensg_id (str): ENSG id
        hgnc_data (dict): Dict of HGNC data from HGNC dump
    Returns:
        str: HGNC id
    """

    for hgnc_id in hgnc_data:
        if ensg_id == hgnc_data[hgnc_id]["ext_ensembl_id"]:
            return hgnc_id

    return


def get_nirvana_data_dict(
    nirvana_gff: str, all_symbols: list, hgnc_data: dict, symbol_data: dict,
    alias_data: dict, prev_data: dict
):
    """ Return dict of parsed data for Nirvana
    Args:
        nirvana_refseq (str): GFF filepath for nirvana
        hgnc_data (dict): Dict of HGNC data
        alias_data (dict): Dict of alias HGNC data
        prev_data (dict): Dict of previous symbol HGNC data

    Returns:
        dict: Dict of gene2transcripts2exons
    """

    nirvana_tx_dict = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(None)
            )
        )
    )

    symbol2hgnc = {}

    with gzip.open(nirvana_gff) as nir_fh:
        for index, line in enumerate(nir_fh):
            fields = line.decode("utf-8").strip().split("\t")
            record_type = fields[2]
            info_field = fields[8]
            info_fields = info_field.split("; ")
            info_dict = {}

            # skip lines where the entity type is gene, UTR, CDS
            if record_type in ["UTR", "CDS"]:
                continue

            for field in info_fields:
                key, value = field.split(" ")
                value = value.strip("\"").strip("\'")
                info_dict[key] = value

            gff_gene_name = info_dict["gene_name"]

            if record_type == "gene":
                if "ensembl_gene_id" in info_dict:
                    hgnc_id = find_id_using_ensg(
                        info_dict["ensembl_gene_id"], hgnc_data
                    )

                    if hgnc_id is None:
                        hgnc_id, ensg_id = assign_ids_to_symbol(
                            gff_gene_name, all_symbols, symbol_data,
                            alias_data, prev_data
                        )
                else:
                    hgnc_id, ensg_id = assign_ids_to_symbol(
                        gff_gene_name, all_symbols, symbol_data,
                        alias_data, prev_data
                    )

                if hgnc_id is None:
                    continue
                else:
                    symbol2hgnc[gff_gene_name] = hgnc_id

                continue

            gff_transcript = info_dict["transcript_id"]

            if gff_gene_name not in symbol2hgnc:
                continue
            else:
                hgnc_id = symbol2hgnc[gff_gene_name]

            if record_type == "transcript":
                if "tag" in info_dict:
                    nirvana_tx_dict[hgnc_id][gff_transcript]["canonical"] = True
                else:
                    nirvana_tx_dict[hgnc_id][gff_transcript]["canonical"] = False

    return nirvana_tx_dict


def assign_ids_to_symbol(
    gene_symbol: str, all_symbols: list, hgnc_data: dict, alias_data: dict,
    prev_data: dict
):
    """ Get the hgnc and ensg ids from a gene symbol

    Args:
        gene_symbol (str): Gene symbol
        hgnc_data (dict): Dict of hgnc data
        alias_data (dict): Dict of alias hgnc data
        prev_data (dict): Dict of previous symbols hgnc data

    Returns:
        tuple: hgnc id and ensg id
    """

    hgnc_id = None
    ensg_id = None

    amount_symbol_in_data = len([
        symbol
        for symbol in all_symbols
        if symbol.upper() == gene_symbol.upper()
    ])

    if amount_symbol_in_data > 1:
        return f"{gene_symbol}_hgnc_id_TBD", f"{gene_symbol}_ensg_id_TBD"

    if gene_symbol in hgnc_data:
        hgnc_id = hgnc_data[gene_symbol]["hgnc_id"]
        ensg_id = hgnc_data[gene_symbol]["ensg_id"]
    elif gene_symbol in alias_data:
        hgnc_id = alias_data[gene_symbol]["hgnc_id"]
        ensg_id = alias_data[gene_symbol]["ensg_id"]
    elif gene_symbol in prev_data:
        hgnc_id = prev_data[gene_symbol]["hgnc_id"]
        ensg_id = prev_data[gene_symbol]["ensg_id"]

    return hgnc_id, ensg_id


def assign_transcript(
    session, meta, hgnc_id: str, nirvana_dict: dict
):
    """ Return transcript data and clinical transcript from hgmd/nirvana

    Args:
        session (SQLAlchemy session): SQLAlchemy session obj
        meta (SQLAlchemy meta): SQLAlchemy meta obj
        hgnc_id (str): HGNC id
        ensg_id (str): ENSG id
        nirvana_dict (dict): Dict of parsed data from gff nirvana

    Returns:
        tuple: Dict of transcript data and clinical transcript refseq
    """

    markgene_tb = meta.tables["markname"]
    gene2refseq_tb = meta.tables["gene2refseq"]

    transcript_data = {}
    clinical_transcript = None

    # get the hgmd transcripts using the HGNC id provided
    hgmd_transcripts = session.query(
        gene2refseq_tb.c.refcore, gene2refseq_tb.c.refversion
    ).join(
        markgene_tb, markgene_tb.c.gene_id == gene2refseq_tb.c.hgmdID
    ).filter(
        markgene_tb.c.hgncID == hgnc_id[5:]
    ).all()

    # check if the gene is in nirvana
    if hgnc_id in nirvana_dict:
        # get transcripts from the gff
        nirvana_data = nirvana_dict[hgnc_id]

        for nirvana_tx in nirvana_data:
            transcript_data.setdefault(nirvana_tx, {})
            base_nirvana, version_nirvana = nirvana_tx.split(".")

            # store all transcripts in the dict
            transcript_data[nirvana_tx]["match"] = None

            # add the canonical status
            if nirvana_data[nirvana_tx]["canonical"] is True:
                transcript_data[nirvana_tx]["canonical"] = True
            else:
                transcript_data[nirvana_tx]["canonical"] = False

            # if the gene is present in hgmd
            if hgmd_transcripts:
                for base_hgmd_tx, hgmx_tx_version in hgmd_transcripts:
                    # check if one of the hgmd transcripts is equal to one of
                    # the nirvana transcripts
                    if f"{base_hgmd_tx}.{hgmx_tx_version}" == nirvana_tx:
                        transcript_data[nirvana_tx]["match"] = "exact"
                    elif base_hgmd_tx == base_nirvana:
                        transcript_data[nirvana_tx]["match"] = "partial"
                    else:
                        # none of the hgmd transcripts match the nirvana ones
                        # assign the canonical transcript with the canonical
                        # match
                        if nirvana_data[nirvana_tx]["canonical"] is True:
                            transcript_data[nirvana_tx]["match"] = "canonical"

            else:
                # if it's not present in hgmd the canonical transcript because
                # the clinical transcript
                if nirvana_data[nirvana_tx]["canonical"] is True:
                    transcript_data[nirvana_tx]["match"] = "canonical"

        # okay assignment of clinical transcript
        # firstly get the various matches
        final_check = [transcript_data[tx]["match"] for tx in transcript_data]

        # if not exact match, find partial or canonical matches
        if "exact" not in final_check:
            # priority order: check if there a partial match before checking if
            # there is a canonical match
            if "partial" in final_check:
                partial_match_txs = [
                    tx
                    for tx in transcript_data
                    if transcript_data[tx]["match"] == "partial"
                ]

                # check how many partial matches are present
                if len(partial_match_txs) > 1:
                    # get the highest version number for choosing the clinical
                    # transcript
                    max_version = max(
                        [tx.split(".")[1]for tx in partial_match_txs]
                    )
                    clinical_transcript = [
                        tx
                        for tx in partial_match_txs
                        if tx.split(".")[1] == max_version
                    ][0]
                else:
                    clinical_transcript = partial_match_txs[0]
            elif "canonical" in final_check:
                clinical_transcript = [
                    tx
                    for tx in transcript_data
                    if transcript_data[tx]["match"] == "canonical"
                ][0]
                clinical_transcript = f"{clinical_transcript}_TO_REVIEW"
        # exact match so clear clinical transcript
        else:
            clinical_transcript = [
                tx
                for tx in transcript_data
                if transcript_data[tx]["match"] == "exact"
            ][0]

    return transcript_data, clinical_transcript


def get_header_index(header_name: str, headers: list):
    """ Get the index of a given header

    Args:
        header_name (str): Header name
        headers (list): List of headers

    Returns:
        int: Index for given header name
    """

    return headers.index(header_name)


def parse_test_directory(file: str):
    """ Parse the data in the National test directory

    Args:
        file (str): XLS of the National test directory

    Returns:
        tuple: Dict of clin_ind_id2clin_ind and dict of test_id2targets
    """

    clinind_data = defaultdict(lambda: defaultdict(str))

    xls = xlrd.open_workbook(file)
    sheet_with_tests = xls.sheet_by_name("R&ID indications")

    ci_dict = {}

    for row in range(sheet_with_tests.nrows):
        if row >= 2:
            (
                ci_id, ci, criteria, test_code,
                targets, method, clinical_group, comment
            ) = sheet_with_tests.row_values(row)

            if ci != "" and ci_id != "":
                ci_dict[ci_id] = ci
            else:
                ci_id, code = test_code.split(".")
                ci = ci_dict[ci_id]

            test_code = test_code.strip()

            if "panel" in method or "WES" in method or "Single gene" in method:
                clinind_data[test_code]["targets"] = targets.strip()
                clinind_data[test_code]["method"] = method.strip()
                clinind_data[test_code]["name"] = ci.strip()
                clinind_data[test_code]["version"] = file

    return clinind_data


def clean_targets(clinind_data: dict):
    """ Replace the methods from the XLS to abbreviation:
    WES -> P
    Panel -> P
    Single Gene -> G

    Args:
        clinind_data (dict): Dict of data from the test directory

    Returns:
        dict: Dict of dict for test2targets
    """

    clean_clinind_data = defaultdict(lambda: defaultdict(list))

    ci_to_remove = []

    for test_code in clinind_data:
        data = clinind_data[test_code]
        clinind = data["name"]
        targets = data["targets"]
        method = data["method"]
        version = data["version"]

        clean_clinind_data[test_code]["name"] = clinind
        clean_clinind_data[test_code]["version"] = version

        if "WES" in method:
            clean_clinind_data[test_code]["method"] = "P"

        elif "panel" in method:
            match = regex.search(r"(.*)[pP]anel", method)
            type_panel = match.groups()[0][0]
            clean_clinind_data[test_code]["method"] = f"{type_panel}P"

        elif "gene" in method:
            clean_clinind_data[test_code]["method"] = "G"

        cleaned_method = clean_clinind_data[test_code]["method"]

        clean_clinind_data[test_code]["gemini_name"] = (
            f"{test_code}_{clinind}_{cleaned_method}"
        )

        for indiv_target in targets.split(";"):
            indiv_target = indiv_target.strip()

            if "Relevant" not in indiv_target:
                # Panels can have "As dictated by blabla" "As indicated by"
                # so I remove those
                if indiv_target.startswith("As "):
                    ci_to_remove.append(test_code)

                # check if the target has parentheses with numbers in there
                match = regex.search(r"(?P<panel_id>\(\d+\))", indiv_target)

                # it's a panel, parentheses detected, really reliable
                if match:
                    target_to_add = match.group("panel_id").strip("()")
                    clean_clinind_data[test_code]["panels"].append(
                        target_to_add
                    )

                # it's a single gene
                else:
                    target_to_add = hq_get_id(
                        indiv_target.strip(), verbose=False
                    )

                    if target_to_add is not None:
                        clean_clinind_data[test_code]["genes"].append(
                            target_to_add
                        )
                    else:
                        # only case where this happens is a
                        # As dictated by clinical indication case
                        ci_to_remove.append(test_code)
            else:
                ci_to_remove.append(test_code)

        # convert default dict to dict because accessing absent key later on
        # will create the key with a list as default breaking what I'm doing
        # later
        clean_clinind_data[test_code] = dict(clean_clinind_data[test_code])

    # remove the clinical indication that have the Relevant panel/gene text
    for key in ci_to_remove:
        clean_clinind_data.pop(key, None)

    return clean_clinind_data


def gather_single_genes(clin_ind2targets: dict):
    """ Return list of all single genes tests

    Args:
        clin_ind2targets (dict): Dict of clinical indications

    Returns:
        list: List of single genes used in clinical indications
    """

    single_genes = []

    for test_code in clin_ind2targets:
        if "genes" in clin_ind2targets[test_code]:
            genes = clin_ind2targets[test_code]["genes"]
            single_genes += genes

    return single_genes


def write_new_output_folder(output_dump: str, output_suffix: str = ""):
    """ Return new folder to output files in

    Args:
        output_dump (str): Type of output folder
        output_suffix (str, optional): Suffix to be added to subfolder. Defaults to "".

    Returns:
        str: Folder path to the final output
    """

    output_date = get_date()
    output_index = 1

    output_folder = f"{output_dump}/{output_date}-{output_index}"

    if output_suffix:
        output_folder = f"{output_folder}_{output_suffix}"

    # don't want to overwrite files so create folders using the output index
    while Path(output_folder).is_dir():
        output_index += 1
        output_folder = f"{output_dump}/{output_date}-{output_index}"

        if output_suffix:
            output_folder = f"{output_folder}_{output_suffix}"

    Path(output_folder).mkdir(parents=True)

    return output_folder


def main(**args):
    if args["command"] == "g2t":
        print("Parsing and processing data")

        with open(args["gene_file"]) as f:
            genes = f.read().splitlines()

        hgnc_data, symbol_dict, alias_dict, prev_dict = parse_hgnc_dump(
            args["hgnc"]
        )
        all_symbols = (
            list(symbol_dict.keys()) + list(alias_dict.keys()) +
            list(prev_dict.keys())
        )

        transcript_data = get_nirvana_data_dict(
            args["gff"], all_symbols, hgnc_data, symbol_dict, alias_dict,
            prev_dict
        )

        session, meta = connect_to_db(
            "hgmd_ro", "hgmdreadonly", "localhost", args["database"]
        )

        folder = write_new_output_folder("genes2transcripts")
        logger = create_log(folder)
        print("Writing output file")

        with open(f"{folder}/{get_date()}_g2t.tsv", "w") as f:
            for gene in genes:
                if not gene.startswith("HGNC:"):
                    hgnc_id, ensg_id = assign_ids_to_symbol(
                        gene, all_symbols, symbol_dict, alias_dict, prev_dict
                    )

                    if hgnc_id.endswith("TBD"):
                        msg = (
                            f"{gene} has hgnc ids in multiple sources in HGNC "
                            "(main, alias, previous symbols)"
                        )
                        print(msg)
                        logger.warning(msg)
                        f.write(f"{gene}\t\t\n")
                else:
                    hgnc_id = gene

                tx_data, clinical_transcript = assign_transcript(
                    session, meta, hgnc_id, transcript_data
                )

                for tx in tx_data:
                    clinical_tx_status = None

                    if tx == clinical_transcript:
                        clinical_tx_status = "clinical_transcript"
                    else:
                        clinical_tx_status = "not_clinical_transcript"

                    if f"{tx}_TO_REVIEW" == clinical_transcript:
                        msg = (
                            f"{hgnc_id} - {tx}: Review the clinical transcript"
                        )
                        print(msg)
                        logger.warning(msg)
                        clinical_tx_status = "to_review"

                    if tx_data[tx]["canonical"] is True:
                        status = "canonical"
                    else:
                        status = "not_canonical"

                    f.write(
                        f"{hgnc_id}\t{tx}\t{clinical_tx_status}\t{status}\n"
                    )

        print(f"Written file is '{folder}/{get_date()}_g2t.tsv'")

    elif args["command"] == "gene_list":
        print("Parsing and processing data")

        genes = set()

        for folder in args["panel_folder"]:
            for file in os.listdir(folder):
                if "superpanel" not in file:
                    with open(f"{folder}/{file}") as f:
                        for line in f:
                            line = line.strip().split("\t")

                            if line[4] == "gene":
                                genes.add(line[-1])

        clinind_data = parse_test_directory(args["test_directory"])
        clean_clinind_data = clean_targets(clinind_data)

        # get the single genes in the test directory
        single_genes = gather_single_genes(clean_clinind_data)

        for gene in single_genes:
            genes.add(gene)

        folder = write_new_output_folder("gene_list")
        print("Writing output file")

        with open(f"{folder}/{get_date()}_genes", "w") as f:
            for gene in genes:
                f.write(f"{gene}\n")

        print(f"Written file is '{folder}/{get_date()}_genes.tsv'")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparser = parser.add_subparsers(dest="command")

    gene_file = subparser.add_parser(
        "gene_list", help=(
            "Generate file containing the hgnc ids of every gene in the "
            "panel folder given and in the test directory"
        )
    )
    gene_file.add_argument(
        "panel_folder", nargs="+", help=(
            "Folder containing the panel files to get all the genes"
        )
    )
    gene_file.add_argument("test_directory", help="National test directory")

    g2t = subparser.add_parser("g2t", help="Generate genes2transcripts file")
    g2t.add_argument("gene_file", help="Gene file")
    g2t.add_argument("hgnc", help="HGNC dump")
    g2t.add_argument("gff", help="Nirvana gff")
    g2t.add_argument("database", help="HGMD database to connect to")

    args = vars(parser.parse_args())
    main(**args)
