import argparse
from collections import defaultdict
import gzip

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.schema import MetaData


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


def find_id_using_ensg(ensg_id, hgnc_data):
    for hgnc_id in hgnc_data:
        if ensg_id == hgnc_data[hgnc_id]["ext_ensembl_id"]:
            return hgnc_id

    return


def get_nirvana_data_dict(
    nirvana_gff: str, hgnc_data: dict, symbol_data: dict, alias_data: dict,
    prev_data: dict
):
    """ Return dict of parsed data for Nirvana
    Args:
        nirvana_refseq (str): GFF file for nirvana
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
                    hgnc_id = find_id_using_ensg(info_dict["ensembl_gene_id"], hgnc_data)

                    if hgnc_id is None:
                        hgnc_id, ensg_id = assign_ids_to_symbol(
                            gff_gene_name, symbol_data, alias_data, prev_data
                        )
                else:
                    hgnc_id, ensg_id = assign_ids_to_symbol(
                        gff_gene_name, symbol_data, alias_data, prev_data
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
    gene_symbol: str, hgnc_data: dict, alias_data: dict, prev_data: dict
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
                # if it's not present in hgmd the canonical transcript because
                # the clinical transcript
                if nirvana_data[nirvana_tx]["canonical"] is True:
                    transcript_data[nirvana_tx]["match"] = "canonical"

        # okay assignment of clinical transcript
        # firstly get the various matches
        final_check = [transcript_data[tx]["match"] for tx in transcript_data]

        # if not exact match, find partial or canonical matches 
        if "exact" not in final_check:
            if "partial" in final_check:
                clinical_transcript = [
                    tx
                    for tx in transcript_data
                    if transcript_data[tx]["match"] == "partial"
                ][0]
            elif "canonical" in final_check:
                clinical_transcript = [
                    tx
                    for tx in transcript_data
                    if transcript_data[tx]["match"] == "canonical"
                ][0]
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

    return [
        i
        for i, ele in enumerate(headers)
        if ele == header_name
    ][0]


def main(args):
    with open(args.gene_file) as f:
        genes = f.read().splitlines()

    hgnc_data, symbol_dict, alias_dict, prev_dict = parse_hgnc_dump(args.hgnc)
    transcript_data = get_nirvana_data_dict(
        args.gff, hgnc_data, symbol_dict, alias_dict, prev_dict
    )

    session, meta = connect_to_db(
        "hgmd_ro", "hgmdreadonly", "localhost", "hgmd_2020_3"
    )

    for gene in genes:
        hgnc_id, ensg_id = assign_ids_to_symbol(
            gene, symbol_dict, alias_dict, prev_dict
        )
        tx_data, clinical_transcript = assign_transcript(
            session, meta, hgnc_id, transcript_data
        )

        print(f"{gene}\t{clinical_transcript}")

        # for tx in data:
        #     print(f"{gene}\t{tx}\t{data[tx]['match']}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gene_file", help="Gene file")
    parser.add_argument("--hgnc", help="HGNC dump")
    parser.add_argument("--gff", help="Nirvana gff")

    args = parser.parse_args()
    main(args)
