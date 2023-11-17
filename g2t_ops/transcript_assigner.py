from concurrent.futures import ThreadPoolExecutor, as_completed
import datetime
import os
import gffutils
import sqlite3
from g2t_ops import utils


def get_date_for_db() -> str:
    """ Return date in YYYY-MM-DD format for insertion in database

    Returns:
        str: String with date
    """

    return datetime.datetime.now().strftime("%Y-%m-%d")


def parse_mane_file(mane_file, hgnc_dump) -> dict:
    """ Parse a MANE CSV file downloaded from
    http://tark.ensembl.org/web/mane_GRCh37_list/

    Args:
        mane_file (str): Path to MANE CSV file
        hgnc_dump (dict): Dict containing the HGNC data

    Returns:
        dict: Dict containing MANE data
    """

    mane_data = {}

    with open(mane_file) as f:
        with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
            future_res = {
                executor.submit(
                    find_hgnc_id_for_mane, line, hgnc_dump
                ): line for line in f
            }

            for res in as_completed(future_res):
                if res.result():
                    mane_data = {**mane_data, **res.result()}

    return mane_data


def find_hgnc_id_for_mane(line, hgnc_dump) -> dict:

    # each column has quotes around them so i'm striping them i.e.
    # "A1BG","MANE SELECT","ENST00000263100.8","NM_130786.4"
    cleaned_data = [ele.strip("\"") for ele in line.strip().split(",")]
    gene_symbol, mane_type, ensembl_38, refseq = cleaned_data[0:4]

    if mane_type == "MANE SELECT":
        hgnc_data = utils.find_hgnc_id(gene_symbol, hgnc_dump)
        hgnc_id = hgnc_data[1]

        if hgnc_id:
            return {hgnc_id: refseq}
        else:
            return None


def parse_mane_gff(gff):
    """ Parse the gff data

    Args:
        gff (str): Path to GFF

    Returns:
        FeatureDB: FeatureDB object for the gff
    """

    # try to create sqlite db
    try:
        db = gffutils.create_db(
            gff, "MANE_refseq.sqlite", verbose=True,
            merge_strategy="create_unique"
        )
    except sqlite3.OperationalError as e:
        # use existing db
        db = gffutils.FeatureDB("MANE_refseq.sqlite")

    return db


def get_mane_transcripts_from_b38_gff(db):
    """
    Make dictionary with each transcript ID in the GFF and its HGNC ID
    Args:
        FeatureDB: FeatureDB object for the gff
    Returns:
        mane_data (dict): a dictionary in the format {hgnc_id: refseq}
    """
    mane_data = {}
    # For each exon in the gff
    for feature in db.features_of_type("exon"):
        # Extract HGNC ID
        hgnc_list = [
            i
            for i in feature.attributes["Dbxref"]
            if "HGNC" in i
        ]

        if hgnc_list:
            hgnc_id = hgnc_list[0].split(":", 1)[-1]
        else:
            hgnc_id = "None provided"

        # Extract MANE tag (MANE Select / MANE Plus Clinical) and transcript ID
        mane_tag = feature.attributes["tag"][0]
        transcript_id = feature.attributes["transcript_id"][0]

        # Add unique transcript IDs and MANE status to MANE data dict
        if hgnc_id in mane_data.keys():
            if transcript_id not in mane_data[hgnc_id].keys():
                mane_data[hgnc_id][transcript_id] = mane_tag
        else:
            mane_data[hgnc_id] = {transcript_id: mane_tag}
    print(
        f" There are {len(list(mane_data.values()))} MANE transcripts in the"
        " Refseq GFF"
    )
    return mane_data


def find_HGMD_transcript(session, meta, hgnc_id) -> str:
    """ Find the HGMD transcript using a HGNC id

    Args:
        session (session object): Session SQLAlchemy object
        meta (meta object): Meta SQLAlchemy object
        hgnc_id (str): HGNC id

    Raises:
        Exception: if 2 or more transcripts are found for the given HGNC id

    Returns:
        str: HGMD transcript
    """

    markgene_tb = meta.tables["markname"]
    gene2refseq_tb = meta.tables["gene2refseq"]

    # get the hgmd transcripts using the HGNC id provided
    hgmd_transcripts = session.query(
        gene2refseq_tb.c.refcore, gene2refseq_tb.c.refversion
    ).join(
        markgene_tb, markgene_tb.c.gene_id == gene2refseq_tb.c.hgmdID
    ).filter(markgene_tb.c.hgncID == hgnc_id[5:]).all()

    hgmd_transcript = [
        f"{hgmd_base}.{hgmd_version}"
        for hgmd_base, hgmd_version in hgmd_transcripts
    ]

    if len(hgmd_transcript) == 0:
        return
    elif len(hgmd_transcript) == 1:
        return hgmd_transcript[0]
    else:
        # HGNC ids should only be linked to one HGMD transcript
        raise Exception(
            f"{hgnc_id} has the following transcripts in the HGMD database: "
            f"{', '.join(hgmd_transcript)}"
        )


def get_transcripts(hgnc_ids, exon_data) -> dict:
    data = {}

    for hgnc_id in hgnc_ids:
        if hgnc_id in exon_data:
            data[hgnc_id] = exon_data[hgnc_id]
        else:
            print(f"{hgnc_id} is not present in the GFF")

    return data


def label_mane_transcripts(data, gene, mane_transcripts, tx, source):
    '''
    Label MANE transcripts
    inputs:
        data: dictionary of gene and transcripts
        gene: HGNC ID of gene for query transcript
        mane_transcripts: MANE transcript for gene (can be list or dictionary)
        depending on source of MANE transcript data
        tx: query transcript
        source: source of MANE transcript data, either "refseq_gff" or
        "ensembl_csv"
    outputs:
        data: dictionary of gene and transcripts with HGMD labels
    '''
    tx_base, tx_version = tx.split(".")

    if source == "refseq_gff":
        mane_transcripts_keys = list(mane_transcripts.keys())

        for mane_transcript in mane_transcripts_keys:
            mane_base, mane_version = mane_transcript.split(".")

            # Compare g2t transcript with MANE transcript and add
            # 'clinical transcript' label if they match. Multiple
            # clinical transcripts are possible as MANE can have
            # both a Select and Plus Clinical for the same gene
            if tx_base == mane_base:
                mane_status = mane_transcripts[mane_transcript]
                if "clinical_transcript" in data[gene]:
                    data[gene]["clinical_transcript"].append(
                        [tx, mane_status]
                    )
                else:
                    data[gene]["clinical_transcript"] = [[tx, mane_status]]

    if source == "ensembl_csv":
        mane_base, mane_version = mane_transcripts.split(".")

        # compare transcripts without the versions
        if tx_base == mane_base:
            data[gene]["clinical_transcript"] = [[tx, "MANE"]]

    return data


def label_hgmd_transcripts(data, gene, hgmd_transcript, tx):
    '''
    Label HGMD transcripts if not previously labelled as MANE
    inputs:
        data: dictionary of gene and transcripts
        gene: HGNC ID of gene for query transcript
        hgmd_transcript: HGMD transcript for gene
        tx: query transcript
    outputs:
        data: dictionary of gene and transcripts with HGMD labels
    '''
    tx_base, tx_version = tx.split(".")
    hgmd_base, hgmd_version = hgmd_transcript.split(".")
    if tx_base == hgmd_base:
        if "clinical_transcript" not in data[gene]:
            # This HGMD transcript has already been labelled as
            # MANE, no need to duplicate it
            data[gene]["clinical_transcript"] = [[tx, "HGMD"]]
    return data


def assign_transcripts(session, meta, mane_select_data, g2t_data, source) -> dict:
    """ Assign a clinical transcript status for all the genes in the g2t file

    Args:
        session (session object): Session SQLAlchemy object
        meta (meta object): Meta SQLAlchemy object
        mane_select_data (dict): Dict containing the MANE select data
        g2t_data (dict): Dict containing the g2t data

    Returns:
        dict: Dict of dict with the gene as key, status as subkey and the
        transcripts as values
    """

    data = {}
    print("There are " + str(len(list(g2t_data.values()))))
    for gene, transcripts in g2t_data.items():
        data.setdefault(gene, {})
        data[gene].setdefault("no_clinical_transcript", [])

        for tx in transcripts:
            # Check if tx is MANE
            if gene in mane_select_data:
                mane_transcripts = mane_select_data[gene]
                data = label_mane_transcripts(
                    data, gene, mane_transcripts, tx, source
                    )

            # Check if tx is HGMD
            hgmd_transcript = find_HGMD_transcript(session, meta, gene)
            if hgmd_transcript:
                data = label_hgmd_transcripts(data, gene, hgmd_transcript, tx)

            # Check the query transcript is not already listed. If not,
            # add it as a non-clinical transcript
            if "clinical_transcript" in data[gene]:
                if any(tx in sublist for sublist in data[gene]["clinical_transcript"]) is False:
                    data[gene]["no_clinical_transcript"].append([tx, "None"])
            else:
                data[gene]["no_clinical_transcript"].append([tx, "None"])

    return data


def write_g2t(data, output_path):
    """ Write g2t style file to know which transcript have been assigned to
    which genes

    Args:
        data (dict): Dict of dict with the gene as key, status as subkey and
        the transcripts as values
        output_path (str): Path to output
    """

    with open(f"{output_path}/g2t_file.tsv", "w") as f:
        for gene in data:
            for status, txs in data[gene].items():
                for one_tx in txs:
                    tx, source = one_tx
                    f.write(f"{gene}\t{tx}\t{status}\t{source}\n")


def write_sql_queries(data, output_path):
    """ Write SQL queries to run on the panel database to add genes and
    transcripts output

    Args:
        data (dict): Dict of dict with the gene as key, status as subkey and
        the transcripts as values
        output_path (str): Path to output
    """

    with open(f"{output_path}/sql_queries_to_import.sql", "w") as f:
        for gene in data:
            f.write(f"INSERT INTO gene (hgnc_id) VALUES (\"{gene}\");\n")
            f.write("SET @gene_id = (SELECT LAST_INSERT_ID());\n")
            f.write(
                f"INSERT INTO feature (gene_id, feature_type_id) "
                "VALUES (@gene_id, 1);\n"
            )

            for status, txs in data[gene].items():
                for tx in txs:
                    tx_base, tx_version = tx[0].split(".")
                    f.write((
                        "INSERT INTO transcript (refseq_base, version, canonical) "
                        f"VALUES (\"{tx_base}\", \"{tx_version}\", 0);\n"
                    ))
                    f.write("SET @transcript_id = (SELECT LAST_INSERT_ID());\n")
                    f.write((
                        "INSERT INTO genes2transcripts "
                        "(clinical_transcript, date, gene_id, reference_id, transcript_id) "
                        f"VALUES (0, \"{get_date_for_db()}\", @gene_id, 1, @transcript_id);\n"
                    ))


def write_transcript_status(data, output_path):
    """ Write a file with genes and their clinical status

    Args:
        data (dict): Dict of dict with the gene as key, status as subkey and
        the transcripts as values
        output_path (str): Path to output
    """

    with open(f"{output_path}/gene_transcript_status.txt", "w") as f:
        for gene, transcript_data in data.items():
            # check if a gene has "clinical_transcript" subkey
            all_clinical_transcripts = any([
                True if status == "clinical_transcript" else False
                for status in transcript_data
            ])

            if all_clinical_transcripts:
                f.write(f"{gene}\tTrue\n")
            else:
                f.write(f"{gene}\tFalse\n")
