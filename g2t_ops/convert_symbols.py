from g2t_ops import utils


def convert_symbols(symbols: list, hgnc_dump):
    """ Convert list of gene symbols into HGNC ids

    Args:
        symbols (list): List of gene symbols to convert
        hgnc_dump (pd.Dataframe): HGNC dataframe

    Returns:
        list: List for which each element is composed of gene symbol, HGNC id
        and status
    """

    hgnc_ids = []

    for symbol in symbols:
        hgnc_id = utils.find_hgnc_id(symbol, hgnc_dump)
        hgnc_ids.append(hgnc_id)

    return hgnc_ids


def write_hgnc_ids(data, output_path):
    """ Write gene symbols and their HGNC id

    Args:
        data (list): List for which each element is composed of gene symbol,
        HGNC id and status
        output_path (str): Output path folder
    """

    with open(f"{output_path}/hgnc_ids.txt", "w") as f:
        for symbol, hgnc_id, status in data:
            f.write(f"{symbol}\t{hgnc_id}\t{status}\n")
