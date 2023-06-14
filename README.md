# g2t_ops

## Prerequisites

To run this code, you need to have:

- A HGNC dump downloaded from https://www.genenames.org/download/custom/
- A MANE Select file downloaded from http://dev-tark.ensembl.org/web/mane_GRCh37_list/
- A local HGMD database

DNAnexus has some HGMD dumps in dev projects. To setup the HGMD database:

- Download the dump (.dump.gz)
- Unzip it
- Import it

Commands to do it:

```bash
dx download ${hgmd_dump}
gunzip ${hgmd_dump}
mysql -u ${user} -p < ${unzipped_hgmd_dump}
```

Create a user for the HGMD database:

```sql
CREATE IF NOT EXISTS USER 'hgmd_ro'@'localhost' IDENTIFIED BY 'hgmdreadonly';
GRANT SELECT ON hgmd_2020_3 . * TO 'hgmd_ro'@'localhost';
```

Setup environment

```bash
python3 -m venv ${path_to_env}
pip install -r requirements.txt
```

## How to run

```bash
source ${path_to_env}/bin/activate

python main.py --hgnc_dump ${hgnc_dump} convert_symbols -symbols ${symbol} ${symbol} ...
python main.py --hgnc_dump ${hgnc_dump} convert_symbols -symbol_file ${symbol_file}

python main.py --hgnc_dump ${hgnc_dump} assign_transcripts -mane ${MANE_Select_file} -hgmd ${database_name} ${database_usr} ${database_pwd} -hgnc_ids ${hgnc_id} ${hgnc_id} ...
python main.py --hgnc_dump ${hgnc_dump} assign_transcripts -mane ${MANE_Select_file} -hgmd ${database_name} ${database_usr} ${database_pwd} -hgnc_file ${hgnc_file}
```

By default, the output created by this script will be located `./YYMMDD_results`. If that folder already exists, a new folder will be created with an number to differenciate folders.
Otherwise a custom output folder can be specified like so:

```bash
python main.py -o ${output_folder} --hgnc_dump ${hgnc_dump} ...
```

## Outputs

The `convert_symbols` command will output one file called `hgnc_ids.txt` with the following format:

```
gene_symbol hgnc_id outcome
```

The `gene_symbol` column will contain all the symbols passed to g2t_ops.
The `hgnc_id` column will contain either the HGNC id found or `None` depending on the outcome of the search in the HGNC dump.
The `outcome` column will contain the outcome of the search i.e. in which column the HGNC id was found (Main, Previous, Alias) or if it is not, why it wasn't found (Multiple Previous, Multiple Alias).

The `transcript_assigner` command will output 3 files called:

- g2t_file.tsv
- sql_queries_to_import.sql
- gene_transcript_status.txt

The first is formatted as such:

```
gene    transcript  status  source
```

For every gene passed, its transcripts gathered from either a given exon file or a default exon file in 001_Reference. The `status` column will inform whether the transcript is the clinical transcript or not. And the `source` column informs how the clinical transcript status was selected (MANE/HGMD).

The second contains the SQL queries needed to import the genes and the transcripts.

The third file contains the genes and if a clinical transcript was found for it.
