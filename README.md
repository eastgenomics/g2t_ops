# g2t_ops

## Prerequisites

To run this code, you need to have:

- A MANE Select file downloaded from Ensembl http://tark.ensembl.org/web/mane_GRCh37_list/ (for build 37), or a MANE Select file from RefSeq https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/ (for build 38)
- A HGNC dump downloaded from https://www.genenames.org/download/custom/ (default columns). This is not needed if a RefSeq MANE gff is being used, because the HGNC IDs are included in the gff.
- A local HGMD database

DNAnexus has some HGMD dumps in dev projects. To setup the HGMD database:

- Download the dump (.dump.gz)
- Unzip it
- Import it

Commands to do it:

```bash
dx download ${hgmd_dump}
gunzip ${hgmd_dump}

mysql
# in mysql environment
mysql> create database ${hgmd_database_name};
# exit mysql environment with ctrl + D

mysql -D ${hgmd_database_name} < ${hgmd_dump}
```

Create a user for the HGMD database:

```sql
CREATE USER IF NOT EXISTS 'hgmd_ro'@'localhost' IDENTIFIED BY 'hgmdreadonly';
GRANT SELECT ON ${hgmd_database_name} . * TO 'hgmd_ro'@'localhost';
```

Setup environment

```bash
python3 -m venv ${path_to_env}
pip install -r requirements.txt
```

## How to run

The symbol file and hgnc files are single column files with a single element per row.

```
BRCA1
BRCA2
TAZ
TENS1
```

```
HGNC:1100
HGNC:1101
HGNC:12027s
```

```bash
source ${path_to_env}/bin/activate

python main.py --hgnc_dump ${hgnc_dump} convert_symbols -symbols ${symbol} ${symbol} ...
python main.py --hgnc_dump ${hgnc_dump} convert_symbols -symbol_file ${symbol_file}

# for build 37, Ensembl MANE Select csv 
python main.py --hgnc_dump ${hgnc_dump} assign_transcripts -mane ${Ensembl_MANE_Select_file} -hgmd ${database_name} ${database_usr} ${database_pwd} -hgnc_ids ${hgnc_id} ${hgnc_id} ...
python main.py --hgnc_dump ${hgnc_dump} assign_transcripts -mane ${Ensembl_MANE_Select_file} -hgmd ${database_name} ${database_usr} ${database_pwd} -hgnc_file ${hgnc_file}

# for build 38, RefSeq MANE gff
python main.py assign_transcripts --mane_gff ${MANE_RefSeq_gff} -hgmd ${database_name} ${database_usr} ${database_pwd} -hgnc_ids ${hgnc_id} ${hgnc_id} ...
python main.py assign_transcripts --mane_gff ${MANE_RefSeq_gff} -hgmd ${database_name} ${database_usr} ${database_pwd} -hgnc_file ${hgnc_file}
```

By default, the output created by this script will be located `./YYMMDD_results`. If that folder already exists, a new folder will be created with an number to differentiate folders.
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
These columns represent:
HGNC ID | Transcript | Clinical transcript status (i.e. could the code find a match in MANE or HGMD) | Source of the clinical transcript (i.e. MANE or HGMD)
-- | -- | -- | --

For every gene passed, its transcripts gathered from either a given exon file or a default exon file in 001_Reference. The `status` column will inform whether the transcript is the clinical transcript or not. And the `source` column informs how the clinical transcript status was selected (MANE/HGMD).

The second file contains the SQL queries needed to import the genes and the transcripts.

The third file contains the genes and if a clinical transcript was found for it.
