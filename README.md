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
