# generate_g2t

Need to have HGMD database with those credentials: "hgmd_ro", "hgmdreadonly", "localhost", "hgmd_2020_3"

## Import hgmd database

```bash
dx download project-Fz4Q15Q42Z9YjYk110b3vGYQ:file-Fz4Q46842Z9z2Q6ZBjy7jVPY
gunzip hgmd_pro-2020.3.dump.gz
mysql -u ${user} -p < hgmd_pro-2020.3.dump
```

```sql
CREATE IF NOT EXISTS USER 'hgmd_ro'@'localhost' IDENTIFIED BY 'hgmdreadonly';
GRANT SELECT ON hgmd_2020_3 . * TO 'hgmd_ro'@'localhost';
```

## How to run

```bash
# Generate gene file
python g2t.py gene_list ${panel_folder} ${test_directory_xls}

# Generate g2t file from the gene file previously generated
python g2t.py g2t ${hgnc_file} ${gff} ${hgmd_database_name} -gf ${gene_file}

# Generate g2t file from panel form and older g2t
python g2t.py g2t ${hgnc_file} ${gff} ${hgmd_database_name} -pf ${panel_form} -g2t ${g2t_file}
```

The script also outputs a log file with each file generated in the output folder: `genes2transcripts/${date}-${index}/` and `gene_list/${date}-${index}/`
