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

```python
python g2t.py gene_file --hgnc hgnc_dump --nirvana nirvana_gff
```
