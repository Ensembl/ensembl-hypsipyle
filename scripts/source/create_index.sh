#!/usr/bin/env bash
set -euo pipefail

VCF_DIR="${1:-.}"
DB_FILE="${2:-variants.duckdb}"

if [ ! -d "$VCF_DIR" ]; then
  echo "Error: VCF directory not found: $VCF_DIR" >&2
  exit 1
fi

find "$VCF_DIR" -maxdepth 1 -type f -name "*.vcf.gz" | parallel -j 8 \
  "bcftools query -f '%CHROM\t%POS\t%ID\t%ALT\n' {} | awk -v f='{}' 'BEGIN{OFS=\"\t\"; n=split(f,p,\"/\"); base=p[n]} {print \$1, \$2, \$3, \$4, base}'" \
| duckdb "$DB_FILE" -c "
CREATE OR REPLACE TABLE variants AS
SELECT *
FROM read_csv('/dev/stdin',
    delim='\t',
    columns={
        'chr': 'VARCHAR',
        'start': 'INTEGER',
        'id': 'VARCHAR',
        'allele': 'VARCHAR',
        'source_file': 'VARCHAR'
    }
)
ORDER BY chr, start, id;

CREATE INDEX IF NOT EXISTS idx_variants_id ON variants(id);
CREATE INDEX IF NOT EXISTS idx_variants_chr_start ON variants(chr, start);
"


