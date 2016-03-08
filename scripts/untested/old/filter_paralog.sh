bgzip $1
tabix -p vcf $1.gz
vcf-annotate --filter c=$2,$3 $1.gz | mawk '!/SnpC/' > $4.vcf
vcf-annotate --filter c=$2,$3 $1.gz | mawk '!/#/' | mawk '/SnpC/' | wc -l
gunzip $1
