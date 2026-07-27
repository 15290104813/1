
/data2/marui/pipeline/SMR/smr-1.3.1-linux-x86_64/smr --beqtl-summary /data2/marui/pipeline/SMR/00.eQTL/cis-eQTL-SMR/eQTLGen --query 5.0e-8 --genes /data2/marui/project/11.KXF2025011502/04_SMR/01.eqtl.list --out /data2/marui/project/11.KXF2025011502/04_SMR/besd_eqtl/multiple --make-besd

/data2/marui/pipeline/SMR/smr-1.3.1-linux-x86_64/smr --bfile /data2/marui/pipeline/SMR/00.g1000_eur/g1000_eur --gwas-summary /data2/marui/project/11.KXF2025011502/00_rawdata/00.outcome.ma --beqtl-summary /data2/marui/project/11.KXF2025011502/04_SMR/besd_eqtl/multiple --out /data2/marui/project/11.KXF2025011502/04_SMR/results_eqtl/multiple --thread-num 10 

/data2/marui/pipeline/SMR/smr-1.3.1-linux-x86_64/smr --bfile /data2/marui/pipeline/SMR/00.g1000_eur/g1000_eur --gwas-summary /data2/marui/project/11.KXF2025011502/00_rawdata/00.outcome.ma --beqtl-summary /data2/marui/pipeline/SMR/00.eQTL/cis-eQTL-SMR/eQTLGen --out /data2/marui/project/11.KXF2025011502/04_SMR/eqtl --plot --probe ENSG00000124813 --probe-wind 500 --gene-list /data2/marui/pipeline/SMR/00.glist_hg19_refseq.txt

/data2/marui/pipeline/SMR/smr-1.3.1-linux-x86_64/smr --bfile /data2/marui/pipeline/SMR/00.g1000_eur/g1000_eur --gwas-summary /data2/marui/project/11.KXF2025011502/00_rawdata/00.outcome.ma --beqtl-summary /data2/marui/pipeline/SMR/00.eQTL/cis-eQTL-SMR/eQTLGen --out /data2/marui/project/11.KXF2025011502/04_SMR/eqtl --plot --probe ENSG00000204217 --probe-wind 500 --gene-list /data2/marui/pipeline/SMR/00.glist_hg19_refseq.txt
