import pandas as pd

configfile: 'config.yaml'

SAMPLE = pd.read_table(config["samples"]).set_index("isolate", drop=False)

rule all:
    input:
        expand('{sample}/prokka/{sample}.ffn', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.faa', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.fna', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.fsa', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.gbk', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.gff', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.sqn', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.tbl', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.tsv', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.txt', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.err', sample=SAMPLE['isolate']),
        expand('{sample}/prokka/{sample}.log', sample=SAMPLE['isolate']),
        expand('{sample}/abricate_amr.txt', sample=SAMPLE['isolate']),
        expand('{sample}/abricate_plasmid.txt', sample=SAMPLE['isolate']),
        expand('{sample}/mlst.txt', sample=SAMPLE['isolate']),
        #expand('{sample}/abricate_ecoh.txt', sample=SAMPLE['isolate']),
        #expand('{sample}/abricate_ecolivf.txt', sample=SAMPLE['isolate']),
        #expand('{sample}/ariba_virfinder/assemblies.fa.gz', sample=SAMPLE['isolate']),
        #expand('{sample}/ariba_virfinder/report.tsv', sample=SAMPLE['isolate']),
        #'ecoh_summary.txt',
        'amr_summary.txt',
        'vf_summary.txt',
        'plasmid_summary.txt',
        'roary/gene_presence_absence.csv'

rule shovill:
    input:
        forward=lambda wildcards: SAMPLE.loc[wildcards.sample, 'forward'],
        reverse=lambda wildcards: SAMPLE.loc[wildcards.sample, 'reverse']
    output:
        contigs='{sample}/contigs.fa',
        graph='{sample}/contigs.gfa'
    params:
        dir='{sample}'
    shell:
        'shovill --outdir {params.dir} --ram 32 --cpus 32 --R1 {input.forward} --R2 {input.reverse} --trim --force'

rule abricate_amr:
    input:
        '{sample}/contigs.fa'
    output:
        '{sample}/abricate_amr.txt'
    log:
        'logs/abricate/{sample}.log'
    shell:
        '(abricate --threads 4 --mincov 60 --db ncbi {input} > {output}) 2> {log}'

rule abricate_plasmid:
    input:
        '{sample}/contigs.fa'
    output:
        '{sample}/abricate_plasmid.txt'
    log:
        'logs/abricate/{sample}.log'
    shell:
        '(abricate --threads 4 --mincov 60 --db plasmidfinder {input} > {output}) 2> {log}'

#        '{sample}/contigs.fa'
#    output:
#        '{sample}/abricate_ecoh.txt'
#    log:
#        'logs/abricate_serotyping/{sample}.log'
#    shell:
#       '(abricate --threads 4 --mincov 60 --db ecoh {input} > {output}) 2> {log}'

#rule abricate_ecolivf:
#    input:
#        '{sample}/contigs.fa'
#    output:
#        '{sample}/abricate_ecolivf.txt'
#    log:
#        'logs/abricate_ecolivf/{sample}.log'
#    shell:
#        '(abricate --threads 4 --mincov 60 --db ecoli_vf {input} > {output}) 2> {log}'

rule abricate_vf:
    input:
        '{sample}/contigs.fa'
    output:
        '{sample}/abricate_vf.txt'
    log:
        'logs/abricate_vf/{sample}.log'
    shell:
        '(abricate --threads 4 --mincov 60 --db vfdb {input} > {output}) 2> {log}'

#rule abricate_ecoh_summary:
#    input:
#        expand('{sample}/abricate_ecoh.txt', sample=SAMPLE)
#    output:
#        'ecoh_summary.txt'
#    shell:
#        'abricate --summary {input} > {output}'

rule abricate_amr_summary:
    input:
        expand('{sample}/abricate_amr.txt', sample=SAMPLE.isolate)
    output:
        'amr_summary.txt'
    shell:
        'abricate --summary {input} > {output}'

rule abricate_vf_summary:
    input:
        expand('{sample}/abricate_vf.txt', sample=SAMPLE.isolate)
    output:
        'vf_summary.txt'
    shell:
        'abricate --summary {input} > {output}'

rule abricate_plasmid_summary:
    input:
        expand('{sample}/abricate_plasmid.txt', sample=SAMPLE.isolate)
    output:
        'plasmid_summary.txt'
    shell:
        'abricate --summary {input} > {output}'

rule prokka:
    input:
        '{sample}/contigs.fa'
    output:
        function='{sample}/prokka/{sample}.ffn',
        amino='{sample}/prokka/{sample}.faa',
        nucl='{sample}/prokka/{sample}.fna',
        fsa='{sample}/prokka/{sample}.fsa',
        genbank='{sample}/prokka/{sample}.gbk',
        gff='{sample}/prokka/{sample}.gff',
        sqn='{sample}/prokka/{sample}.sqn',
        tbl='{sample}/prokka/{sample}.tbl',
        tsv='{sample}/prokka/{sample}.tsv',
        summary='{sample}/prokka/{sample}.txt',
        error='{sample}/prokka/{sample}.err',
        log_file='{sample}/prokka/{sample}.log'
    params:
        prefix='{sample}',
        outdir='{sample}/prokka'
    shell:
        'prokka --kingdom Bacteria --prefix {params.prefix} --outdir {params.outdir} --cpus 32 --force {input}'

#rule ariba_run:
#    input:
#        fwd='{sample}_fwd_trimmed.fastq.gz',
#        rev='{sample}_rev_trimmed.fastq.gz'
#    output:
#        gz='{sample}/ariba_virfinder/assemblies.fa.gz',
#        rep='{sample}/ariba_virfinder/report.tsv'
#    params:
#        ariba_db = '/media/amr_storage/ariba_ref/virfinder',
#        ariba_out = '{sample}/ariba_virfinder'
#    shell:
#        'ariba run --threads 12 {params.ariba_db} {input.fwd} {input.rev} --force {params.ariba_out}'

rule mlst:
    input:
        '{sample}/contigs.fa'
    output:
        '{sample}/mlst.txt'
    shell:
        'mlst --threads 4 {input} > {output}'

rule pangenome:
    input:
        expand('{sample}/prokka/{sample}.gff', sample=SAMPLE.isolate)
    output:
        pan='roary/gene_presence_absence.csv'
    shell:
        'roary -f roary -p 32 -t 11 {input}'
