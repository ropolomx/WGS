import pandas as pd

configfile: "config.yaml"

SAMPLE = pd.read_table(config["samples"].set_index("samples", drop=False)

rule all:
    input:
        expand('{sample}/{sample}_prokka/{sample}.ffn', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.faa', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.fna', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.fsa', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.gbk', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.gff', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.sqn', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.tbl', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.tsv', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.txt', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.err', sample=SAMPLE),
	expand('{sample}/{sample}_prokka/{sample}.log', sample=SAMPLE),
        expand('{sample}/abricate_amr.txt', sample=SAMPLE),
        expand('{sample}/abricate_plasmid.txt', sample=SAMPLE),
        expand('{sample}/mlst.txt', sample=SAMPLE),
        #expand('{sample}/abricate_ecoh.txt', sample=SAMPLE),
        #expand('{sample}/abricate_ecolivf.txt', sample=SAMPLE),
        #expand('{sample}/ariba_virfinder/assemblies.fa.gz', sample=SAMPLE),
        #expand('{sample}/ariba_virfinder/report.tsv', sample=SAMPLE),
        #'ecoh_summary.txt',
        'amr_summary.txt',
        'vf_summary.txt',
        'plasmid_summary.txt'

rule trimmomatic:
    input:
        fwd=lambda wilcards: SAMPLE.loc[wilcards.sample, 'forward'],
        rev=lambda wilcards: SAMPLE.loc[wilcards.sample, 'reverse']

    output:
        fwd_paired='{sample}/{sample}_fwd_trimmed.fastq.gz',
        rev_paired='{sample}_rev_trimmed.fastq.gz'

    log:
        'logs/trimmomatic/{sample}.log'

    shell:
        '(trimmomatic PE -threads 12 -phred33 {input.fwd} {input.rev} {output.fwd_paired} /dev/null {output.rev_paired} /dev/null ILLUMINACLIP:/home/rodrigo/miniconda3/envs/wgs/share/trimmomatic-0.38-1/adapters/NexteraPE-PE.fa:1:25:11 LEADING:10 TRAILING:10 MINLEN:80) 2> {log}'

rule shovill:
    input:
        fwd='{sample}_fwd_trimmed.fastq.gz',
        rev='{sample}_rev_trimmed.fastq.gz'

    output:
        contigs='{sample}/contigs.fa',
        graph='{sample}/contigs.gfa'
    params:
        dir='{sample}'
    shell:
        'shovill --outdir {params.dir} --ram 32 --cpus 8 --R1 {input.fwd} --R2 {input.rev} --trim --force'

rule abricate_amr:
    input:
        '{sample}/{sample}_shovill/contigs.fa'
    output:
        '{sample}/abricate_amr.txt'
    log:
        'logs/abricate/{sample}.log'
    shell:
        '(abricate --threads 4 --mincov 60 --db ncbi {input} > {output}) 2> {log}'

rule abricate_plasmid:
    input:
        '{sample}/{sample}_shovill/contigs.fa'
    output:
        '{sample}/abricate_plasmid.txt'
    log:
        'logs/abricate/{sample}.log'
    shell:
        '(abricate --threads 4 --mincov 60 --db plasmidfinder {input} > {output}) 2> {log}'


#rule abricate_serotyping:
#    input:
#        '{sample}/{sample}_shovill/contigs.fa'
#    output:
#        '{sample}/abricate_ecoh.txt'
#    log:
#        'logs/abricate_serotyping/{sample}.log'
#    shell:
#       '(abricate --threads 4 --mincov 60 --db ecoh {input} > {output}) 2> {log}'

#rule abricate_ecolivf:
#    input:
#        '{sample}/{sample}_shovill/contigs.fa'
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
        '(abricate --threads 4 --mincov 60 --db ecoli_vf {input} > {output}) 2> {log}'

#rule abricate_ecoh_summary:
#    input:
#        expand('{sample}/abricate_ecoh.txt', sample=SAMPLE)
#    output:
#        'ecoh_summary.txt'
#    shell:
#        'abricate --summary {input} > {output}'

rule abricate_amr_summary:
    input:
        expand('{sample}/abricate_amr.txt', sample=SAMPLE)
    output:
        'amr_summary.txt'
    shell:
        'abricate --summary {input} > {output}'

rule abricate_vf_summary:
    input:
        expand('{sample}/abricate_ecolivf.txt', sample=SAMPLE)
    output:
        'vf_summary.txt'
    shell:
        'abricate --summary {input} > {output}'

rule abricate_plasmid_summary:
    input:
        expand('{sample}/abricate_plasmid.txt', sample=SAMPLE)
    output:
        'plasmid_summary.txt'
    shell:
        'abricate --summary {input} > {output}'

rule prokka:
    input:
        '{sample}/contigs.fa'
    output:
        function='{sample}/{sample}_prokka/{sample}.ffn',
	amino='{sample}/{sample}_prokka/{sample}.faa',
	nucl='{sample}/{sample}_prokka/{sample}.fna',
	fsa='{sample}/{sample}_prokka/{sample}.fsa',
	genbank='{sample}/{sample}_prokka/{sample}.gbk',
	gff='{sample}/{sample}_prokka/{sample}.gff',
	sqn='{sample}/{sample}_prokka/{sample}.sqn',
	tbl='{sample}/{sample}_prokka/{sample}.tbl',
	tsv='{sample}/{sample}_prokka/{sample}.tsv',
	summary='{sample}/{sample}_prokka/{sample}.txt',
	error='{sample}/{sample}_prokka/{sample}.err',
	log='{sample}/{sample}_prokka/{sample}.log'
    params:
        prefix='{sample}',
        outdir='{sample}/{sample}_prokka'
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
