configfile: "config.yaml"

SAMPLE = pd.read_table(config["samples"].set_index("samples", drop=False)

FWD = expand('{sample}/{sample}_unmapped_fwd.fastq', sample=SAMPLE)
REV = expand('{sample}/{sample}_unmapped_rev.fastq', sample=SAMPLE)

rule all:
    input:
        expand('{sample}/{sample}_prokka/{sample}.ffn', sample=SAMPLE),
        expand('{sample}/abricate_amr.txt', sample=SAMPLE),
        expand('{sample}/abricate_plasmid.txt', sample=SAMPLE),
        expand('{sample}/abricate_ecoh.txt', sample=SAMPLE),
        expand('{sample}/abricate_ecolivf.txt', sample=SAMPLE),
        expand('{sample}/ariba_virfinder/assemblies.fa.gz', sample=SAMPLE),
        expand('{sample}/ariba_virfinder/report.tsv', sample=SAMPLE),
        'ecoh_summary.txt',
        'amr_summary.txt',
        'vf_summary.txt',
        'plasmid_summary.txt'

rule trimmomatic:
    input:
        fwd=lambda wilcards: SAMPLE.loc[wilcards.sample, 'forward'],
        rev=lambda wilcards: SAMPLE.loc[wilcards.sample, 'reverse'],

    output:
        fwd_paired='{sample}_fwd_trimmed.fastq.gz',
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
        contigs='{sample}/{sample}_shovill/contigs.fa',
        graph='{sample}/{sample}_shovill/contigs.gfa'
    params:
        dir='{sample}/{sample}_shovill'
    shell:
        'shovill --outdir {params.dir} --ram 22 --cpus 12 --R1 {input.fwd} --R2 {input.rev} --force'

rule abricate_amr:
    input:
        '{sample}/{sample}_shovill/contigs.fa'
    output:
        '{sample}/abricate_amr.txt'
    log:
        'logs/abricate/{sample}.log'
    shell:
        '(abricate --threads 4 --mincov 60 --db card {input} > {output}) 2> {log}'

rule abricate_plasmid:
    input:
        '{sample}/{sample}_shovill/contigs.fa'
    output:
        '{sample}/abricate_plasmid.txt'
    log:
        'logs/abricate/{sample}.log'
    shell:
        '(abricate --threads 4 --mincov 60 --db plasmidfinder {input} > {output}) 2> {log}'


rule abricate_serotyping:
    input:
        '{sample}/{sample}_shovill/contigs.fa'
    output:
        '{sample}/abricate_ecoh.txt'
    log:
        'logs/abricate_serotyping/{sample}.log'
    shell:
       '(abricate --threads 4 --mincov 60 --db ecoh {input} > {output}) 2> {log}'

rule abricate_ecolivf:
    input:
        '{sample}/{sample}_shovill/contigs.fa'
    output:
        '{sample}/abricate_ecolivf.txt'
    log:
        'logs/abricate_ecolivf/{sample}.log'
    shell:
        '(abricate --threads 4 --mincov 60 --db ecoli_vf {input} > {output}) 2> {log}'

rule abricate_ecoh_summary:
    input:
        expand('{sample}/abricate_ecoh.txt', sample=SAMPLE)
    output:
        'ecoh_summary.txt'
    shell:
        'abricate --summary {input} > {output}'

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
        '{sample}/{sample}_shovill/contigs.fa'
    output:
        '{sample}/{sample}_prokka/{sample}.ffn'
    params:
        prefix='{sample}',
        outdir='{sample}/{sample}_prokka'
    shell:
        'prokka --kingdom Bacteria --prefix {params.prefix} --outdir {params.outdir} --cpus 12 --force {input}'

rule ariba_run:
    input:
        fwd='{sample}_fwd_trimmed.fastq.gz',
        rev='{sample}_rev_trimmed.fastq.gz'
    output:
        gz='{sample}/ariba_virfinder/assemblies.fa.gz',
        rep='{sample}/ariba_virfinder/report.tsv'
    params:
        ariba_db = '/media/amr_storage/ariba_ref/virfinder',
        ariba_out = '{sample}/ariba_virfinder'
    shell:
        'ariba run --threads 12 {params.ariba_db} {input.fwd} {input.rev} --force {params.ariba_out}'

rule mlst:
    input:
	'{sample}/{sample}_shovill/contigs.fa'
    output:
        rep='{sample}/ariba_virfinder/report.tsv'
    params:
        ariba_db = '/media/amr_storage/ariba_ref/virfinder',
        ariba_out = '{sample}/ariba_virfinder'
    shell:
        'ariba run --threads 12 {params.ariba_db} {input.fwd} {input.rev} --force {params.ariba_out}'
