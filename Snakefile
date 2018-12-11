#configfile: 'config.yaml'

#SAMPLE = config['samples']

SAMPLE = [
"11_S2",
"12_S3",
"14_S4",
"15_S5",
"22_S6",
"25_S7",
"30_S4",
"32_S5",
"37_S1",
"38_S2",
"3_S1",
"42_S3",
"50_S10",
"56_S6",
"60_S11",
"64_S7",
"65_S8",
"67_S9",
"68_S10",
"69_S11",
"72_S12",
"76_S13",
"77_S8",
"80_S9"
]

FWD = expand('{sample}/{sample}_unmapped_fwd.fastq', sample=SAMPLE)
REV = expand('{sample}/{sample}_unmapped_rev.fastq', sample=SAMPLE)

rule all:
    input:
        expand('{sample}/{sample}_prokka/{sample}.ffn', sample=SAMPLE),
        expand('{sample}/abricate_amr.txt', sample=SAMPLE),
        expand('{sample}/abricate_plasmid.txt', sample=SAMPLE),
        expand('{sample}/abricate_ecoh.txt', sample=SAMPLE),
        expand('{sample}/abricate_ecolivf.txt', sample=SAMPLE),
        'ecoh_summary.txt',
        'amr_summary.txt',
        'vf_summary.txt',
        'plasmid_summary.txt'

rule trimmomatic:
    input:
        fwd='{sample}_L001_R1.fastq.gz',
        rev='{sample}_L001_R2.fastq.gz'

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
        '{sample}/{sample}_shovill/contigs.fa',
        '{sample}/{sample}_shovill/contigs.gfa'
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
