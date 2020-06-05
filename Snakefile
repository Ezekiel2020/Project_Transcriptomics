samps = "RNA-seq-sample-data"
results = "out"
wild_stems = 'Collibri_standard_protocol-HBR-Collibri-100_ng-2_S1 Collibri_standard_protocol-HBR-Collibri-100_ng-3_S2 Collibri_standard_protocol-UHRR-Collibri-100_ng-2_S3 Collibri_standard_protocol-UHRR-Collibri-100_ng-3_S4 KAPA_mRNA_HyperPrep_-HBR-KAPA-100_ng_total_RNA-2_S5 KAPA_mRNA_HyperPrep_-HBR-KAPA-100_ng_total_RNA-3_S6 KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-2_S7 KAPA_mRNA_HyperPrep_-UHRR-KAPA-100_ng_total_RNA-3_S8'
R = ['1', '2']

def get_stranded(wildcards):
    stem = wildcards.stem
    if 'KAPA' in stem:
        return '-s 2'
    elif 'Collibri' in stem:
        return '-s 1'

rule all:
    input:
        results + "/multiqc_report_raw.html",
        results + "/multiqc_report_raw_data.zip",
        expand(samps + "/{stem}_R1_trimmed.fastq.gz", stem=wild_stems),
        expand(samps + "/{stem}_R2_trimmed.fastq.gz", stem=wild_stems),
        results + "/multiqc_report_trimmed.html",
        results + "/multiqc_report_trimmed_data.zip",
        expand(results + "/BAMs/MAPPED_{stem}_Aligned.sortedByCoord.out.bam", stem=wild_stems),
        expand(results + "/STRD/{stem}_str.log", stem=wild_stems),
        expand(results + "/counts/counts_{stem}.log", stem=wild_stems),
        results + "/counts/bla.log"

rule fastqc:
    input:
        samps + "/{stem}.fastq.gz"
    output:
        results + "/fastqc/{stem}_fastqc.zip",
        results + "/fastqc/{stem}_fastqc.html"
    params:
        reslt = results + "/fastqc/"
    shell:
        "fastqc {input} -o {params.reslt}"

rule multiqc_raw:
    input:
        expand(results + "/fastqc/{stem}_L001_R{R}_001_fastqc.zip", stem = wild_stems, R=R)
    output:
        results + "/multiqc_report_raw.html",
        results + "/multiqc_report_raw_data.zip"
    shell:
        "multiqc {input} -o " + results + " -n multiqc_report_raw -z"

rule multiqc_trimmed:
    input:
        expand(results + "/fastqc/{stem}_R{R}_trimmed_fastqc.zip", stem = wild_stems, R=R)
    output:
        results + "/multiqc_report_trimmed.html",
        results + "/multiqc_report_trimmed_data.zip"
    shell:
        "multiqc {input} -o " + results + " -n multiqc_report_trimmed -z"

rule trim_fastq:
    input:
        samps + "/{stem}_L001_R1_001.fastq.gz",
        samps + "/{stem}_L001_R2_001.fastq.gz"
    output:
        samps + "/{stem}_R1_trimmed.fastq.gz",
        samps + "/{stem}_R2_trimmed.fastq.gz"
    params:
        adapters = "adapters.fa"
    shell:
        "bbduk.sh ref={params.adapters} ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10 " +
        "in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]}"
rule STAR_index:
    output:
        directory(results + "/STAR_genome")
    threads: 2
    params:
        ref = "play_data_ref_annot/chr19_20Mb.fa",
        gtf = 'play_data_ref_annot/chr19_20Mb.gtf'
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {output} " +
        "--genomeFastaFiles {params.ref} --sjdbGTFfile {params.gtf}"

rule STAR:
    input:
        samps + "/{stem}_R1_trimmed.fastq.gz",
        samps + "/{stem}_R2_trimmed.fastq.gz",
        results + "/STAR_genome"
    output:
        results + "/BAMs/MAPPED_{stem}_Aligned.sortedByCoord.out.bam"
    threads: 2
    params:
        pref = results + "/BAMs/MAPPED_{stem}_"
    shell:
        "STAR --runThreadN {threads} --genomeDir {input[2]} " +
        "--readFilesIn {input[0]} {input[1]} --readFilesCommand zcat " +
        "--outSAMtype BAM SortedByCoordinate --outFileNamePrefix {params.pref}"

rule sort:
    input:
        results + "/BAMs/MAPPED_{stem}_Aligned.sortedByCoord.out.bam"
    output:
        results + "/BAMs/MAPPED_{stem}_NameSorted.bam"
    shell:
        "samtools sort -n -o {output} {input}"

rule get_strand:
    input:
        results + "/BAMs/MAPPED_{stem}_NameSorted.bam"
    output:
        results + "/STRD/{stem}_str.log"
    params:
        bed = 'play_data_ref_annot/chr19_20Mb.bed'
    shell:
        "infer_experiment.py -r {params.bed} -i {input} > {output}"

rule featureCounts:
    input:
        results + "/BAMs/MAPPED_{stem}_NameSorted.bam"
    output:
        results + "/counts/counts_{stem}.log"
    params:
        strd = get_stranded,
        gtf = "play_data_ref_annot/chr19_20Mb.bed"
    shell:
        "featureCounts -p {params.strd} -t exon -g gene_id -a {params.gtf} -o {output} {input}"

rule deseq2:
    input:
        expand(results + "/counts/counts_{stem}.log", stem=wild_stems)
    output:
        results + "/counts/bla.log"
    shell:
        "Rscript Analysis.r {input}"
