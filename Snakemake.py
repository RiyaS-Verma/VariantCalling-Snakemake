from glob import glob
from numpy import unique

reads = glob('{}/*'.format(config['readsDir']))#gets unique filesnames
samples = []
for i in reads:
  sampleName = i.replace('{}/'.format(config['readsDir']), '') #format adds the files into {}.
  sampleName = sampleName.replace('{}'.format(config['read1Suffix']), '')
  sampleName = sampleName.replace('{}'.format(config['read2Suffix']), '')
  samples.append(sampleName)
samples = unique(samples)

#get trio numbers
trios = []
for j in samples:
    trioName = j.split('_',)[0]
    trios.append(trioName)
trios = unique(trios)
print(trios)

wildcard_constraints:
  trio='[a-zA-Z0-9]+',
  sample='\w+'

rule all:
  input:
    #expand('ruletoexpand/{sample}.outputsyntax', sample=samples)

rule fastqc:
  input:
    r1 = config['readsDir'] + '/{sample}' + config['read1Suffix'],
    r2 = config['readsDir'] + '/{sample}' + config['read2Suffix']
  output:
    o1 = 'fastqc/{sample}_1_fastqc.html',
    o2 = 'fastqc/{sample}_2_fastqc.html'
  params:
    'fastqc'
  shell:
    'fastqc {input.r1} {input.r2} -o {params}'

rule trimAdapters:
    input:
      r1 = config['readsDir'] + '/{sample}' + config['read1Suffix'],
      r2 = config['readsDir'] + '/{sample}' + config['read2Suffix']
    output:
      'trimmed_reads/{sample}_val_1.fq.gz',
      'trimmed_reads/{sample}_val_2.fq.gz',
    params:
      outDir = 'trimmed_reads',
      suffix = '{sample}',
      minPhred = config['minPhred'],
      minOverlap = config['minOverlap']
    shell:
      'trim_galore --paired --quality {params.minPhred} '
      '--stringency {params.minOverlap} --basename {params.suffix} '
      '--output_dir {params.outDir} {input.r1} {input.r2}'

rule alignReads:
  input:
    r1 = 'trimmed_reads/{sample}_val_1.fq.gz',
    r2 = 'trimmed_reads/{sample}_val_2.fq.gz'
  output:
    sam = temp('aligned_reads/{sample}.sam'),
    bam = 'aligned_reads/{sample}.bam'
  params:
    reference = config['reference']
  threads:
    config['nThreads']
  log:
    'logs/{sample}_bwamem2.log'
  shell:
    'bwa-mem2 mem -t {threads} {params.reference} '
    '{input.r1} {input.r2} > {output.sam} 2> {log}; '
    'samtools sort -@ {threads} -o {output.bam} {output.sam}; '

rule alignReads_bw2: #bowtie2
  input:
    r1 = 'trimmed_reads/{sample}_val_1.fq.gz',
    r2 = 'trimmed_reads/{sample}_val_2.fq.gz'
  output:
    sam = temp('aligned_reads_bw2/{sample}.sam'),
    bam = 'aligned_reads_bw2/{sample}.bam'
  params:
    reference = config['reference_bowtie']
  threads:
    config['nThreads']
  log:
    'logs/{sample}_bowtie2.log'
  shell:
    'bowtie2 -p {threads} -x {params.reference} ' #make sure bowtie indexes is in path 
    '-1 {input.r1} -2 {input.r2} -S {output.sam} 2> {log}; ' 
    'samtools sort -@ {threads} -o {output.bam} {output.sam}; '

#check if reads are sorted

rule MarkDuplicates:
  input:
    bam = 'aligned_reads_bw2/{sample}.bam'
  output:
    md_bam = 'MarkDuplicates/{sample}_md.bam',
    md_txt = 'MarkDuplicates/{sample}_md_metrics.txt'
  shell:
    'picard MarkDuplicates '
    'I={input.bam} O={output.md_bam} M={output.md_txt} '

rule AddRGtoBam:
  input:
    md_bam = 'MarkDuplicates/{sample}_md.bam'
  output:
    rg_bam = 'AddRGtoBam/{sample}_rg.bam'
  params:
    RGLB = config['RGLB'],
    RGPL = config['RGPL'],
    RGPU = config['RGPU'],
  shell:
    'picard AddOrReplaceReadGroups '
    'I={input.md_bam} O={output.rg_bam} '
    'RGID={wildcards.sample} RGLB={params.RGLB} RGPL={params.RGPL} '
    'RGPU={params.RGPU} RGSM={wildcards.sample}; '
    'samtools index {output.rg_bam}'

rule BaseRecalibrate:
  input:
    rg_bam = 'AddRGtoBam/{sample}_rg.bam'
  output:
    rc_table = temp('BaseRecalibrate/{sample}_recal.table'),
    rc_bam = 'BaseRecalibrate/{sample}_recal.bam'
  params:
    reference = config['reference'],
    KnownSites = config['KnownSites']
    JavaOptions = config['JavaOptions']
  shell:
    'gatk --java-options {params.JavaOptions} BaseRecalibrator '
    '-I {input.rg_bam} -R {params.reference} '
    '--known-sites {params.KnownSites} '
    '-O {output.rc_table}; '
    'gatk --java-options {params.JavaOptions} ApplyBQSR '
    '-I {input.rg_bam} -R {params.reference} '
    '--bqsr-recal-file {output.rc_table} '
    '-O {output.rc_bam}'

rule FreeBayes:
  input:
    rc_bam = 'BaseRecalibrate/{sample}_recal.bam'
  output:
    vcf = 'FreeBayes/{sample}.vcf'
  params:
    reference = config['reference']
  log:
    'logs/FreeBayes/{sample}.log'
  threads:
    config['nThreads']
  shell:
    'freebayes -f {params.reference} '
    '{input.rc_bam} >{output.vcf} 2>{log}'

rule HaplotypeCaller:
  input:
    rg_bam = 'BaseRecalibrate/{sample}_recal.bam'
  output:
    vcf = 'HaplotypeCaller/{sample}.g.vcf.gz'
  params:
    reference = config['reference']
  threads:
    config['nThreads']
    config['JavaOptions
  shell:
    'gatk --java-options {params.JavaOptions} HaplotypeCaller'
    '-R {params.reference} '
    '-I {input.recal_bam} '
    '-O {output.vcf} '
    '-ERC GVCF '
    '-G StandardAnnotation '
    '-G AS_StandardAnnotation'

rule Strelka:
  input:
    rc_bam = 'BaseRecalibrate/{sample}_recal.bam'
  output:
    # Strelka results - either use directory or complete file path
    directory('strelka/{sample}')
  log:
    'logs/strelka/germline/{sample}.log'
  params:
    reference = config['reference'],
    workDir = config['workDir']+'/Strelka/{sample}'
  threads:
    config['nThreads']
  conda:
    'envs/strelka.yaml'
  shell:
    'configureStrelkaGermlineWorkflow.py '
    '--bam {input.rc_bam} '
    '--referenceFasta {params.reference} '
    '--exome '
    '--runDir {params.workDir} '
    '&& {params.workDir}/runWorkflow.py '
    '-m local '
    '-j {threads} '
    '2>{log}'

rule GenomicsDBImport:
  input:
    proband = 'HaplotypeCaller/{trio}' + config['probandsuff'] + '.g.vcf.gz',
    mother = 'HaplotypeCaller/{trio}' + config['mothersuff'] + '.g.vcf.gz',
    father = 'HaplotypeCaller/{trio}' + config['fathersuff'] + '.g.vcf.gz'
  output:
    directory('GenomicsDBImport/{trio}')
  params:
    reference = config['reference'],
    workDir = config['workDir'] + '/GenomicsDBImport/{trio}'
  shell:
    'gatk --java-options "-Xmx10g -Xms10g" GenomicsDBImport '
    '-V {input.proband} '
    '-V {input.father} '
    '-V {input.mother} '
    '--tmp-dir /scratch2/rsverma '
    '--genomicsdb-workspace-path {params.workDir} '
    '--consolidate true '
    '-L /path/to/interval/list'

rule GenotypeGVCFs: #genotype based on trios
  input:
    vcf = config['workDir'] +'/GenomicsDBImport/{trio}'
  output:
    vcf = 'GenotypeGVCFs/{trio}.vcf.gz'
  params:
    reference = config['reference']
  shell:
    'gatk --java-options "-Xmx4g" GenotypeGVCFs '
    '-R {params.reference} '
    '-V gendb://{input.vcf} '
    '-O {output.vcf}'

rule ValidateVariants:
  input:
    vcf = 'GenotypeGVCFs/{trio}.vcf.gz'
  params:
    reference = config['reference'],
    knownsites = config['KnownSites']
  log:
    'logs/ValidateVariants/{trio}'
  shell:
    'gatk ValidateVariants '
    '-R {params.reference} '
    '-V {input.vcf} '
    '--dbsnp {params.knownsites} > {log}'

rule VariantRecalibrator:
  input:
    vcf = 'GenotypeGVCFs/{trio}.vcf.gz'
  output:
    recal = temp('VariantRecalibrator/{trio}.AS.recal'),
    tranches = temp('VariantRecalibrator/{trio}.AS.tranches'), #if error, replace trio.AS.tranches to output.AS.tranches
    plot = 'VariantRecalibrator/{trio}.plot.AS.R',
    vqsr = 'VariantRecalibrator/{trio}_VQSR.vcf.gz'
  params:
    reference = config['reference'],
    hapmap = config['hapmap'],
    omni = config['omni'],
    1000G = config['1000G'],
    dbsnp = config['dbsnp']
  shell:
    'gatk VariantRecalibrator '
    '-R {params.reference} '
    '-V {input.vcf} '
    '-AS '
    '--resource:hapmap,known=false,training=true,truth=true,prior=15.0 {params.hapmap} '
    '--resource:omni,known=false,training=true,truth=false,prior=12.0 {params.omni} '
    '--resource:1000G,known=false,training=true,truth=false,prior=10.0 {params.1000G} '
    '--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {params.dbsnp} '
    '-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR '
    '-mode SNP '
    '-O {output.recal} '
    '--tranches-file {output.tranches} '
    '--rscript-file {output.plot}; '
    'gatk ApplyVQSR '
    '-R {params.reference} '
    '-V {input.vcf} '
    '-O {output.vqsr} '
    '-AS '
    '--truth-sensitivity-filter-level 99.0 '
    '--tranches-file {output.tranches} '
    ' --recal-file {output.recal} '
    '-mode SNP'

rule CalculateGenotypePosteriors:
  input:
    vqsr = 'VariantRecalibrator/{trio}_VQSR.vcf.gz'
  output:
    genoposteriors = 'GenotypePosteriors/{trio}_GenoPosteriors.vcf.gz'
  params:
    ped = config['ped']
    JavaOptions = config['JavaOptions']
  shell:
    'gatk --java-options {params.JavaOptions} CalculateGenotypePosterior '
    '-V {input.vqsr} '
    '-O {output.genoposteriors} '
    '-ped {params.ped} '

rule VariantFiltration:
  input:
    genoposteriors = 'GenotypePosteriors/{trio}_GenoPosteriors.vcf.gz'
  output:
    varfilt = 'VariantFiltration/{trio}_VariantFiltered.vcf.gz'
  params:
    reference = config['reference']
  shell:
    'gatk VariantFiltration '
    '-R {params.reference} '
    '-V {input.genoposteriors} '
    '-O {output.varfilt} '
    '--filter-name "LowDP" '
    '--filter-expression "DP < 10" '
    '--filter-name "QUAL300" '
    '--filter-expression "QUAL < 300.0" '
    '--filter-name "FS60" '
    '--filter-expression "FS > 60.0" '
    '--filter-name "MQ40" '
    '--filter-expression "MQ < 40.0" '
    '--filter-name "SOR4" '
    '--filter-expression "SOR > 4.0" '
    '--filter-name "MQRankSum-8" '
    '--filter-expression "MQRankSum < -8.0" '
    '--filter-name "ReadPosRankSum-8" '
    '--filter-expression "ReadPosRankSum < -8.0" '
    '--genotype-filter-name "LowGQ" '
    '--genotype-filter-expression "GQ < 20" '
    '--tmp-dir //scratch2/rsverma'

rule VariantAnnotator:
  input:
    varfilt = 'VariantFiltration/{trio}_VariantFiltered.vcf.gz'
  output:
    varanno = 'VariantAnnotator/{trio}_DeNovo.vcf.gz'
  params:
    reference = config['reference'],
    ped = config['ped'],
    dbsnp = config['dbsnp']
  shell:
    'gatk VariantAnnotator '
    '-V {input.varfilt} '
    '-A PossibleDeNovo '
    '-ped {params.ped} '
    '-O {output.varanno} '
    '-dbsnp {params.dbsnp}'

rule snpEff:
  input:
    varanno = config['workDir'] + 'VariantAnnotator/{trio}_DeNovo.vcf.gz'
  output:
    snpEff = config['workDir'] + 'snpEff/{trio}.ann.vcf'
    html = config['workDir'] + 'snpEff/{trio}.html'
  params:
     snpEffconfig = config['snpEffconfig']
  shell:
    'snpEff -Xmx8g '
    '-c {params.snpEffconfig} '
    '-v -no-intron '
    '-lof GRCh38.99 '
    '-s {output.html} '
    '{input.varanno} > {output.snpEff}'

rule snpSift:
  input:
    snpEff = 'snpEff/{trio}.ann.vcf'
  output:
    snpSift = 'snpSift/{trio}_hiconfDeNovo.vcf'
  shell:
    'cat {input.snpEff}| '
    'SnpSift filter "(exists hiConfDeNovo)" > {output.snpSift}'


rule VariantSelector:
  input:
    snpSift = 'snpSift/{trio}_hiconfDeNovo.vcf'
  output:
    selected = 'VariantSelector/{trio}_HiQual_DeNovo.vcf.gz'
  shell:
    'gatk SelectVariants '
    '-V {input.snpSift} '
    '--exclude-filtered '
    '-O {output.selected}'

rule VariantsToTable:
  input:
    selected = 'VariantSelector/{trio}_HiQual_DeNovo.vcf.gz'
  output:
    table = 'VariantsToTable/{trio}.table'
  shell:
    'gatk VariantsToTable '
    '-V {input.selected} '
    '-F CHROM -F POS -F QUAL -F DP -F HET '
    '-F HOM-REF -F HOM-VAR -F VAR -F NSAMPLES '
    '-F NCALLED -F MULTI-ALLELIC -GF GQ '
    '-O {output.table}'

rule snpSiftTable:
  input:
    selected = 'VariantSelector/{trio}_HiQual_DeNovo.vcf.gz'
  output:
    table = 'snpSiftTable/{trio}_DeNovo_snpSift.txt'
  shell:
    'SnpSift extractFields -s "," -e "." '
    '{input.selected} CHROM POS REF ALT FILTER DP '
    '"ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" '
    '"ANN[*].BIOTYPE" "ANN[*].ERRORS" > {output.table}'
