nextflow.enable.dsl = 2

//  container "https://depot.galaxyproject.org/singularity/"

// storeDir skips the process if the file exists at the destination
// publishDir runs the process and overwrites the file at the destination
// so you can manipulate these commands like the IF EXISTS operator in SQL

// container: if mentioned AND profile -singulatiry is given on the command line
// AND we have nextflow.config with singularity instructions inside,
// it is called from the link.
// if singularity image is not in the singularity cache directory mentioned
// in nextflow.config, it is downloaded to it

process prefetch {
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"
  input:
    val sra_id
  output:
    path "${sra_id}/${sra_id}.sra" // this folder structure is the special one that is produced by prefetch
  script:
    """
    prefetch ${sra_id} # produces the SRR*/SRR*.sra
    """
}

process fasterq_dump {
  container "https://depot.galaxyproject.org/singularity/sra-tools:2.11.0--pl5262h314213e_0"
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  input:
    path infile_sra // this file path will be linked to a working dir of this process
  output:
    path "${infile_sra.getSimpleName()}*.fastq"
  script:
    """
    fasterq-dump "${infile_sra}" # translates to fastq-dump SRR*.sra
    """
}

process fastqc {
  container "https://depot.galaxyproject.org/singularity/fastqc%3A0.11.9--hdfd78af_1"
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  input:
    path infile_fastq // this file path will be linked to a working dir of this process
  output:
    path "${params.accession}*"
    path "${params.accession}*_fastqc*.zip", emit: zipped
  script:
    """
    fastqc ${params.accession}*.fastq # runs the file(s) associated with a SRR number
    """
}

process fFox {
//  container "https://depot.galaxyproject.org/singularity/fastqc%3A0.11.9--hdfd78af_1"
//  publishDir "${params.outdir}", mode: "copy", overwrite: true
  input:
    path infile_fastqc // this file path will be linked to a working dir of this process
//  output:
//    path "*.html"
  script:
    """
    firefox ${params.accession}*.html # opens the html report of your FastQC
    """
}

process fastP {
  container "https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0"
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  input:
    path infile_fastqc
    val accession
  output:
    path "${params.accession}*_fastp*.fastq", emit: fastq_cut
    path "${params.accession}*_fastp*.html", emit: fastpreport
    path "${params.accession}*"
  script:
    if (infile_fastqc instanceof List) {
    """
    fastp -i ${params.accession}_1.fastq -o ${params.accession}_fastp_1.fastq -I ${params.accession}_2.fastq -O ${params.accession}_fastp_2.fastq -h ${params.accession}_fastp.html -j ${params.accession}_fastp.json
    """
    } else {
      """
      fastp -i ${params.accession}.fastq -o ${params.accession}_fastp.fastq -h ${params.accession}_fastp.html -j ${params.accession}_fastp.json
      """
      }}

process multiqc {
  container "https://depot.galaxyproject.org/singularity/multiqc:1.9--py_1"
  publishDir "${params.outdir}", mode: "copy", overwrite: true
  input:
    path report_files
  //output:
  //  path "multi*"
  script:
    """
    multiqc .
    """
}

workflow {
  // start process prefetch with --accession channel assigned in cmd line (--accession)
  // throw the prefetch output path into a channel called sraresult
  // pipe the sraresult into the process fasterq_dump
  sraresult = prefetch(params.accession)
  fastqs = fasterq_dump(sraresult)
  fastqcs = fastqc(fastqs.flatten())
  fastpout = fastP(fastqs, params.accession)
  //fastpout.fastpreport.view()
  //fastqcs.zipped.view()
  multiqc_input = fastpout.fastpreport.concat(fastqcs.zipped)
  multiqc_input.view()
  multiqc(multiqc_input.collect())
  //newprocess(fastpout.fastq_cut)
  //ffox = fFox(fastqcs)
}
