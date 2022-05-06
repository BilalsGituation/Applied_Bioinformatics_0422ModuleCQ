nextflow.enable.dsl=2

process count_word {
  publishDir "/tmp/", mode: "copy", overwrite: true
  input:
    path infile
    val word
  output:
    path "${infile}.wordcount", emit: wordcount
  script:
    """
    grep -o -i ${word} ${infile} > grepresult
    cat grepresult | wc -l > ${infile}.wordcount
    """
}

workflow {
  inchannel = channel.fromPath("/home/cq/LinuxClass/nonsense.txt")
  wordcount = count_word(inchannel, "lol")

}
