rule consensus_merge:  # TODO: ta mal
    input: lambda w: build_targets_filtering(w, f"output/variants/consensus/{w.study}/{{}}/{{}}/{{}}/{{}}_{{}}_{{}}/sample.fasta", ("sample_accession", "instrument_platform", "run_accession", "library_layout", "fastq_ftp", "library_strategy"), study_accession=w.study)
    output: temp("output/pangolin/consensus_merge/{study}/sequences.fasta")
    resources:
        runtime = "15m",
        mem_gb = 2
    shell: "cat {input:q} > {output:q}"


rule pangolin_assignment:
    threads: 8
    conda: "../envs/lineages.yaml"
    shadow: "minimal"
    input:
        fasta = "output/pangolin/consensus_merge/{study}/sequences.fasta"
    output:
        table = temp("output/pangolin/pangolin_assignment/{study}/pangolin.csv")
    resources:
        mem_gb = 8
    log: "output/logs/pangolin/pangolin_assignment/{study}.txt"
    shell: "pangolin {input.fasta:q} --outfile {output.table:q} --threads {threads} >{log:q} 2>&1"


rule pangolin_assignment_merge:
    input: lambda w: build_targets(w, "output/pangolin/pangolin_assignment/{}/pangolin.csv", ("study_accession",))
    output: "output/pangolin/pangolin.csv"
    resources:
        runtime = "15m",
        mem_gb = 2
    shell: "head -1 {input[0]:q} >{output:q} && tail +2 -q {input:q} >>{output:q}"
