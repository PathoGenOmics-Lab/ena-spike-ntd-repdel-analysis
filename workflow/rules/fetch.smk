rule efetch_fasta_reference:
    params:
        accession = "NC_045512.2"
    output: "output/reference/sequence.fasta"
    resources:
        runtime = "10m"
    shell: 'curl -s -o {output:q} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={params.accession}&rettype=fasta"'


rule efetch_gff3_reference:
    params:
        accession = "NC_045512.2"
    output: "output/reference/features.gff3"
    resources:
        runtime = "10m"
    shell: 'curl -s -o {output:q} "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={params.accession}&rettype=gff3"'


rule filter_gff3_reference:
    input: "output/reference/features.gff3"
    params:
        selection = {"type": ["gene"]}
    output: "output/reference/features.filtered.gff3"
    resources:
        runtime = "10m"
    run:
        import pandas as pd
        # see: https://gmod.org/wiki/GFF3
        COLUMNS = ("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "tags")
        gff3 = pd.read_csv(input[0], sep="\t", names=COLUMNS, comment="#")
        for column, values in params.selection.items():
            gff3 = gff3[gff3[column].isin(values)]
        gff3.to_csv(output[0], sep="\t", index=False)
