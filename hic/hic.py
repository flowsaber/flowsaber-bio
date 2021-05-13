from flowsaber.api import *
from flowsaber.tasks.bio.fetch_refgenie import fetch_refgenie


@shell(publish_dirs=["results/fastqc"])
def fastqc(self, fastq_pair: List[str]):
    """fastqc --threads {self.config.cpu} -o `pwd` {fastq_pair[0]} {fastq_pair[1]}"""
    return "*"


@shell(publish_dirs=["results/multiqc"])
def mutliqc(fastqc_results: List[File]):
    """
    for f in {files}; do ln -s -f $f ./; done
    multiqc ./
    """
    files = " ".join(str(f) for f in fastqc_results)
    return "*"


@shell
def digest(fasta: File, chromsizes: File, enzyme: str):
    """cooler digest --out {fasta.stem}.{enzyme}.bed {chromsizes} {fasta} {enzyme}
    """
    return "*.bed"


@shell
def map_parse_sort(self, assembly: str, index: File, chromsizes: File, fastq_pair: Tuple[File, File],
                   parse_kwargs: str = ""):
    """touch {bam}
       bwa mem -SP5M -t {self.config.cpu} {index} {fq1} {fq2} \
       | tee >(samtools sort --threads {self.config.cpu} -O BAM > {bam}) \
       | pairtools parse --assembly {assembly} -c {chromsizes} {parse_kwargs} \
       | pairtools sort --nproc {self.config.cpu} --tmpdir /tmp -o {pair}
    """
    fq1, fq2 = fastq_pair
    name = fq1.name
    bam = f"{name}.bam"
    pair = f"{name}.pair.gz"

    return bam, pair


@shell(publish_dirs=["results/pairs/experiments"])
def filter_pair(self, pair: File, frag: File = None, select: str = "True"):
    """zcat -f {pair} \
        | pairtools select --output-rest {unknown_rest_pair} \"(chrom1!='!') and (chrom2!='!')\" \
        {restrict_cmd} \
        | pairtools select --output-rest {select_rest_pair} -o {filtered_pair} \"{select}\"

        pairtools merge -o {rest_pair} {unknown_rest_pair} {select_rest_pair}
    """
    unknown_rest_pair = f"{pair.stem}.unknown.res.pair.gz"
    select_rest_pair = f"{pair.stem}.select.rest.pair.gz"
    rest_pair = f"{pair.stem}.rest.pair.gz"
    filtered_pair = f"{pair.stem}.filtered.pair.gz"

    restrict_cmd = f" | pairtools restrict -f {frag} " if frag else ""

    return filtered_pair, rest_pair


@shell(publish_dirs=["results/pairs/merged"])
def merge_dedup_pair(self, group_pairs):
    """pairtools merge --nproc {self.config.cpu} {pairs} \
        | pairtools dedup --max-mismatch 2 --mark-dups --output {merged_pair} --output-dups {dup_pair}
    """
    group, pairs = group_pairs
    pairs = ' '.join(str(pair) for pair in pairs)
    merged_pair = f"{group}.dedup.pair.gz"
    dup_pair = f"{group}.dups.pair.gz"
    return merged_pair, dup_pair


@shell(publish_dirs=["results/pairs/stats"])
def stat_pair(pair: File):
    """pairtools stats --output {pair.stem}.stat {pair}"""
    return "*.stat"


@shell(publish_dirs=["results/cools/minres"])
def pair_to_cool(self, pair: File, chromsizes: File, min_resolution: int = 5000, assembly: str = None):
    """cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 \
        {options} {chromsizes}:{min_resolution} {pair} {cool}
    """
    cool = f"{pair.stem}_{min_resolution / 1000}_k.cool"
    options = f" --assembly {assembly}" if assembly else ""
    return cool


@shell(publish_dirs=["results/cools"])
def zoomify_cool(self, cool: File, resolutions="5000N"):
    """cooler zoomify --nproc {self.config.cpu} --balance --resolutions {resolutions} --out {cool.stem}.mcool {cool}"""
    return f"*.mcool"


@shell
def link_to_group(sample: Tuple[str, str, str]):  # symbol link will not change name
    """
    ln {fq1} {grp}.{fq1.name}
    ln {fq2} {grp}.{fq2.name}
    """
    grp = sample[0]
    fq1, fq2 = File(sample[1]), File(sample[2])
    return "*"


@flow
def download_genome(self, assembly):
    # download data
    index = fetch_refgenie(genome=assembly.first(), assets="bwa_index") | view | flatten \
            | filter_(lambda x: str(x).endswith('.fa')) | constant
    fasta_files = fetch_refgenie(genome=assembly.first(), assets="fasta") | flatten
    fasta = fasta_files.filter(by=lambda x: str(x).endswith('.fa')) | first | view | constant
    chromsizes = fasta_files.filter(by=lambda x: str(x).endswith('.sizes')) | first | view | constant

    return fasta, chromsizes, index


@task
def check_inputs(context):
    params = ['assembly', 'enzyme', 'parse_kwargs', 'select', 'samples']
    for param in params:
        assert context.get(param), f"The param: {param} should be valid in input context."
    return context


@flow
def hic_flow():
    """

    Parameters
    ----------
    samples
    assembly
    enzyme
    select
    parse_kwargs

    """

    def group_pair(pair: File):
        return pair.name[:pair.name.index('.')]

    # get user config
    context = get_context() | view | check_inputs
    assembly = context['assembly'] | view | constant
    enzyme = context['enzyme'] | view | constant
    parse_kwargs = context['parse_kwargs'] | view | constant
    select = context['select'] | view | constant
    # download files
    fasta, chromsizes, index = download_genome(assembly=assembly)

    # main step
    fastq_pairs = context['samples'] | flatten | view | link_to_group | view
    fragments = digest(fasta.first(), chromsizes, enzyme) | view | constant
    bams, exp_pairs = map_parse_sort(assembly, index, chromsizes, fastq_pairs, parse_kwargs) | view | split(2)
    filtered_pairs, rest_pairs = filter_pair(exp_pairs, fragments, select) | split(2)
    merged_pairs, dup_pairs = merge_dedup_pair(filtered_pairs.group(by=group_pair)) | split(2)
    mcools = pair_to_cool(merged_pairs, chromsizes, assembly=assembly) | zoomify_cool
    # stats
    multiqc_results = fastq_pairs | fastqc | flatten | collect | mutliqc
    pair_stats = [exp_pairs, merged_pairs] | mix | stat_pair

    return bams, exp_pairs, merged_pairs, mcools, multiqc_results, pair_stats


context = {
    'task_config': {
        # 'cache': None,
        'cpu': 7,
    },
    "flow_workdir": '/store/qzhong/hictools_test',
    "samples": [
        ('group1', 'SRR1658777_1.fastq', 'SRR1658777_2.fastq'),
        ('group1', 'SRR1658777_copy_1.fastq', 'SRR1658777_copy_2.fastq'),
        ('group2', 'SRR1658777_11.fastq', "SRR1658777_22.fastq"),
    ],
    'assembly': "hg19",
    "enzyme": "NcoI",
    "select": "((pair_type=='UU') or (pair_type=='UR') or (pair_type=='RU'))",
    "parse_kwargs": "--min-mapq 0 --drop-sam --drop-seq --drop-readid",
}

f = hic_flow()
if __name__ == "__main__":
    run(f, context=context)
