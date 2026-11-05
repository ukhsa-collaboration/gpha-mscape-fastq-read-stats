# fastq_read_stats

## What is this?

An argparse tool that generates basic statistics given a fastq.gz as input and can create a json analysis table for ingest into onyx via create_onyx_json.py. 

## How do I use this?

### 1. fastq_read_stats
After pip installing, you can run it directly from the commandline:

```bash
fastq_read_stats input.fastq.gz > input.fastq.stats
```

Which returns something like in the output table:

```
input_filename  read_id mean_phred_score        read_length     gc_content      compression_ratio       paired_base_complexity
input.fastq.gz   SRR32793206.1   0.05724021314818394     43      0.4883720930232558      1.1627906976744187      0.6904761904761905
input.fastq.gz   SRR32793206.2   0.05154328183292301     119     0.44537815126050423     0.6554621848739496      0.728813559322033
```

Alternatively, you can also call it directly from a script:

```python
import fastq_read_stats.fastq_read_stats as fqrs
import gzip
from Bio import SeqIO

## read a gzipped file
with gzip.open("input.fastq.gz", "rt") as infile:
    ## load fastq into records
    recs = SeqIO.parse(infile, "fastq")

    ## loop through the records
    for record in recs:
        ## generate the stats
        stats = fqrs.get_quality_stats(record)
        ## print them so we can see
        print(stats)
```

Which returns something like:

```python
QualityStats(
    read_id="SRR32793206.1",
    mean_phred_score=0.05724021314818394,
    read_length=43,
    gc_content=0.4883720930232558,
    compression_ratio=1.1627906976744187,
    paired_base_complexity=0.6904761904761905,
)
QualityStats(
    read_id="SRR32793206.2",
    mean_phred_score=0.05154328183292301,
    read_length=119,
    gc_content=0.44537815126050423,
    compression_ratio=0.6554621848739496,
    paired_base_complexity=0.7288135593220338,
)
```
### 2. create_onyx_json
A simple script to generate onyx analysis fields from unclassified read QC results. It reads in the table and generates some summary statistics for the unclassified reads for storage in the analysis table. 

```bash
create_onyx_json -i C-FAKEID -s mscape -o output_dir -p 'pipeline name, version, homepage' -r True -c orange_box_version -f C-FAKEID.stats
```

Example Output
```
identifiers
name"gpha-mscape-fastq-read-stats-nf"
description"This is an analysis to perform read QC for unclassified reads"
analysis_date"2026-10-19"
pipeline_name"gpha-mscape-fastq-read-stats-nf"
pipeline_version" v1.5.0"
pipeline_url" https://github.com/ukhsa-collaboration/gpha-mscape-fastq-read-stats-nf"
outputs"No s3 outputs"
methods"{"unclassified_version": "1.0", "onyx_versions_hash": "53e3c51ec1c1c40047d8329f720f7e18d9f093367b955603b713c9dc2919e66c", "versions": [{"name": "classifier_version", "version": "2.1.2"}, {"name": "classifier_db_date", "version": "2026-02-26"}, {"name": "ncbi_taxonomy_date", "version": "2026-03-01"}, {"name": "scylla_version", "version": "v2.3.0"}, {"name": "sylph_db_version", "version": "gtdb-r232-c200-dbv1"}, {"name": "alignment_db_version", "version": "0.4.0"}, {"name": "orange_box_version", "version": "True"}]}"
result 4855
result_metrics"{"mean_phred_score": {"min": 0.012, "max": 0.12, "mean": 0.06, "std": 0.028, "median": 0.052}, "read_length": {"min": 253.0, "max": 5513.0, "mean": 2573.0469618949537, "std": 699.6761743315348, "median": 2508.0}, "gc_content": {"min": 0.273, "max": 0.642, "mean": 0.447, "std": 0.093, "median": 0.446}, "compression_ratio": {"min": 0.263, "max": 0.490, "mean": 0.334, "std": 0.012, "median": 0.33}, "paired_base_complexity": {"min": 0.611, "max": 0.795, "mean": 0.731, "std": 0.029, "median": 0.730}, "unclassified_read_count": 4855}"

```
## What do these stats mean?

| Field                  | Definition                                                                                                                                                                                         | Format |
| ---------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------ |
| read_id                | ID from the read header in the fastq                                                                                                                                                               | str    |
| mean_phred_score       | Average of the per-base [Phred scores](https://en.wikipedia.org/wiki/Phred_quality_score) for a read, representing the average probability that bases are incorrectly called, expressed as a ratio | float  |
| read_length            | Length of the read in bases                                                                                                                                                                        | int    |
| gc_content             | Ratio of count of GC bases to total length of read                                                                                                                                                 | float  |
| compression_ratio      | Ratio of length in bytes of a gzip compressed verison of a read to the uncompressed version. Proxy for information content                                                                         | float  |
| paired_base_complexity | Ratio of bases that differ from the previous base to total length. Same formula as [fastp](https://github.com/OpenGene/fastp)                                                                      | float  |

## What can I do with these metrics?

You can attempt to interpret them...

A high compression ratio indicates that the read is highly repetitive, and the information content is low.

A low compression ratio indicates that the read is nearly random - random noise is relatively incompressible - though you may need to factor in the read length as a very short read will likely not compress very much.

An average Phred score of 0.2 would indicate that for any given base, there is a 20% chance that it is incorrectly called. Expressed a different way, the expected identity of the read to the "true sequence" is around 80%.

GC content is complex to interpret, but coding regions generally have a higher GC content than the rest of the genome. It may be informative to plot the GC content for all reads in a sample and examine the outliers.

Paired base complexity is a somewhat unintuitive metric, with annoying edge cases - for example, "ACACACACAC" is clearly repetitive, but as each successive base is different from the previous it still gets a high complexity score. Take it with a pinch of salt.
