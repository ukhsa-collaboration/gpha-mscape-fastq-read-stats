# fastq_read_stats

## What is this?

An argparse tool that generates basic statistics given a fastq.gz as input.

## How do I use this?

After pip installing, you can run it directly from the commandline:

```bash
fastq_read_stats input.fastq.gz
```

Which returns something like:

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
QualityStats(read_id='SRR32793206.1', mean_phred_score=0.05724021314818394, read_length=43, gc_content=0.4883720930232558, compression_ratio=1.1627906976744187, paired_base_complexity=0.6904761904761905)
QualityStats(read_id='SRR32793206.2', mean_phred_score=0.05154328183292301, read_length=119, gc_content=0.44537815126050423, compression_ratio=0.6554621848739496, paired_base_complexity=0.7288135593220338)
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
