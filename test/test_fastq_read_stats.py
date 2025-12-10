import pytest
import fastq_read_stats.fastq_read_stats as fqrs


def test_phred_to_ratio():
    assert fqrs.phred_to_ratio(10) == 0.1
    assert fqrs.phred_to_ratio(20) == 0.01
    assert fqrs.phred_to_ratio(30) == 0.001
    assert fqrs.phred_to_ratio(40) == 0.0001

def test_get_gc_content():
    assert fqrs.get_gc_content("GGGAAA") == 0.5
    assert fqrs.get_gc_content("G") == 1
    assert fqrs.get_gc_content("AAAAAA") == 0

def test_get_compression_ratio():
    """
    Test that the values are sane.
    Can't necessarily guarantee they will always
    be identical
    """
    assert fqrs.get_compression_ratio("A" * 9001) < 1

def test_get_paired_base_complexity():
    assert fqrs.get_paired_base_complexity("AAAA") == 0
    assert fqrs.get_paired_base_complexity("ACGT") == 1
    assert fqrs.get_paired_base_complexity("AAAGG") == 0.25
