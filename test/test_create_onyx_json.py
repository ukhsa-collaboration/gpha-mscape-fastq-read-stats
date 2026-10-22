import json
from pathlib import Path

import pandas as pd
import pytest

import create_qc_json as qc


def make_qc_tsv(path: Path, n_rows: int = 2):
    df = pd.DataFrame(
        {
            "mean_phred_score": [30 + i for i in range(n_rows)],
            "read_length": [100 + i * 10 for i in range(n_rows)],
            "gc_content": [0.5 + (i * 0.05) for i in range(n_rows)],
            "compression_ratio": [1.0 + (i * 0.1) for i in range(n_rows)],
            "paired_base_complexity": [0.2 + (i * 0.05) for i in range(n_rows)],
        }
    )
    df.to_csv(path, sep="\t", index=False)
    return path


def test_get_unclassified_results_single_end():
    tmp = Path("tmp_test_single.tsv")
    make_qc_tsv(tmp, n_rows=2)

    try:
        result = qc.get_unclassified_results([tmp])

        assert "read1" in result
        assert "read2" not in result
        assert result["read1"]["unclassified_read_count"] == 2
        assert "mean_phred_score" in result["read1"]
        assert "read_length" in result["read1"]
    finally:
        tmp.unlink(missing_ok=True)


def test_get_unclassified_results_paired_end():
    r1 = Path("tmp_test_r1.tsv")
    r2 = Path("tmp_test_r2.tsv")
    make_qc_tsv(r1, n_rows=2)
    make_qc_tsv(r2, n_rows=3)

    try:
        result = qc.get_unclassified_results([r1, r2])

        assert "read1" in result
        assert "read2" in result
        assert "read1" in result
        assert result["read1"]["unclassified_read_count"] == 2
        assert result["read2"]["unclassified_read_count"] == 3
        assert "unclassified_read_count" not in result
    finally:
        r1.unlink(missing_ok=True)
        r2.unlink(missing_ok=True)


def test_get_unclassified_results_raises_for_missing_file():
    missing = Path("this_file_does_not_exist.tsv")

    with pytest.raises(FileNotFoundError):
        qc.get_unclassified_results([missing])


def test_get_unclassified_results_raises_for_more_than_two_files():
    files = [Path("a.tsv"), Path("b.tsv"), Path("c.tsv")]

    for f in files:
        make_qc_tsv(f, n_rows=1)

    try:
        with pytest.raises(ValueError):
            qc.get_unclassified_results(files)
    finally:
        for f in files:
            f.unlink(missing_ok=True)


def test_make_pipeline_dict():
    pipeline = "fastqc,1.2.3,https://example.com"
    result = qc.make_pipeline_dict(pipeline)

    assert result == {
        "name": "fastqc",
        "version": "1.2.3",
        "homePage": "https://example.com",
    }


def test_create_analysis_fields_sets_top_result_from_read_counts(monkeypatch):
    class FakeOnyxAnalysis:
        def __init__(self):
            self.pipeline_name = None
            self.pipeline_version = None
            self.pipeline_url = None
            self.outputs = None

        def add_analysis_details(self, **kwargs):
            pass

        def add_methods(self, methods_dict):
            return False

        def add_versions_to_methods(self, tool_versions, onyx_versions):
            return False

        def add_results(self, top_result, results_dict):
            assert top_result == 5
            return False

        def add_server_records(self, **kwargs):
            pass

        def check_analysis_object(self, publish_analysis=False):
            return (False, False)

        def write_analysis_to_json(self, result_file):
            return result_file

    monkeypatch.setattr(qc.oa, "OnyxAnalysis", FakeOnyxAnalysis)
    monkeypatch.setattr(
        qc.oa,
        "get_data_and_versions_from_onyx",
        lambda sample_id, server, fields: (None, {}, 0),
    )

    result_dict = {
        "read1": {"unclassified_read_count": 2},
        "read2": {"unclassified_read_count": 3},
    }

    analysis, exitcode = qc.create_analysis_fields(
        record_id="SAMPLE1",
        unclassified_qc_results=result_dict,
        server="devscape",
        pipeline_info={"name": "pipe", "version": "1.0", "homePage": "https://example.com"},
        orange_box_version="v1",
    )

    assert isinstance(analysis, FakeOnyxAnalysis)
    assert exitcode == 0


def test_main_store_onyx_writes_json(monkeypatch, tmp_path):
    # Arrange: fake args
    args = {
        "climbid": "SAMPLE_TEST",
        "server": "devscape",
        "output": str(tmp_path),
        "pipeline_info": "pipe,1.0,https://example.com",
        "readQC_result": "True",
        "readQC_file": [],
        "orange_box_version": "v1",
        "store_onyx": True,
        "test_onyx": False,
        "prod_onyx": False,
    }

    # Create valid readQC file
    qcf = tmp_path / "qc.tsv"
    make_qc_tsv(qcf, n_rows=2)

    args["readQC_file"] = [str(qcf)]

    monkeypatch.setattr(qc, "get_args", lambda: type("Args", (), args)())
    monkeypatch.setattr(
        qc,
        "set_up_logger",
        lambda *args, **kwargs: None,
    )

    class FakeOnyxAnalysis:
        def __init__(self):
            self.pipeline_name = None
            self.pipeline_version = None
            self.pipeline_url = None
            self.outputs = None

        def add_analysis_details(self, **kwargs):
            pass

        def add_methods(self, methods_dict):
            return False

        def add_versions_to_methods(self, tool_versions, onyx_versions):
            return False

        def add_results(self, top_result, results_dict):
            return False

        def add_server_records(self, **kwargs):
            pass

        def check_analysis_object(self, publish_analysis=False):
            return (False, False)

        def write_analysis_to_json(self, result_file):
            payload = {"status": "ok"}
            with open(result_file, "w") as fh:
                json.dump(payload, fh)
            return str(result_file)

    monkeypatch.setattr(qc.oa, "OnyxAnalysis", FakeOnyxAnalysis)
    monkeypatch.setattr(
        qc.oa,
        "get_data_and_versions_from_onyx",
        lambda sample_id, server, fields: (None, {}, 0),
    )

    rc = qc.main()

    assert rc == 0
    out_json = tmp_path / "SAMPLE_TEST.unclassified_readQC.analysis_fields.json"
    assert out_json.exists()