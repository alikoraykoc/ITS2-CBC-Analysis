import os
import subprocess

def test_cbc_analysis_runs():
    changes = os.path.join("tests", "data", "changes_output.txt")
    xfasta = os.path.join("tests", "data", "test_aligned.xfasta")
    out_file = "test_results.tsv"


    result = subprocess.run(
        [
            "python", "cbc_analysis.py",
            "--changes", changes,
            "--xfasta", xfasta,
            "--seq1", "Tetrix_japonica",
            "--seq2", "tetrix_bolivari_ITS2",
            "--out", out_file
        ],
        capture_output=True,
        text=True
    )


    assert result.returncode == 0, f"Script failed: {result.stderr}"


    assert os.path.exists(out_file), "Output file not created"


    with open(out_file) as f:
        lines = f.readlines()
        assert len(lines) > 1, "Output file seems empty"