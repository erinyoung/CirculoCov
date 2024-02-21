import pytest
import subprocess

def run_command(cmd):
    stdout=subprocess.PIPE
    stderr=subprocess.PIPE
    
    proc = subprocess.Popen(cmd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmd}\n{err}")

def test_circulocov():
    """test circulocov"""    
    cmd = "circulocov -n tests/data/test_nanopore.fastq.gz -i tests/data/test_R1.fastq.gz tests/data/test_R2.fastq.gz -g tests/data/test.fasta -o pytest_both -t 1"
    run_command(cmd)

def test_circulocov_all():
    """test circulocov"""    
    cmd = "circulocov -n tests/data/test_nanopore.fastq.gz -i tests/data/test_R1.fastq.gz tests/data/test_R2.fastq.gz -g tests/data/test.fasta -o pytest_all -a -t 1"
    run_command(cmd)

def test_circulocov_nano():
    """test circulocov with nanopore"""    
    cmd = "circulocov -n tests/data/test_nanopore.fastq.gz -g tests/data/test.fasta -o pytest_nano -t 1"
    run_command(cmd)

def test_circulocov_ill():
    """test circulocov"""    
    cmd = "circulocov -i tests/data/test_R1.fastq.gz tests/data/test_R2.fastq.gz -g tests/data/test.fasta -o pytest_illumina -t 1"
    run_command(cmd)

def test_version():
    """test circulocov version"""    
    cmd = "circulocov -v"
    run_command(cmd)

def test_help():
    """test circulocov help"""    
    cmd = "circulocov -h"
    run_command(cmd)