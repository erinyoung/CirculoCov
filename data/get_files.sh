# get nanopore reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR233/016/SRR23387316/SRR23387316_1.fastq.gz

# get illumina reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR230/067/SRR23080367/SRR23080367_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR230/067/SRR23080367/SRR23080367_2.fastq.gz

# get assembly
wget https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_028622915.1 -O GCA_028622915.1.fasta

# get assembly that doesn't map
wget https://www.ebi.ac.uk/ena/browser/api/fasta/GCA_021664795.2?download=true -O GCA_021664795.2.fasta 

# get minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 | tar -jxvf -

wget https://github.com/arq5x/bedtools2/releases/download/v2.31.0/bedtools.static -O bedtools

chmod +x bedtools

export PATH=$(pwd)/minimap2-2.26_x64-linux:$(pwd):$PATH