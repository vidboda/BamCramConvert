To install BamCramConvert and the bam2cram-check submodule either run:

    git clone --recurse-submodules https://github.com/beboche/BamCramConvert.git

or

    git clone https://github.com/beboche/BamCramConvert.git

then cd to BamCramConvert and

    cd BamCramConvert
    git submodule init
    git submodule update 

BamCramConvert requires a decently recent version of samtools
bam2cram-check (option -c) requires python > 3.5 (python path can be passed as an argument with -p) ans samtools > 1.3 (the same, using -st) 
