VCFY
====
It generates a VCF file with simulated random variants based on the given probability model.

Tools
-----
### vcfy

    Usage: vcfy [OPTIONS] REFERENCE
    
      Generate VCF file with simulated variants in specified range [low, high)
      for the given region of the reference genome. In case that the region is
      not specified the first region is used. If no range is provided, it is
      assumed that the variants are scattered throughout the region.
    
      For more information, consult with the README file.
    
    Options:
      -o, --output FILENAME      Write to this file instead of standard output.
      -m, --mutation-rate FLOAT  Base mutation rate.  [required]
      -r, --region TEXT          Region ID (default=first region in the reference)
      -l, --low INTEGER          Range lower bound (default=first locus in the
                                 region)
      -h, --high INTEGER         Range upper bound (default=last locus in the
                                 region)
      --help                     Show this message and exit.

### ksnper

    Usage: ksnper [OPTIONS] [VCF]
    
      Report the number of SNPs in all k-mers. Specify the k and the VCF file,
      it reports number of SNPS occurred in each k-mer.
    
    Options:
      -o, --output FILENAME           Write to this file instead of standard
                                      output.
      -r, --reference FILENAME        Reference genome FASTA file. It will be
                                      inferred from VCF header, if not specified.
      -k INTEGER                      The value of k.  [required]
      -c                              Set if the input VCF is compressed
      -d, --dialect [unix|excel-tab|excel]
                                      Use this CSV dialect.  [default: unix]
      --help                          Show this message and exit.
