The sequencing alignment pipeline was followed by:
[FASTQ]
  │
  └── alignment.sh → [SAM] + [Unaligned FASTQ] + [Mapping log]
          │
          └── bamcount.sh → [q60 BAM] → [pairend BAM] → [mt_info vs no_mt_info BAM]
                                          │
                                          ├─→ [spliced/unspliced BAM] → bamcount.sh → [count tables]
                                          └─→ juncbed.sh → [junction bed + annotations]
