**Motif**

This class will store each individual time a binding motif occurs.

It will store: 
- The motif it corresponds to
- The color assigned to that motif
- Where this motif was found (x coordinate)
- Which gene/sequence it was found on (y coordinate)

**Intron**

This class will store each intron found in the FASTA file.

It will store:
- The start and end position of the intro (x coordinates)
- The specific gene/sequence it corresponds to (y coordinate)
- The width that has been determied for introns (it will be thinner than exons)

**Exon**

The class will store each exon found in the FASTA file.

It will store:
- The start and end position of the intro (x coordinates)
- The specific gene/sequence it corresponds to (y coordinate)
- The width that has been determied for introns (it will be thicker than introns)

**Gene/Sequence Name**

This class will store the name of each FASTA sequence

It wil store:
- The name of the sequence/gene found in the header
- The coordinates (x,y) where the name will be written in the image output 