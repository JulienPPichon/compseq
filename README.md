Analyse the conservation of reversed complementary motifs in an alignment file.


usage: compseq.py [-h] [-i FILENAME] [-o OUTPUT] [-d MINIMAL_DISTANCE]
                  [-k KMER_LENGTH] [-f] [-s [SELECT_SEQUENCE ...]]
                  [-r [REJECT_SEQUENCE ...]]


optional arguments:
-h, --help            show this help message and exit

-i FILENAME, --filename FILENAME : Take nexus alignement file as input.

-o OUTPUT, --output OUTPUT
Name of the tsv output file containing the
individual's name, the kmer sequence, the start
position of the kmer sequence, the start position of
its reverse complementary sequence, number of
individuals possessing this kmer sequence and its
reversed complementary at these positions, number of
individuals possessing a kmer sequence and its
reverted complementary at these positions, kmer
length.

-d MINIMAL_DISTANCE, --minimal_distance MINIMAL_DISTANCE
Minimal distance between a kmer and its reverse
complementary sequence.

-k KMER_LENGTH, --kmer_length KMER_LENGTH
Complementary sequences length.

-f, --filter_kmer     Filter out complementary sequences belonging to a
greater kmer length. (e. g: AGGGA TCCCT. AGGG or GGGA
will not count because AGGGA have a reversed
complementary sequence).

-s [SELECT_SEQUENCE ...], --select_sequence [SELECT_SEQUENCE ...]
Only the specified sequences will count for the conservation.

-r [REJECT_SEQUENCE ...], --reject_sequence [REJECT_SEQUENCE ...]
The specified sequences will not count for the conservation.
