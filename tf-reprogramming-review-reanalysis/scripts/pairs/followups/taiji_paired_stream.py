"""Check every interleaved spot and feed paired reads to Bowtie2 unmodified."""
from collections import Counter
from pathlib import Path
import sys

source, dest = sys.stdin.buffer, sys.stdout.buffer
reads = 0
size_histogram = Counter()
first, last = None, None
chunk = bytearray()
while header1 := source.readline():
    sequence1, plus1, qual1, header2, sequence2, plus2, qual2 = (
        source.readline() for _ in range(7)
    )
    key1, key2 = header1.strip(), header2.strip()
    if not key1.startswith(b"@") or not key1.endswith(b"/1") or key2 != key1[:-1]+b"2":
        raise ValueError(f"Unpaired FASTQ interleaved record {reads+1}: {header1!r}, {header2!r}")
    if not plus1.startswith(b"+") or not plus2.startswith(b"+"):
        raise ValueError(f"Invalid separator at {header1!r}")
    if not sequence1 or not sequence2 or len(sequence1) != len(qual1) or len(sequence2) != len(qual2):
        raise ValueError(f"Bad sequence or quality length at {header1!r}")
    reads += 1
    last = key1.decode("ascii") if reads == 1 else key1
    if reads == 1:
        first = key1
    size_histogram[(len(sequence1)-1, len(sequence2)-1)] += 1
    chunk.extend(b"".join((header1, sequence1, plus1, qual1, header2, sequence2, plus2, qual2)))
    if len(chunk) >= 1024 * 1024:
        dest.write(chunk)
        chunk.clear()
if not reads:
    raise ValueError("Zero paired FASTQ spots")
dest.write(chunk)
dest.flush()
if isinstance(last, str):
    last = last.encode()
receipt = Path(sys.argv[1])
with receipt.open("w") as out:
    out.write("mate1_length\tmate2_length\tspots\tfirst_qname\tlast_qname\n")
    for (n1, n2), count in sorted(size_histogram.items()):
        out.write(f"{n1}\t{n2}\t{count}\t{first.decode()}\t{last.decode()}\n")
