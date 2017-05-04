from src.coevolution import coevol

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


from src.network import write_network_xml
from Bio.Seq import translate

output = coevol(['P1.nex'], itype = "cds", rCut=0.8, method = "Intra", pdbfile="5k6k.pdb", chain1="A", atom_interval1="215-765")

write_network_xml(output)

with open("network_out.txt", "w") as fp:
    for key, val in output.items():
        keyObj = key.split("-")
        fp.write("%s\t%s\n" %(keyObj[0], keyObj[1]))

