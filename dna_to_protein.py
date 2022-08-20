import sys

def translatsioon(rna):
    geenid = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }
    valk = ''
    for i in range(0, len(rna), 3):
            if len(rna[i:i + 3]) == 3:
                koodon = rna[i:i + 3]
                valk+= geenid[koodon]
    return valk

def protein_finder(rna):
    valk = translatsioon(rna)
    start = valk.find('M')
    stop = valk[start::].find('*')
    if len(valk[start:start + stop]) >= 100:
        return 'ORF: ' + valk[start:start + stop]
    return valk[start:start + stop]

dna_fail = open('bakteri_dna.FASTA',encoding='utf-8')
dna = dna_fail.read()
dna = dna.split('\n',maxsplit=1)
dna[1] = dna[1].replace('\n','')
dna_fail.close()
rna_1 = dna[1].replace("T", "U")
rna_2 = rna_1[1::]
rna_3 = rna_2[1::]
rna_4 = rna_1.replace('A', 'u').replace('U', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1]
rna_5 = rna_4[1::]
rna_6 = rna_5[1::]

fail = open('bakteri_valgud.FASTA','w')
fail.writelines('>Järjestus_J1'+"\n")
fail.writelines(protein_finder(rna_1)+"\n")
fail.writelines('>Järjestus_J2'+"\n")
fail.writelines(protein_finder(rna_2)+"\n")
fail.writelines('>Järjestus_J3'+"\n")
fail.writelines(protein_finder(rna_3)+"\n")
fail.writelines('>Järjestus_J4'+"\n")
fail.writelines(protein_finder(rna_4)+"\n")
fail.writelines('>Järjestus_J5'+"\n")
fail.writelines(protein_finder(rna_5)+"\n")
fail.writelines('>Järjestus_J6'+"\n")
fail.writelines(protein_finder(rna_6)+"\n")
fail.close()

sys.stdout.write('>Järjestus_J1'+"\n")
sys.stdout.write(protein_finder(rna_1)+"\n")
sys.stdout.write('>Järjestus_J2'+"\n")
sys.stdout.write(protein_finder(rna_2)+"\n")
sys.stdout.write('>Järjestus_J3'+"\n")
sys.stdout.write(protein_finder(rna_3)+"\n")
sys.stdout.write('>Järjestus_J4'+"\n")
sys.stdout.write(protein_finder(rna_4)+"\n")
sys.stdout.write('>Järjestus_J5'+"\n")
sys.stdout.write(protein_finder(rna_5)+"\n")
sys.stdout.write('>Järjestus_J6'+"\n")
sys.stdout.write(protein_finder(rna_6)+"\n")



