from ast import Is
from multiprocessing.spawn import prepare
from pickle import TRUE


def ConvertTomRNA(dna):
    mRNA = ""
    dna = dna.upper()
    for i in dna:
        if i == "A":
            mRNA += "U"
        elif i == "T":
            mRNA += "A"
        elif i == "G":
            mRNA += "C"
        elif i == "C":
            mRNA += "G"
        else:
            mRNA += "O"
    return mRNA


def IsCorrect(kodon, mRNATest, protein):
    if kodon == mRNATest:
        return protein + "-"
    return ""


def GetProtein(kodon):
    amino = ""
    amino += IsCorrect(kodon, "AUG", "Met")
    amino += IsCorrect(kodon, "AUA", "lle")
    amino += IsCorrect(kodon, "AUC", "lle")
    amino += IsCorrect(kodon, "AUU", "lle")

    amino += IsCorrect(kodon, "ACG", "Thr")
    amino += IsCorrect(kodon, "ACA", "Thr")
    amino += IsCorrect(kodon, "ACC", "Thr")
    amino += IsCorrect(kodon, "ACU", "Thr")

    amino += IsCorrect(kodon, "AAG", "Lys")
    amino += IsCorrect(kodon, "AAA", "Lys")
    amino += IsCorrect(kodon, "AAC", "Asn")
    amino += IsCorrect(kodon, "AAU", "Asn")

    amino += IsCorrect(kodon, "AGG", "Arg")
    amino += IsCorrect(kodon, "AGA", "Arg")
    amino += IsCorrect(kodon, "AGC", "Ser")
    amino += IsCorrect(kodon, "AGU", "Ser")

    amino += IsCorrect(kodon, "GUU", "Val")
    amino += IsCorrect(kodon, "GUC", "Val")
    amino += IsCorrect(kodon, "GUA", "Val")
    amino += IsCorrect(kodon, "GUG", "Val")

    amino += IsCorrect(kodon, "GCU", "Ala")
    amino += IsCorrect(kodon, "GCC", "Ala")
    amino += IsCorrect(kodon, "GCA", "Ala")
    amino += IsCorrect(kodon, "GCG", "Ala")

    amino += IsCorrect(kodon, "GAU", "Asp")
    amino += IsCorrect(kodon, "GAC", "Asp")
    amino += IsCorrect(kodon, "GAA", "Glu")
    amino += IsCorrect(kodon, "GAG", "Glu")

    amino += IsCorrect(kodon, "GGU", "Gly")
    amino += IsCorrect(kodon, "GGC", "Gly")
    amino += IsCorrect(kodon, "GGA", "Gly")
    amino += IsCorrect(kodon, "GGG", "Gly")

    amino += IsCorrect(kodon, "UUU", "Phe")
    amino += IsCorrect(kodon, "UUC", "Phe")
    amino += IsCorrect(kodon, "UUA", "Leu")
    amino += IsCorrect(kodon, "UUG", "Leu")

    amino += IsCorrect(kodon, "UCU", "Ser")
    amino += IsCorrect(kodon, "UCC", "Ser")
    amino += IsCorrect(kodon, "UCA", "Ser")
    amino += IsCorrect(kodon, "UCG", "Ser")

    amino += IsCorrect(kodon, "UAU", "Tyr")
    amino += IsCorrect(kodon, "UAC", "Tyr")
    amino += IsCorrect(kodon, "UAA", "sto")
    amino += IsCorrect(kodon, "UAG", "sto")

    amino += IsCorrect(kodon, "UGU", "Cys")
    amino += IsCorrect(kodon, "UGC", "Cys")
    amino += IsCorrect(kodon, "UGA", "sto")
    amino += IsCorrect(kodon, "UGG", "Trp")

    amino += IsCorrect(kodon, "CUU", "Leu")
    amino += IsCorrect(kodon, "CUC", "Leu")
    amino += IsCorrect(kodon, "CUA", "Leu")
    amino += IsCorrect(kodon, "CUG", "Leu")

    amino += IsCorrect(kodon, "CCU", "Pro")
    amino += IsCorrect(kodon, "CCC", "Pro")
    amino += IsCorrect(kodon, "CCA", "Pro")
    amino += IsCorrect(kodon, "CCG", "Pro")

    amino += IsCorrect(kodon, "CAU", "His")
    amino += IsCorrect(kodon, "CAC", "His")
    amino += IsCorrect(kodon, "CAA", "Gln")
    amino += IsCorrect(kodon, "CAG", "Gln")

    amino += IsCorrect(kodon, "CGU", "Arg")
    amino += IsCorrect(kodon, "CGC", "Arg")
    amino += IsCorrect(kodon, "CGA", "Arg")
    amino += IsCorrect(kodon, "CGG", "Arg")

    return amino


def ConvertToProtein(mRNA):
    if len(mRNA) % 3 != 0:
        print("ERROR mRNA består inte bara av kodon %3")
        return

    protein = ""
    for i in range(0, len(mRNA), 3):
        kodon = mRNA[i] + mRNA[i + 1] + mRNA[i + 2]
        protein += GetProtein(kodon)
    return protein


def FixProtein(prePro):
    protein = ""
    foundStart = False
    for i in range(0, len(prePro), 4):
        amino = prePro[i] + prePro[i + 1] + prePro[i + 2] + prePro[i + 3]
        if amino == "Met-":
            foundStart = TRUE
            continue

        if foundStart:
            if amino == "sto-":
                return protein[:-1]
            else:
                protein += amino
    return protein[:-1]


def main():
    dna = input("DNA: ")
    mRNA = ConvertTomRNA(dna)
    protein = ConvertToProtein(mRNA)
    protein = FixProtein(protein)

    print("\n\n\n\n"+protein+"\n\n")


main()
