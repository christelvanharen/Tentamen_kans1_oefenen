import re

def main():
    fasta = "TAIR10_pep_20101214.fa"
    headers, seqs = lees_fasta(fasta)
    print(headers)
    print(seqs)


def lees_fasta(bestand):
    try:
        bestand = open(bestand)
        headers = []
        seqs = []
        seq = ""
        for regel in bestand:
            regel = regel.strip()
            if ">" in regel:
                if seq != "":
                    seqs.append(seq)
                    seq = ""
                headers.append(regel)
            else:
                seq += regel.strip()
        seqs.append(seq)
        return headers, seqs
    except FileNotFoundError:
        print("Kan bestand niet vinden, probeer opnieuw")
    except IOError:
        print("Het bestand is niet in de juiste format, probeer het anders.")
    except TypeError:
        print("Er is een onjuist bestand ingegeven")

def referentie(headers, seqs):
    try:
        for seq in seqs:
        overeenkomst = re.search("[LIVMFYC].[HY].D[LIVMFY]K..N["
                                 "LIVMFYCT]{3}", seq)


main()
