import re

def main():
    fasta = "TAIR10_pep_20101214.fa"
    headers, seqs = lees_fasta(fasta)
    # print(headers)
    # print(seqs)
    headers_voldaan, seqs_voldoen = referentie(headers, seqs)
    gff = "TAIR10_GFF3_genes.gff"
    data = lees_gff(gff)
    # print(data)
    info_alles = gff_headers(headers_voldaan, data)

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
    seqs_voldoen = []
    headers_voldoen = []
    for seq in seqs:
        overeenkomst = re.search("[LIVMFYC].[HY].D[LIVMFY]K..N[LIVMFYCT]{3}", seq)
        # print(overeenkomst)
        if overeenkomst is not None:
            seqs_voldoen.append(seq)
            bijbehorende_headers = headers[seqs.index(seq)]
            headers_voldoen.append(bijbehorende_headers[1:12])
    print("Headers:", len(headers_voldoen))
    print("Sequenties:", len(seqs_voldoen))

    return headers_voldoen, seqs_voldoen

def lees_gff(gff):
    """
    Een functie die het hele GFF3 file leest en alle
    info van de genen returnd in vorm van een lijst.
    :param gff: is TAIR10_GFF3_genes.gff
    :return: Een lijst met de data van de genen
    """
    # Afvangen fout bestand
    if gff[-3:] != "gff":
        gff = input("Bestand is geen GFF3 bestand. Geef goed bestand")
        # Roept functie opnieuw aan
        # Om opnieuw het bestand te controleren
        lees_gff(gff)
    bestand = open(gff)
    data = []
    for regel in bestand:
        regel = regel.strip()
        data.append(regel.split("\t"))

    return data

def gff_headers(headers_voldoen, data):
    info_alles = []
    for header in headers_voldoen:
        exonen = 0
        gen_info = []
        for i in data:
            if header in i[8]:
                chr = i[0]
                if i[2] == "exon":
                    exonen += 1
                if i[2] == "mRNA":
                    start = i[3]
                    eind = i[4]
                    lengte_mrna = int(eind) - int(start)

        gen_info.append(header)
        gen_info.append(chr)
        gen_info.append(lengte_mrna)
        gen_info.append(exonen)

        info_alles.append(gen_info)
    print(info_alles)

    return info_alles

main()
