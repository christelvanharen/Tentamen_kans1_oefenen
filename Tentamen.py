import re
import numpy as np
import matplotlib.pyplot as plt
import tkinter
from tkinter.ttk import Combobox

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
    grafiek(info_alles)
    Gui(info_alles)

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
    # print(info_alles)

    return info_alles

def grafiek(info_alles):
    x_as = []
    y_as = []
    for onderdeel_lijst in info_alles:
        # print(onderdeel_lijst[1])
        if onderdeel_lijst[1] not in x_as:
            # print(onderdeel_lijst[1])
            x_as.append(onderdeel_lijst[1])
        y_as.append(onderdeel_lijst[1])
    print("X-as", x_as)
    print("Y-as", y_as)

    chr1 = y_as.count("Chr1")
    chr2 = y_as.count("Chr2")
    chr3 = y_as.count("Chr3")
    chr4 = y_as.count("Chr4")
    chr5 = y_as.count("Chr5")

    print(chr1, chr2, chr3, chr4, chr5)


    # Make a fake dataset:
    height = [chr1, chr2, chr3, chr4, chr5]
    y_pos = np.arange(len(x_as))

    # Create bars
    plt.bar(y_pos, height)

    # Create names on the x-axis
    plt.xticks(y_pos, x_as)
    plt.title("Het aantal genen met Serine/Threonine kinase active "
              "site per chromosoom")
    plt.xlabel("Chromosoom nummer")
    plt.ylabel("Aantal genen met Serine/Threonine kinase active site")

    # Show graphic
    plt.show()

class Gui:

    def __init__(self, info_alles):
        # maakt een window en geeft een titel
        self.lijst = info_alles

        self.main_window = tkinter.Tk()
        self.main_window.title("Uitleg")

        # maakt Frames ( gedeeltes in de main_window).
        self.top_frame = tkinter.Frame(self.main_window)
        self.mid_frame = tkinter.Frame(self.main_window)
        self.bottom_frame = tkinter.Frame(self.main_window)

        #lijst in gui

        self.drop = Combobox(self.bottom_frame,
                             values = self.lijst)

        self.drop.set("Select")

        self.drop.pack()

        # knop die functie doe iets aanroept
        self.ok_button = tkinter.Button(self.bottom_frame,
                                        text='OK',
                                        command=self.doe_iets)
        self.ok_button.pack()

        # pack Frames (zichtbaar maken frames in main window)
        self.top_frame.pack()
        self.mid_frame.pack()
        self.bottom_frame.pack()

        #zorg dat Gui 'aan' gaat
        tkinter.mainloop()

    def doe_iets(self):
        entry = self.drop.get()
        print(entry)

main()
