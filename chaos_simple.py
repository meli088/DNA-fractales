import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

def lire_genome(file_path):
    # Lecture d'un fichier GenBank
    record = SeqIO.read(file_path, "genbank")
    sequence = str(record.seq)
    if len(sequence) < 10000:
        raise ValueError("La séquence doit contenir plus de 10000 bases pour être valide. ")
    return sequence

def chaos_game_simple(sequence):
    # Attribution des nucléotides à chaque sommets
    nucleotide_coords = {'C': (0, 0), 'A': (0, 1), 'T': (1, 1), 'G': (1, 0)}
    
    # Positionnement au milieu du carré
    x, y = 0.5, 0.5
    coordonnees = [(x, y)]
    
    # Mise à jour et stockage des coordonnées
    for nuc in sequence:
        if nuc in nucleotide_coords:
            coin_x, coin_y = nucleotide_coords[nuc]
            x = (x + coin_x) / 2
            y = (y + coin_y) / 2
            coordonnees.append((x, y))
    
    return coordonnees

def plot_chaos_game(coordonnees, title):
    # Récupération des coordonnées pour pouvoir placer les point
    i, j = zip(*coordonnees)
    
    # Placement des points et traçage d'un carré pour la représentation
    plt.figure(figsize=(8, 8))
    plt.scatter(i, j, color='black', s=0.02)  # Utilisation d'une taille de point très faible
    plt.title(title)
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.show()

if __name__ == "__main__":
    fichier = input("Entrez le nom du fichier GenBank : ")
    sequence = lire_genome(fichier)
    coordonnees = chaos_game_simple(sequence)
    plot_chaos_game(coordonnees, f"Représentation en du Chaos Game Simple pour {fichier}")
