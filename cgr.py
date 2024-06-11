from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def lire_genome(fichier):
    # Lire un fichier GenBank et retourner la séquence d'ADN.
    record = SeqIO.read(fichier, "genbank")
    sequence = record.seq
    if len(sequence) < 10000:
        raise ValueError("La séquence doit contenir plus de 10 000 paires de bases.")
    return sequence


def count_kmers(sequence, k):
    # Compte les occurrences de k-mers dans la séquence donnée.
    d = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        if 'N' not in kmer:  # Ignore les k-mers contenant le nucléotide inconnu 'N'
            if kmer in d:
                d[kmer] += 1
            else:
                d[kmer] = 1
    return d


def chaos_game_representation(counts, k):
    # Génère la matrice pour la représentation du jeu du chaos pour les k-mers donnés
    array_size = 2**k  
    chaosmat = np.zeros((array_size, array_size))  # La taille de la matrice est 2^k par 2^k
    
    nucleotide_positions = {'C': (0, 0), 'A': (0, 1), 'T': (1, 1), 'G': (1, 0)}
    
    for kmer, count in counts.items():
        posx, posy = 0.5, 0.5 
        for i, char in enumerate(kmer):
            corner_x, corner_y = nucleotide_positions[char]
            posx += corner_x * (2**(k-(i+1)))
            posy += corner_y * (2**(k-(i+1)))

        # Convertir position float en index
        index_x = int(posx)
        index_y = int(posy)

        # Update la matrice
        chaosmat[index_y, index_x] += count

    return chaosmat


def plot_chaos_game(chaos_matrix_done, title):
    plt.figure(figsize=(8, 8))
    plt.title(title)
    plt.imshow(chaos_matrix_done, cmap='gray_r', origin='lower')
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    fichier = input("Entrez le nom du fichier GenBank : ")
    k = int(input("Entrez la taille souhaitée pour les k-mers (1-6) : "))
    sequence = lire_genome(fichier)
    kmer_counts = count_kmers(sequence, k)
    chaos_matrix_done = chaos_game_representation(kmer_counts, k)
    plot_chaos_game(chaos_matrix_done, f"CGR pour {fichier} et k-mers de taille {k}")