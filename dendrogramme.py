import scipy.cluster.hierarchy as sci
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt

def lire_genome(fichier):
    # Lire un fichier GenBank et retourner la séquence d'ADN.
    record = SeqIO.read(fichier, "genbank")
    sequence = str(record.seq)
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
        posx, posy = 0, 0 
        for i, char in enumerate(kmer):
            corner_x, corner_y = nucleotide_positions[char]
            posx += corner_x * (2**(k-(i+1)))
            posy += corner_y * (2**(k-(i+1)))

        # Convertir position float en index
        index_x = int(posx) % array_size
        index_y = int(posy) % array_size

        # Update la matrice
        chaosmat[index_y, index_x] += count

    return chaosmat

def dendogramme(chaos_matrix, samples):
    #Calcul des distances et exploitation de celles-ci pour la modélisation de l'arbre
    distance = pdist(chaos_matrix)
    arbre = linkage(distance)
    
    #Affichage de l'arbre
    plt.figure(figsize=(8, 8))
    dendrogram(arbre, labels=samples)
    plt.title('Dendrogramme de clustering hiérarchique')
    plt.xlabel('Échantillons')
    plt.ylabel('Distance')
    plt.show()

if __name__ == "__main__":
    chaos_vectors = []  # Liste pour stocker les vecteurs de chaos
    samples = []  # Liste pour les noms des échantillons
    nb_sequence = int(input("Entrez le nombre de séquences à analyser : "))
    k = int(input("Entrez la taille souhaitée pour les k-mers (1-6) : "))
    
    for i in range(nb_sequence):
        fichier = input("Entrez le nom du fichier GenBank : ")
        try:
            sequence = lire_genome(fichier)
        except Exception as e:
            print(f"Erreur de lecture du fichier {fichier}: {e}")
            continue
        
        kmer_counts = count_kmers(sequence, k)
        chaos_matrix_done = chaos_game_representation(kmer_counts, k)
        chaos_vector = chaos_matrix_done.flatten()
        
        # Ajouter le vecteur de chaos à la liste
        chaos_vectors.append(chaos_vector)
        
        # Ajouter le nom du fichier à la liste des échantillons
        samples.append(fichier)
    
    # Convertir la liste des vecteurs de chaos en une matrice NumPy
    big_matrix = np.array(chaos_vectors)
    
    # Créer le dendrogramme
    dendogramme(big_matrix, samples)
