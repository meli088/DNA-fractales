Liste des commandes modules à utiliser pour installer les modules nécessaire au lancement du code :
Pour le module Bio :
pip install biopython 
ou 
conda install conda-forge::biopython

Pour le module scipy :
python3 -m pip install scipy 
ou 
conda install scipy

Pour le module matplotlib :
pip install matplotlib 
ou 
conda install conda-forge::matplotlib

1) Pour lancer le chaos game simple dans un terminal :
python3 chaos_simple.py
entrer le nom du fichier avec son extension
Exemple :
python3 chaos_simple.py
Entrez le nom du fichier GenBank: NC_001133.gbk
(! Mettre la fenêtre du graphique en plein écran pour les séquences de grandes tailles !)


2) Pour lancer la CGR aller dans un terminal :
python3 cgr.py
Et entrer le nom du fichier avec son extension,
Enfin entrer la taille des k-mers souhaitée.

Par exemple : 
python3 cgr.py
Entrez le nom du fichier GenBank: NC_001133.gbk
Entrez la taille souhaitée pour les k-mers (1-6): 4

3) Pour effectuer un dendogramme :
python3 dendrogramme.py
Entrer le nombre d'échantillons que vous voulez entrer
Entrer la taille des k-mers souhaitée
Enfin, pour le nombre entré au-dessus, le terminal vous demandera le nom d'un fichier GenBank

Exemple :
python3 dendrogramme.py
Entrez le nombre de séquences à analyser : 3
Entrez la taille souhaitée pour les k-mers (1-6) : 3
Entrez le nom du fichier GenBank : BA000017.gbk
Entrez le nom du fichier GenBank : NC_001133.gbk
Entrez le nom du fichier GenBank : NC_004353.gbk
