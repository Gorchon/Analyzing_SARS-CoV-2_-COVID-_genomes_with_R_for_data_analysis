#Cargar las librerías
library(Biostrings);
library(DECIPHER);
library(ade4);
library(seqinr);
library(adegenet);
library(ape);
library(ggtree);
library(viridis);
library(ggplot2);
# Leer genomas de virus
# covid muestra 1
print("Parte 1, 2 y 3 ejecutandose");
first <- read.GenBank("MN908947.3")
cat("Estructura del objeto DNAbin para el primer ejemplar de covid:\n")
str(first)
cat("\n")

# covid muestra 2
second <- read.GenBank("MT215194")
cat("Estructura del objeto DNAbin para el segundo ejemplar de covid:\n")
str(second)
cat("\n")

# covid muestra 3
third <- read.GenBank("MT947570")
cat("Estructura del objeto DNAbin para el tercer ejemplar de covid:\n")
str(third)
cat("\n")

# covid muestra 4
fourth <- read.GenBank("LR963142")
cat("Estructura del objeto DNAbin para el cuarto ejemplar de covid:\n")
str(fourth)
cat("\n")

# covid muestra 5
fifth <- read.GenBank("LR994660")
cat("Estructura del objeto DNAbin para el quinto ejemplar de covid:\n")
str(fifth)
cat("\n")

# covid muestra 6
sixth <- read.GenBank("OL687293")
cat("Estructura del objeto DNAbin para sexto ejemplar de covid:\n")
str(sixth)
cat("\n")

# covid muestra 7
sept <- read.GenBank("ON072482")
cat("Estructura del objeto DNAbin para el septimo ejemplar de covid:\n")
str(sept)
cat("\n")

# covid muestra 8
eight <- read.GenBank("OK339336")
cat("Estructura del objeto DNAbin para octavo ejemplar de covid:\n")
str(eight)
cat("\n")

# covid muestra 9
nine <- read.GenBank("BS007015")
cat("Estructura del objeto DNAbin para el noveno ejemplar del covid:\n")
str(nine)
cat("\n")



#parte 4
print("Parte 4 ejecutandose");
# Cargar el paquete ape
library(ape)

# Concatenar todas las secuencias en un objeto DNAbin
all_sequences <- c(first,second,third,fourth, fifth, sixth, sept , eight , nine )

# Escribir las secuencias en un archivo FASTA
write.dna(all_sequences, file = "all_genomes.fasta", format = "fasta")

#parte 5
print("Parte 5 ejecutandose");
# Cargar el paquete Biostrings


# Leer el archivo FASTA con las secuencias concatenadas
all_sequences <- readDNAStringSet("all_genomes.fasta")

# Imprimir el contenido del objeto DNAStringSet
all_sequences



#parte 6
print("Parte 6 ejecutandose");
# Cargar el paquete DECIPHER


# Orientar los nucleótidos de las secuencias
oriented_sequences <- OrientNucleotides(all_sequences)

# Imprimir el contenido del objeto DNAStringSet orientado
oriented_sequences

#parte 7
print("Parte 7 ejecutandose");
# Alinear las secuencias orientadas
aligned_sequences <- AlignSeqs(oriented_sequences)

# Visualizar el resultado del alineamiento en el navegador
BrowseSeqs(aligned_sequences)

#parte 8
print("Parte 8 ejecutandose");
# Cargar el paquete Biostrings
library(Biostrings)

# Guardar el resultado del alineamiento en un archivo en formato .fasta
writeXStringSet(aligned_sequences, file = "aligned_sequences.fasta", format = "fasta")

#parte 9
print("Parte 9 ejecutandose");
# Cargar el paquete seqinr
library(seqinr)

# Cargar el archivo .fasta
aligned_sequences <- read.alignment("aligned_sequences.fasta", format = "fasta")

# Mostrar el contenido del objeto tipo alignment
aligned_sequences

#parte 10
print("Parte 10 ejecutandose");

# Cargar el paquete seqinr
library(seqinr)
library(ape)

# Cargar el archivo .fasta
aligned_sequences <- read.alignment("aligned_sequences.fasta", format = "fasta")

# Cargar las secuencias
sequences <- read.fasta("aligned_sequences.fasta")

# Crear una matriz de distancia con la función dist.alignment de seqinr
dist_matrix <- dist.alignment(sequences)

# Mostrar la matriz de distancia
print(dist_matrix)

# Obtener una tabla en escala de grises con la función table.paint del paquete ape
# Obtener una tabla en escala de grises con la función table.paint del paquete ape
gray_table <- table.paint(as.matrix(dist_matrix))

# Mostrar la imagen de la tabla en escala de grises
plot(gray_table)

#parte 11
print("Parte 11 ejecutandose");


library(ape)

# Construir el objeto phylo con la función nj
phylo_tree <- nj(as.dist(dist_matrix))

# Mostrar el contenido del objeto phylo
phylo_tree


#parte 12
print("Parte 12 ejecutandose");
# Cargar el paquete seqinr
library(seqinr)

# Cargar el paquete ape
library(ape)

# Cargar el archivo .fasta
aligned_sequences <- read.alignment("aligned_sequences.fasta", format = "fasta")

# Crear una matriz de distancia con la función dist.alignment de seqinr
dist_matrix <- dist.alignment(aligned_sequences)

# Obtener una tabla en escala de grises con la función table.paint del paquete ape
gray_table <- table.paint(as.matrix(dist_matrix), col = grey.colors(256))


# Construir un objeto de tipo phylo con la función nj del paquete ape
tree <- nj(dist_matrix)

# Mostrar el árbol filogenético
plot(tree)


