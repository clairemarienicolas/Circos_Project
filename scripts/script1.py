'''
Usage : 
'''
#Pour lancer facilement depuis TP_Circos/ : python scripts/script1.py data/Sbicolor_313_v3.0.hardmasked.fa config/karyotype.txt
# link to the input data : /home/genouest/mob/asalmon/TP_EDG/data/Os/Osativa_323_v7.0.hardmasked.fa

import sys

def get_input_filename():
    """extract the input filename from the command line"""
    return sys.argv[1]

def get_output_filename():
    """extract the input filename from the command line"""
    return sys.argv[2]

def transform(input_file, output_file):
    """read the """
    #lecture du fichier fasta
    
    chromosomes = {}
    current_chr = None
    couleur = "88,114,107"   # couleur RGB ajoutée dans la dernière colonne

    with open (input_file, "r") as fh:
        for line in fh:
            line = line.strip()

            #si c'est le header fasta avec un chevron 
            if line.startswith(">"):
                header = line[1:].split()[0] #enlève le chrevon et garde le premier "mot"

                #on garde seulement les chromosomes et pas les contigs
                if header.startswith("Chr"):
                    current_chr = header
                    chromosomes[current_chr] = 0
                else:
                    current_chr = None
            
            #si c'est une ligne qui contient une séquence
            elif current_chr is not None:
                chromosomes[current_chr] += len(line)
        
        with open(output_file, "w") as out:
            for chr_name in chromosomes:
                out.write(f"chr - {chr_name} {chr_name} 1 {chromosomes[chr_name]} {couleur}\n")
    
    print(f"Fichier {output_file} créé avec succès.")

def main():
    if len(sys.argv) != 3:
        print("Usage : script1.py input_filename output_filename")
        sys.exit(1)

    input_file = get_input_filename()
    output_file = get_output_filename()

    transform(input_file, output_file)

if __name__ == "__main__":
    main()

