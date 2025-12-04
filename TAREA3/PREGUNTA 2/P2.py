import random

def randomizar_genoma(input_file, output_file):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    # guardamos el header y la secuencia
    header = lines[0].strip()
    secuencia_original = ''.join([line.strip() for line in lines[1:]])
    
    # convertimos a lista para poder mezclar
    lista_nucleotidos = list(secuencia_original)
    random.shuffle(lista_nucleotidos)
    secuencia_random = ''.join(lista_nucleotidos)
    
    with open(output_file, 'w') as f:
        f.write(header + '\n')
        # decidí escribir la secuencia en bloques de 70 caracteres para legibilidad
        ancho_linea = 70
        for i in range(0, len(secuencia_random), ancho_linea):
            f.write(secuencia_random[i:i+ancho_linea] + '\n')

# main
randomizar_genoma('plasmidoEColi.fna', 'plasmidoEColi_randomized.fna')