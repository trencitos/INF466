from Bio import SeqIO
from Bio.Seq import Seq

def calcular_gc(secuencia_adn):
    """calcula el porcentaje de G+C de una secuencia de nucleótidos."""
    # convertimos a string y mayúsculas para evitar errores
    sec = str(secuencia_adn).upper()
    g = sec.count('G')
    c = sec.count('C')
    length = len(sec)
    if length == 0: return 0
    return ((g + c) / length) * 100

def encontrarOrfs(archivo_fasta, min_pro_len=100):
    """
    busca ORFs en un archivo FASTA, retornando una lista de diccionarios.
    considera ORFs que inician con 'M' y terminan con Stop (*).
    """
    resultados = []
    rec = SeqIO.read(archivo_fasta, "fasta")
    secuencia_total = rec.seq

    # tabla genética bacteriana (11) es estándar para E. Coli, 
    tabla = 11 
    for hebra, nuc in [(1, secuencia_total), (-1, secuencia_total.reverse_complement())]:
        
        for frame in range(3):
            seq_frame = nuc[frame:]
            sobrante = len(seq_frame) % 3
            if sobrante > 0:
                seq_frame = seq_frame[:-sobrante]

            trans = seq_frame.translate(table=tabla)
            
            # lógica para encontrar ORFs: Buscar M ... *
            aa_len = len(trans)
            aa_start = 0
            
            while aa_start < aa_len:
                # buscamos la primera 'M' = ATG (nuestro start codon)
                aa_start = trans.find("M", aa_start)
                if aa_start == -1: break
                
                # buscamos el siguiente Stop '*'
                aa_end = trans.find("*", aa_start)
                
                if aa_end != -1:
                    # encontramos un ORF válido (M...*)
                    largo_proteina = aa_end - aa_start + 1 # +1 para incluir el stop si lo queremos
                    
                    nuc_start = frame + (aa_start * 3)
                    nuc_end = frame + (aa_end * 3) + 3 
                    
                    orf_dna = nuc[nuc_start:nuc_end]
                    
                    # criterio mínimo de longitud (opcional, puse 105 nucleótidos en este caso)
                    if len(orf_dna) >= 105: 
                        gc_percent = calcular_gc(orf_dna)
                        
                        resultados.append({
                            "secuencia_orf": str(orf_dna),
                            "largo_nuc": len(orf_dna),
                            "gc_porcentaje": gc_percent,
                            "hebra": "Forward" if hebra == 1 else "Reverse",
                            "frame": frame + 1
                        })
                    
                    # avanzamos el buscador después de este ORF para buscar otros anidados o siguientes
                    aa_start += 1 
                else:
                    break # esto significa que no hay cierre para esta M

    return resultados

# --- main ---
try:
    archivo = "plasmidoEColi.fna"
    todos_orfs = encontrarOrfs(archivo)

    # orden de mayor a menor
    orfs_ordenados = sorted(todos_orfs, key=lambda x: x["gc_porcentaje"], reverse=True)

    top_3_alto_gc = orfs_ordenados[:3]
    top_3_bajo_gc = orfs_ordenados[-3:]

    print(f"--- Análisis de ORFs para {archivo} ---\n")

    print("=== Top 3 %GC más ALTO ===")
    for i, orf in enumerate(top_3_alto_gc, 1):
        print(f"#{i} | Longitud: {orf['largo_nuc']} bp | GC: {orf['gc_porcentaje']:.2f}%")
        print(f"    Secuencia: {orf['secuencia_orf']}...") 

    print("\n=== Top 3 %GC más BAJO ===")
    for i, orf in enumerate(top_3_bajo_gc, 1):
        print(f"#{i} | Longitud: {orf['largo_nuc']} bp | GC: {orf['gc_porcentaje']:.2f}%")
        print(f"    Secuencia: {orf['secuencia_orf']}...") 

except FileNotFoundError:
    print(f"Error: No se encontró el archivo '{archivo}'. Debe estar en la misma carpeta!")
