from Bio import SeqIO

def longest_complement_palindromes(seq, top_k=10):
    # Complementos
    comp = {'A':'T','T':'A','C':'G','G':'C'}
    n = len(seq)

    # Función para verificar el "match complementario"
    def match(i, j):
        return comp.get(seq[i], 'N') == seq[j]

    # Aquí se guardan los resultados como tuplas (largo, inicio, fin)
    top = []

    def add_result(length, start, end):
        top.append((length, start, end))
        top.sort(reverse=True)  # del más largo al más corto
        if len(top) > top_k:
            top.pop()  # dejamos solo top_k

    # ---- Palíndromos de longitud impar ----
    for center in range(n):
        l = r = center
        while l >= 0 and r < n and match(l, r):
            length = r - l + 1
            add_result(length, l, r+1)
            l -= 1
            r += 1

    # ---- Palíndromos de longitud par ----
    for center in range(n-1):
        l, r = center, center+1
        while l >= 0 and r < n and match(l, r):
            length = r - l + 1
            add_result(length, l, r+1)
            l -= 1
            r += 1

    # Formateo de resultados
    final = []
    for length, start, end in top:
        final.append({
            "largo": length,
            "inicio": start,
            "final": end,
            "secuencia": seq[start:end]
        })
    return final

# -------------main-------------
record = SeqIO.read("sequence.fasta", "fasta")
seq = str(record.seq)

results = longest_complement_palindromes(seq, top_k=10)
for r in results:
    print(r)

