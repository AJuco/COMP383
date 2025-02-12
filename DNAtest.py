s = open('rosalind_dna.txt','r').read()
for n in ["A","C","G","T"]:
    print(s.count(n), end=' ')
