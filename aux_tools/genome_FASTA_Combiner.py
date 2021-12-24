import gzip
import glob

files = 'Escherichia_coli*.dna.toplevel.fa_UR.fasta.gz.fasta.gz'
combined_out = gzip.open('/home/nick/Desktop/Ensem/UR-Con-StORFs/E-coli_UR-Con-StORFS_Combined_DNA.fasta.gz', 'wb')#, newline='\n', encoding='utf-8')
count = 0

for file in glob.glob('/home/nick/Desktop/Ensem/UR-Con-StORFs/'+files):
    count +=1
    with gzip.open(file,'rb') as genome:
        print(file)
        for line in genome:
            if line.startswith(b'#'):
                continue
            elif line.startswith(b'>'): # Might need to change file[?]
                #genome = bytes(file.split('/')[2].split('.cds.all.fa.gz')[0],'utf-8')
                genome = bytes(file.split('/')[6].split('.dna.toplevel.fa_UR.fasta.gz.fasta.gz')[0], 'utf-8')
                genome = genome.capitalize()
                line = line.replace(b'>',b'>'+genome+b'|')
                combined_out.write(line)
            else:
                combined_out.write(line)
