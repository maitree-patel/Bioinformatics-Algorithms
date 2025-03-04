import pandas as pd

# fasta Parser: a custom function that parses a fasta file to separate out and store the header and sequence information in a dictionary for downstream analysis
# The parameter include:
    # fastaFile: file path of the fasta file
def fastaParser(fastaFile):
    
    # Opening a fasta file
    file = open(fastaFile)

    # Using the read method to get the contents of the file as a single string
    file_str = file.read()

    # Splitting to get contents per instance(i.e. each header and sequencce)
    seq_lst = file_str.split(">")
    
    # Looping through each instance to get the header and the sequence
    fastaDict = {}

    for seq in seq_lst:
        split1 = seq.split("\n")
        header = split1[0]
        print(header)
        sequence = "".join(split1[1:])
        print(sequence)
        fastaDict[header] = sequence
    
    firstEmptyKey = next(iter(fastaDict.keys()))
    print(firstEmptyKey)
    del fastaDict[firstEmptyKey]

    # Creating a dataframe of the fasta file
    fastaDf = pd.DataFrame(list(fastaDict.items()), columns=['header', 'sequence'])
    return fastaDf
    
fastaFile =  "seq2.fa"
#testing to see if we get a SequenceRecord of the first header and sequence from the fasta file 
fastaParsed = fastaParser(fastaFile)
print(fastaParsed)

# Example GC content calculation from parsed data
GCperc = [f"{round(((row['sequence'].count('g')+row['sequence'].count('c'))/(len(row['sequence'])))*100, 2)}%" for index, row in fastaParsed.iterrows()]
print(f"The GC% of fasta sequences {list(fastaParsed['sequence'])} is {GCperc}")
