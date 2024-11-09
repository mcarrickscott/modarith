#
# Python script for parsing Wycheproof test vectors for ECDSA and EDDSA
# for example python3 parse.py ed448_test.json
# Output should be C digestible, outputting on separate lines the public key, the comment, the message, the signature and finally the outcome, valid or invalid
#
import sys

# extract substring at position ptr, return substring and updated ptr
def extract(s,ptr,id) :
    pk=s.find(id,ptr)
    fpk=s.find('"',pk+len(id))  # skip "
    lpk=s.find('"',fpk+1)+1
    return s[fpk+1:lpk-1],lpk

inputfile=sys.argv[1]

# read whole file into string
with open(inputfile) as f:     # Wycheproof JSON test vectors
    data = f.read()

ptr=0
finished=False
pubkey=''
if inputfile[0:5]=="ecdsa" :
    pubkey='"uncompressed"' # NOTE: Could be "pk" or "uncompressed"
elif inputfile[0:2]=="ed" :
    pubkey='"pk"'           # NOTE: Could be "pk" or "uncompressed"
else :
    sys.exit(0)

while not finished :
    pk,ptr=extract(data,ptr,pubkey) # ptr points to public key. 
    npk=data.find(pubkey,ptr)   # is there a next public key?
    if npk<0 :
        npk=len(data)   # end of file
# print public key, comments, message, signature and outcome    
    while (True) :
        print(pk)
        ss,ptr=extract(data,ptr,'"comment"')
        print(ss)
        ss,ptr=extract(data,ptr,'"msg"')
        print(ss)
        ss,ptr=extract(data,ptr,'"sig"')
        print(ss)
        ss,ptr=extract(data,ptr,'"result"')
        print(ss)
        nmg=data.find('"comment"',ptr)
        if nmg<0:    # there is no next comment, we are done
            finished=True
            break
        if nmg>=npk: # oops - there is a new public key
            break
