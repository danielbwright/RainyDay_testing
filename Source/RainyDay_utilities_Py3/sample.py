import json

parameterfile = "/Volumes/TheCordex/Madison_Data/test.json"

with open(parameterfile,'r') as f:
    cardinfo = json.loads(f.read())

print(cardinfo)