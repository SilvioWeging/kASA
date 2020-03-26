#taken from https://blog.softhints.com/python-convert-json-to-json-lines/

import json, sys

with open(sys.argv[1], 'r') as f:
    json_data = json.load(f)
    
with open(sys.argv[2], 'w') as outfile:
    for entry in json_data:
        json.dump(entry, outfile)
        outfile.write('\n')