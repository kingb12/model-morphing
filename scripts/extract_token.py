import sys
with open(sys.argv[1], 'r') as f:
    lines = f.readlines()
token = None
for line in lines:
    if line.startswith('token'):
        token = line.split('token=')[-1]
with open(sys.argv[2], 'w') as f:
    if token is not None:
        f.write(token)