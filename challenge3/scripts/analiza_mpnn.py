GLP1 = 'HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR'  # GLP-1(7-36) referencia oficial del pack
EXENDIN4 = 'HGEGTFTSDLSKQMEEEAVRLFIEWLKNGG'  # exendina-4(1-30) aprox, por si acaso

def max_window_identity(a, b, k=8):
    best, where = 0, None
    for i in range(len(a) - k + 1):
        for j in range(len(b) - k + 1):
            m = sum(1 for x, y in zip(a[i:i+k], b[j:j+k]) if x == y)
            if m > best: best, where = m, (i+1, j+1)
    return best, where

recs = []
name, seq = None, ''
for line in open('__import__("sys").argv[1]'):
    line = line.strip()
    if line.startswith('>'):
        if name: recs.append((name, seq))
        name, seq = line[1:], ''
    else:
        seq += line
if name: recs.append((name, seq))

print(f'{"sample":8s} {"score":>6s} {"max8mer":>7s} {"pos":>6s}  B3 B6 B28')
for n, s in recs[1:]:
    gs = [p.split('=')[1] for p in n.split(',') if 'global_score' in p][0]
    sid = [p.split('=')[1] for p in n.split(',') if 'sample' in p][0]
    m8, w = max_window_identity(s, GLP1)
    print(f'{sid:8s} {gs:6s} {m8}/8      {w}  {s[2]}  {s[5]}  {s[27]}')
