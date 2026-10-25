#!/usr/bin/env python3
"""
prepare_target.py -- Challenge 3 (GLP1R) target preparation. CPU-only, stdlib-only.

Reads  ../target/6X18.cif  (GLP-1 peptide hormone bound to GLP-1 receptor, cryo-EM 2.1 A;
chains: A=Galpha-s, B=Gbeta, G=Ggamma, N=Nanobody35, P=GLP-1(7-36)NH2, R=GLP1R).

Produces:
  ../target/6X18_R.pdb        receptor chain R, model 1 (auth numbering kept)
  ../target/6X18_RP.pdb       receptor + peptide chains R/P, model 1
  ../docs/epitope_contacts.tsv  per-residue contacts R<->P (heavy atoms <= 4.5 A),
                                with UniProt P43220 numbering mapped by sequence offset
Prints a summary: observed ranges, unmodeled gaps, ECD/TM contact split.

Numbering: the deposited construct adds a signal peptide + FLAG tag
(MKTIIALSYIFCLVFADYKDDDDLEVLFQGPAR) before UniProt residue 24, so auth_seq_id
is offset from UniProt. The script detects the offset by matching the UniProt
segment and reports BOTH numberings. Verify the printed offset before trusting
UniProt numbers.
"""
import sys, os, math

HERE = os.path.dirname(os.path.abspath(__file__))
CIF = os.path.join(HERE, '..', 'target', '6X18.cif')
UNIPROT_SEG  = os.path.join(HERE, '..', 'target', 'glp1r_P43220_24-432.fasta')   # base 24
UNIPROT_FULL = os.path.join(HERE, '..', 'target', 'P43220.fasta')                  # base 1 (UniProt download)
OUT_R = os.path.join(HERE, '..', 'target', '6X18_R.pdb')
OUT_RP = os.path.join(HERE, '..', 'target', '6X18_RP.pdb')
OUT_TSV = os.path.join(HERE, '..', 'docs', 'epitope_contacts.tsv')
CONTACT_CUTOFF = 4.5
ECD_UNIPROT_END = 145   # N-terminal domain ~24-145 (~120 aa per the challenge page); loops beyond are ECLs

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def read_uniprot():
    """Returns (sequence, base): base = UniProt residue number of sequence[0].
    Prefers the challenge-page segment (24-432); falls back to the full-length
    UniProt download so the DGX bootstrap works with either file present."""
    for path, base in ((UNIPROT_SEG, 24), (UNIPROT_FULL, 1)):
        if os.path.exists(path):
            return ''.join(l.strip() for l in open(path) if not l.startswith('>')), base
    raise SystemExit('no UniProt FASTA found in ../target/ (segment or full P43220.fasta)')

def parse_mmcif_atoms(path):
    """Minimal atom_site loop parser. Returns list of dicts for model 1."""
    lines = open(path).read().splitlines()
    # locate atom_site column definitions
    cols, i = [], 0
    while i < len(lines):
        if lines[i].startswith('_atom_site.'):
            cols.append(lines[i].split()[0].split('.')[1]); i += 1
        elif cols:
            break
        else:
            i += 1
    idx = {c: n for n, c in enumerate(cols)}
    need = ['group_PDB','label_atom_id','label_comp_id','label_asym_id','label_seq_id',
            'Cartn_x','Cartn_y','Cartn_z','auth_seq_id','auth_comp_id','auth_asym_id',
            'auth_atom_id','pdbx_PDB_model_num']
    for n in need:
        if n not in idx: sys.exit(f'missing mmCIF column {n}')
    atoms = []
    for ln in lines[i:]:
        if not (ln.startswith('ATOM') or ln.startswith('HETATM')):
            if atoms: break
            continue
        f = ln.split()
        if len(f) < len(cols):  # multiline value; not expected in atom_site rows of this file
            continue
        if f[idx['pdbx_PDB_model_num']] != '1':
            continue
        atoms.append({
            'group': f[idx['group_PDB']],
            'atom': f[idx['auth_atom_id']], 'comp': f[idx['auth_comp_id']],
            'chain': f[idx['auth_asym_id']], 'auth_seq': f[idx['auth_seq_id']],
            'x': float(f[idx['Cartn_x']]), 'y': float(f[idx['Cartn_y']]), 'z': float(f[idx['Cartn_z']]),
        })
    return atoms

def write_pdb(atoms, chains, path):
    n = 0
    with open(path, 'w') as fh:
        for a in atoms:
            if a['chain'] not in chains: continue
            n += 1
            atom_field = (' ' + a['atom']).ljust(4)[:4]   # PDB cols 13-16: one-letter elements start col 14
            fh.write('ATOM  %5d %s %3s %s%4s    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n' % (
                n, atom_field, a['comp'], a['chain'], a['auth_seq'], a['x'], a['y'], a['z'], a['atom'][0]))
        fh.write('END\n')
    return n

def main():
    for out in (OUT_R, OUT_RP, OUT_TSV):
        os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    uni, uni_base = read_uniprot()
    atoms = parse_mmcif_atoms(CIF)
    r = [a for a in atoms if a['chain'] == 'R']
    p = [a for a in atoms if a['chain'] == 'P']
    if not r or not p: sys.exit('chain R or P missing')
    n_r = write_pdb(atoms, {'R'}, OUT_R)
    n_rp = write_pdb(atoms, {'R','P'}, OUT_RP)

    # observed receptor residues in order, for offset mapping
    seen, order = set(), []
    for a in r:
        if a['auth_seq'] not in seen:
            seen.add(a['auth_seq'])
            order.append((a['auth_seq'], AA3.get(a['comp'], 'X')))
    obs_seq = ''.join(s for _, s in order)
    # find offset: UniProt position (1-based within the 24-432 segment) of first observed residue
    probe = obs_seq[:10]
    pos = uni.find(probe)
    if pos < 0:
        offset_note = 'OFFSET NOT FOUND: first observed motif %r not in UniProt segment; UniProt numbers disabled' % probe
        uni_offset = None
    else:
        first_auth = int(order[0][0])
        uni_offset = (uni_base + pos) - first_auth   # uniprot = auth + uni_offset
        offset_note = f'offset check: auth {first_auth} = UniProt {uni_base+pos} (uni = auth {uni_offset:+d}); fasta base {uni_base}; motif {probe}'

    # contacts R<->P (heavy atoms only)
    ph = [a for a in p if a['atom'][0] != 'H' and a['comp'] in AA3]  # amino acids only: waters/ions share chain P in this deposit
    rh = {}
    for a in r:
        if a['atom'][0] == 'H' or a['comp'] not in AA3: continue
        rh.setdefault(a['auth_seq'], {'comp': a['comp'], 'atoms': []})['atoms'].append(a)
    contacts = {}
    for pa in ph:
        pr = f"{pa['comp']}{pa['auth_seq']}"
        if pa['comp'] not in AA3: continue
        for auth_seq, res in rh.items():
            for ra in res['atoms']:
                d2 = (pa['x']-ra['x'])**2 + (pa['y']-ra['y'])**2 + (pa['z']-ra['z'])**2
                if d2 <= CONTACT_CUTOFF**2:
                    d = math.sqrt(d2)
                    cur = contacts.get(auth_seq)
                    if cur is None or d < cur[1]:
                        contacts[auth_seq] = (pr, d, res['comp'])
    rows = []
    for auth_seq in sorted(contacts, key=lambda s: int(s)):
        pr, d, comp = contacts[auth_seq]
        uni_n = int(auth_seq) + uni_offset if uni_offset is not None else None
        region = ''
        if uni_n is not None:
            region = 'ECD' if uni_n <= ECD_UNIPROT_END else 'TM/ECL'
        rows.append((auth_seq, comp, uni_n if uni_n is not None else '?', region or '?', pr, f'{d:.2f}'))
    with open(OUT_TSV, 'w') as fh:
        fh.write('auth_seq\tres\tuniprot\tregion\tpeptide_residue\tmin_dist_A\n')
        for row in rows:
            fh.write('\t'.join(map(str, row)) + '\n')

    ecd = [x for x in rows if x[3] == 'ECD']
    tm = [x for x in rows if x[3] == 'TM/ECL']
    obs = [int(s) for s, _ in order]
    print(f'atoms model1: total={len(atoms)} R={len(r)} P={len(p)}')
    print(f'wrote {OUT_R} ({n_r} atoms), {OUT_RP} ({n_rp} atoms), {OUT_TSV} ({len(rows)} contact residues)')
    print(f'receptor observed auth range: {min(obs)}..{max(obs)} ({len(obs)} residues)')
    print(offset_note)
    print(f'contact residues (<= {CONTACT_CUTOFF} A): ECD={len(ecd)}  TM/ECL={len(tm)}')
    if uni_offset is not None:
        pep_partners = sorted({x[4] for x in rows})
        print('peptide residues in contact:', ', '.join(pep_partners))
        print('first/last UniProt contact:', min(x[2] for x in rows), '-', max(x[2] for x in rows))

if __name__ == '__main__':
    main()
