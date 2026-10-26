
import json, glob, os, time
from Bio.PDB import MMCIFParser, NeighborSearch
F = "/kaggle/working/filter409"
HOT_ECD = [32,39,68,91,121,138,141]
HOT_TM  = [152,190,197,205,214,234,298,299,300,310,372,380]
EPITOPE_TSV = "/kaggle/working/glp1r/docs/epitope_contacts.tsv"
epi = set()
for line in open(EPITOPE_TSV).read().splitlines()[1:]:
    c = line.split("\t")
    if c and c[0].isdigit():
        epi.add(int(c[2]))  # uniprot col
def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)
p = MMCIFParser(QUIET=True)
out = {}
done_p = F+"/contacts409.json"
done = json.load(open(done_p)) if os.path.exists(done_p) else {}
def analyze(name, cifs):
    # boltz writes model index 0..; take first
    s = p.get_structure(name, cifs[0])[0]
    chains = {c.id: c for c in s.get_chains()}
    rec = chains.get("A"); binder = chains.get("B")
    if rec is None or binder is None: return {"error": "chains"}
    batoms = [a for a in binder.get_atoms() if a.element != "H"]
    ns = NeighborSearch(batoms)
    hit = set()
    for res in rec.get_residues():
        if res.id[0] != " ": continue
        rn = res.id[1]  # boltz chain A res 1 == UniProt 24
        uni = rn + 23
        for a in res.get_atoms():
            if a.element == "H": continue
            if ns.search(a.coord, 5.0):
                hit.add(uni); break
    return {"n_contacts": len(hit),
            "hot_ecd": sorted(hit & set(HOT_ECD)),
            "hot_tm": sorted(hit & set(HOT_TM)),
            "epitope_overlap": sorted(hit & epi),
            "n_epitope_overlap": len(hit & epi),
            "contacts": sorted(hit)}
changed = False
for name in sorted(set([os.path.basename(os.path.dirname(os.path.dirname(c))) for c in glob.glob(F+"/out*/*/**/*.cif", recursive=True)])):
    pass
for sub in ("out", "out_seed2"):
    for d in sorted(glob.glob(F+f"/{sub}/*")):
        name = os.path.basename(d)
        key = f"{sub}/{name}"
        if key in done: continue
        cifs = glob.glob(d+"/**/*.cif", recursive=True)
        if not cifs: continue
        try:
            done[key] = analyze(key, cifs)
            log("CONTACTS", key, "n", done[key]["n_contacts"], "hotE", len(done[key]["hot_ecd"]), "hotTM", len(done[key]["hot_tm"]), "epi", done[key]["n_epitope_overlap"])
            changed = True
        except Exception as e:
            done[key] = {"error": str(e)[:200]}
            log("FAIL", key, str(e)[:150])
if changed:
    json.dump(done, open(done_p, "w"), indent=1)
log("CONTACTS DONE", len(done))
