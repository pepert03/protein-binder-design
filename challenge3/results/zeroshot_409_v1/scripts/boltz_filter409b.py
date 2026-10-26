
import json, os, subprocess, time, glob
W = "/kaggle/working/glp1r"
F = "/kaggle/working/filter409"
def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)
def collect(outdir):
    out = {}
    for f in glob.glob(outdir+"/**/*.json", recursive=True):
        if "confidence" in os.path.basename(f):
            d = json.load(open(f))
            out = {k: d[k] for k in ("confidence_score","ptm","iptm","complex_plddt") if k in d}
    return out
def predict(name, yaml_p, outdir, seed=None):
    cmd = f"/kaggle/working/boltzenv/bin/boltz predict {yaml_p} --out_dir {outdir} --cache /root/boltzcache --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 --devices 1"
    if seed is not None:
        cmd += f" --seed {seed}"
    t0=time.time()
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    dt=int(time.time()-t0)
    if r.returncode != 0:
        log("FAIL", name, r.stderr[-300:])
        return {"rc": r.returncode}
    m = collect(outdir)
    log("PRED", outdir.split("/")[-2], name, "t", dt, m)
    return {"rc": 0, "t": dt, **m}
# 1) seed-repeat: top-3 designs by iptm + controls
res = json.load(open(F+"/results.json"))
designs = {k:v for k,v in res.items() if k.startswith(("gA_","gB_")) and v.get("iptm") is not None}
top3 = sorted(designs, key=lambda k: -designs[k]["iptm"])[:3]
targets = top3 + ["nativo", "aleatorio36"]
log("seed-repeat targets:", targets)
s2 = json.load(open(F+"/results_seed2.json")) if os.path.exists(F+"/results_seed2.json") else {}
for name in targets:
    if name in s2: continue
    yp = F+f"/inputs/f409_{name}.yaml"
    s2[name] = {"seq": res[name]["seq"], "seed": 1234, **predict(name, yp, F+f"/out_seed2/{name}", seed=1234)}
    json.dump(s2, open(F+"/results_seed2.json","w"), indent=1)
log("seed-repeat done")
# 2) keep watching for new gA designs (primary-seed predict into out/ + results.json)
done = json.load(open(F+"/done.json"))
while True:
    designs_now = {}
    for tag in ("gA","gB"):
        p = W+f"/{tag}_results_incr.json"
        if os.path.exists(p):
            for r in json.load(open(p)):
                designs_now[f"{tag}_{r['traj']}"] = r["seq"]
    pending = [k for k in designs_now if k not in done]
    if not pending:
        log("idle; waiting 300s")
        time.sleep(300)
        continue
    for name in pending:
        seq = designs_now[name]
        yp = F+f"/inputs/f409_{name}.yaml"
        if not os.path.exists(yp):
            recA = open(F+"/inputs/f409_nativo.yaml").read()
            # build from nativo yaml, replacing chain B sequence
            nat_seq = "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR"
            open(yp,"w").write(recA.replace(nat_seq, seq))
        out = predict(name, yp, F+f"/out/{name}")
        if out.get("rc")==0:
            res[name] = {"seq": seq, "t_s": out["t"], **{k:v for k,v in out.items() if k in ("confidence_score","ptm","iptm","complex_plddt")}}
            json.dump(res, open(F+"/results.json","w"), indent=1)
        done[name] = {"rc": out.get("rc"), "t": out.get("t")}
        json.dump(done, open(F+"/done.json","w"), indent=1)
