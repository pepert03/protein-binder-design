#!/usr/bin/env python3
"""Comprobador de recursos previo a Boltz-2 (CPU, sin GPU, sin descargas).
Ejecutar EN la sesion destino antes de instalar/descargar nada.
Imprime: limite RAM efectivo (cgroup v1/v2 si existe, si no MemTotal),
RAM disponible, cores, disco libre en CWD y $HOME, GPU visible.
NO infiere viabilidad: solo mide y reporta. La decision es humana.
"""
import os, shutil, sys

def cgroup_mem_limit():
    # cgroup v2
    for p in ('/sys/fs/cgroup/memory.max', '/sys/fs/cgroup/memory/memory.limit_in_bytes'):
        try:
            v = open(p).read().strip()
            if v and v != 'max':
                n = int(v)
                if 0 < n < 1 << 60:
                    return n, p
        except FileNotFoundError:
            pass
    return None, None

def meminfo(key):
    for l in open('/proc/meminfo'):
        if l.startswith(key):
            return int(l.split()[1]) * 1024
    return None

def gib(n): return f"{n/2**30:.2f} GiB" if n is not None else "?"

lim, src = cgroup_mem_limit()
total = meminfo('MemTotal:'); avail = meminfo('MemAvailable:')
eff = min(x for x in (lim, total) if x) if (lim or total) else None
cwd_free = shutil.disk_usage('.').free
home_free = shutil.disk_usage(os.path.expanduser('~')).free
gpu = os.path.exists('/proc/driver/nvidia/gpus') or shutil.which('nvidia-smi') is not None

print("== check_resources (Boltz-2 previo) ==")
print(f"cgroup memory limit : {gib(lim)}" + (f"  (fuente: {src})" if src else "  (no hay limite cgroup o no legible)"))
print(f"MemTotal (nodo)     : {gib(total)}   <- NO usar como limite si hay cgroup")
print(f"MemAvailable        : {gib(avail)}")
print(f"LIMITE EFECTIVO RAM : {gib(eff)}  (min(cgroup, MemTotal))")
print(f"cores               : {os.cpu_count()}")
print(f"disco libre CWD     : {gib(cwd_free)}")
print(f"disco libre HOME    : {gib(home_free)}")
print(f"GPU visible         : {gpu}")
# Referencias de dimensionado PUBLICAS (ver prep_boltz.md para URLs y alcance):
# descarga minima ~4.2 GB (boltz2_conf.ckpt 2.29 GB + mols.tar 1.86 GB) + torch/deps
# complejo objetivo: 352 residuos (316 receptor + 36 binder)
need_disk = 10 * 2**30
print(f"disco >=10 GiB CWD+HOME: {cwd_free>=need_disk and home_free>=need_disk}")
print("No se estima tiempo ni se promete ejecucion: medir antes de decidir.")
