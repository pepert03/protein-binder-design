# Ruta verificada Kaggle Boltz 2.2.1 (challenge3) — evidencia para verificador

## Causa raíz fallo pip (25 oct ~18:42 CET)
- boltz 2.2.1 requires_python >=3.10,<3.13: https://pypi.org/pypi/boltz/json
- En py3.13, pip selecciona silenciosamente boltz 2.0.3 (>=3.10 sin cota): https://pypi.org/pypi/boltz/2.0.3/json
- boltz 2.x pinea dm-tree==0.1.8; dm-tree 0.1.8 solo tiene wheels cp37-cp312 (sin cp313): https://pypi.org/pypi/dm-tree/0.1.8/json
- dm-tree 0.1.10 sí tiene cp313, pero boltz exige ==0.1.8: https://pypi.org/pypi/dm-tree/json
- Imagen Kaggle actual = runtime Colab, python 3.13 único, uv>=0.13.4 en PATH: https://github.com/Kaggle/docker-python/blob/main/Dockerfile.tmpl
- Detalle interno del build fallido NO probado (log efímero perdido, iba con -q). Causa presentada como respaldada por metadatos, no confirmada por log.

## Prueba CPU local (entorno aislado, 25 oct ~18:59-19:01 CET)
- uv 0.13.3; `uv python install 3.12` -> CPython 3.12.15 (<1s)
- `uv venv --python 3.12 /tmp/boltzenv`; `uv pip install torch(cpu) boltz==2.2.1` -> 71 paquetes, dm-tree==0.1.8 binario, sin compilar
- `import tree, torch, boltz` OK (torch 2.14.1+cpu); `boltz predict --help` rc=0
- dm-tree==0.1.8 wheel binario también verificado en py3.10 venv (/tmp/dtenv)

## Ejecución Kaggle (25oct 19:05-19:10 CET, GPU T4 x2)
- Celda 2 (install): `uv python install 3.12 && uv venv --python 3.12 /kaggle/working/boltzenv && uv pip install --python .../python --no-cache torch boltz==2.2.1` — ~90 s, 6.3GB usados en /kaggle/working (14.6GB libres)
- Celda 5 (puerta CUDA, PASADA): py 3.12.15 | boltz 2.2.1 | torch 2.14.1+cu130 cuda True ndev 2 | gpu 0/1 Tesla T4 14.6 GiB | driver 580.178.04
- Celda 3 (inputs): recA len 440; sha1 YAMLs: top02 58e5b22da917b798779237e5c8fdb4b972b7aca2, nativo a353048878351afb231c964c1fb921c6ded1f570, aleatorio36 f7f36d775a64bab4bf528e32336511ff25e2584f (idénticos a los verificados previamente por el verificador)
- Celda 4: /kaggle/working/boltzenv/bin/boltz predict por cada YAML, --recycling_steps 3 --sampling_steps 200 --diffusion_samples 1 --devices 1; log /kaggle/working/boltz_run.log; hashes out_hashes.json; empaquetado boltz_results.tar.gz
