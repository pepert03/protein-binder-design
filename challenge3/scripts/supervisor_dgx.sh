#!/bin/bash
# Supervisor local DGX: detiene g3b de forma SEGURA tras completar trayectoria 1.
# Condicion verificable: linea TRAJ traj=1 en log + g3b_results.json contiene "traj": 1
# + mtime del json estable 20 s + g3b_1.pdb no vacio. Copia artefactos a export/ ANTES de matar.
LOG="$HOME/glp1r/g3b_run.log"; JSON="$HOME/glp1r/g3b_results.json"; PDB="$HOME/glp1r/g3b_1.pdb"
while true; do
  if grep -q "'traj': 1" "$LOG" 2>/dev/null && [ -s "$JSON" ] && grep -q '"traj": 1' "$JSON" 2>/dev/null; then
    m1=$(stat -c %Y "$JSON"); sleep 20; m2=$(stat -c %Y "$JSON")
    if [ "$m1" = "$m2" ] && [ -s "$PDB" ]; then
      mkdir -p "$HOME/glp1r/export"
      cp -f "$PDB" "$JSON" "$LOG" "$HOME/glp1r/export/" 2>/dev/null
      PID=$(pgrep -f 'run_glp1r_colabdesign.py' | head -1)
      if [ -n "$PID" ]; then kill "$PID" 2>/dev/null; sleep 15; kill -9 "$PID" 2>/dev/null; fi
      { echo "supervisor: trayectoria 1 completa; artefactos exportados; proceso detenido tras verificacion";
        date -u; ls -la "$HOME/glp1r/export/"; } > "$HOME/glp1r/g3b_SUPERVISOR.txt" 2>&1
      break
    fi
  fi
  # condicion de seguridad: si el proceso murio solo, salir y dejar constancia
  if ! pgrep -f 'run_glp1r_colabdesign.py' > /dev/null; then
    echo "supervisor: proceso ya no existe a $(date -u); sin accion" > "$HOME/glp1r/g3b_SUPERVISOR.txt"; break
  fi
  sleep 20
done
