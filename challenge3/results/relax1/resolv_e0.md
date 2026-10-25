# Resolucion anomalia E0 (relax2 vs relax3 sobre misma entrada)

DIAGNOSTICO: dos preparaciones pdbfixer INDEPENDIENTES de cand2_complex.pdb dan salidas
byte-distintas (hashes distintos). Comparacion atomo a atomo por clave (cadena,resseq,nombre,resname):
- Atomos PESADOS: 0 diferencias de coordenadas (2876 comunes, identicos).
- Hidrogenos: 2837/2837 con coords distintas, repartidos por los 352 residuos.
=> pdbfixer.addMissingHydrogens(7.0) coloca H de forma NO determinista entre procesos.
=> E0 difiere por NonbondedForce de los H. La fuerza de restriccion de A' contribuye 0 en el
punto inicial (verificado: descomposicion por grupo con/sin restriccion identica, CustomExternal=0).

IMPACTO EN METRICAS YA REPORTADAS: ninguna en choques/anclajes/RMSD — se miden solo con atomos
pesados, identicos entre preparaciones. SOLO las energias (nunca usadas en conclusiones) varian.

CORRECCION DE PROTOCOLO (aplicada): preparar UNA vez por sistema, guardar PDB preparado con hash,
y reutilizar la MISMA preparada en toda comparacion. No asumir estocasticidad: verificado por pares.

CORRECCION DE LECTURA (tu nota): contactos C3 601->577 DISMINUYEN (no "conservados").
Cifras tal cual: relax3 C3 21->3 choques, dmin 2.60, contactos -24, binder RMSD 0.92, receptor 0.000000.

ENTRADAS PREPARADAS CANONICAS (usar estas en cualquier reproduccion):
1653e354eb70d5f1b250d673279ce308f0b42c1f97a36642c4215af9544b76c0  prep_LIBRE6_2.pdb
28240cdd39643090b0e9b01e476f93e2528bcee7878509ffb711bf048510a41c  prep_nativo.pdb
28aedc5f6eaf55ae8e9866af9fe8ad4889a999c5956b6f97ede81ecdc36a2ee3  prep_ala6.pdb
2ea1709867528bfa50c443c0f6d3b5dee6fe55f80015ca5db23e9c0fb48561cc  relax3.py
9220be0019d56c36f793b3192ffca1a4d82bf246564e619f2a973bde3cde2158  audita_relax.py
