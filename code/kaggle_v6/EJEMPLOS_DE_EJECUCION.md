# Ejemplos de ejecución

Código de la Version6, 4 de octubre de 2026. La ejecución en Kaggle ya está activa. No lances otra copia mientras siga activa: gastarías cuota y duplicarías el mismo control. Los resultados de frontera del modelo, memoria y pose siguen pendientes.

## 1. Comprobar sintaxis localmente, sin modelo

Descomprime el ZIP y entra en su carpeta. Con Python 3 puedes comprobar los tres archivos sin instalar BoltzGen ni cargar pesos:

```bash
python3 -m py_compile kaggle_version6.py driver.py boundary_assertions.py
```

Si termina sin error, solo confirma sintaxis. No confirma que haya memoria suficiente ni que el modelo recupere la pose.

## 2. Ejecutar la regresión sintética de aserciones

Necesitas Python y NumPy instalado:

```bash
python3 boundary_assertions.py
```

Salida esperada:

```text
Synthetic assertion regression passed. Actual model boundary not executed.
```

Esta prueba usa arrays ficticios. Comprueba que las aserciones rechazan enlaces borrados y geometría resuelta inventada en residuos ausentes. No carga checkpoint, no usa GPU y no valida todavía la instancia real del modelo.

## 3. Inspeccionar el contenido de la celda sin ejecutarla

Este ejemplo usa solo la biblioteca estándar de Python. Lee el paquete incluido en la celda, comprueba el CRC y muestra sus archivos. No instala paquetes, no descarga pesos y no ejecuta el driver:

```bash
python3 - <<'PY'
from pathlib import Path
import ast, base64, io, zipfile
source = Path('kaggle_version6.py').read_text()
tree = ast.parse(source)
payload = next(
    node.args[0].value
    for node in ast.walk(tree)
    if isinstance(node, ast.Call)
    and isinstance(node.func, ast.Attribute)
    and node.func.attr == 'b64decode'
)
with zipfile.ZipFile(io.BytesIO(base64.b64decode(payload))) as archive:
    assert archive.testzip() is None
    print('\n'.join(archive.namelist()))
print('CRC correcto. Ningún modelo ejecutado.')
PY
```

## 4. Repetir el control completo en Kaggle, solo cuando termine la ejecución activa

Este es un cálculo de predicción, no una prueba local ligera. Necesita una sesión Kaggle con Python, Internet, GPU T4 y cuota disponible. Comprueba que el notebook sea privado. El driver usa una GPU aunque la opción del entorno se llame GPU T4 x2; no suma la memoria de las dos GPU.

Abre un notebook privado independiente y pega el contenido completo de `kaggle_version6.py` en una celda. Selecciona GPU T4 y Save & Run All. No ejecutes además Run All en la sesión Draft: crearías una segunda ejecución.

La celda instala BoltzGen 0.3.2, retira paquetes multimedia incompatibles de esa sesión, descarga moléculas y checkpoint públicos con revisión fijada, extrae las referencias y llama al driver. Por eso no debes ejecutarla en un entorno personal que quieras conservar sin cambios. Usa un entorno aislado de Kaggle.

Si prefieres una consola dentro de un entorno Kaggle aislado y ya has descomprimido el ZIP en la carpeta actual:

```bash
python3 kaggle_version6.py
```

No ejecutes `driver.py` directamente desde tu ordenador. Espera rutas absolutas de Kaggle y archivos auxiliares extraídos por la celda.

## Qué ejecuta el control completo

Primero prueba Fab6ARU completo, con 1058 tokens, y después EGF1NQL, con 677. Cada referencia tiene dos brazos explícitos de disulfuros, false y true, con un sample por brazo. Antes de cada predicción comprueba la instancia real del masker y las features que atraviesan esa frontera. Si esa comprobación falla, detiene el flujo. Si falla la predicción del primer brazo de una referencia, omite el segundo para no repetir la misma carga de recursos.

Un sample reduce concurrencia frente a cinco, pero no garantiza que Fab quepa y no sustituye al benchmark anterior de cinco. No se recortan las cadenas ni se cambian los filtros originales de candidatos.

## Dónde mirar el resultado

La celda conserva `CORRECTED_SINGLE_SAMPLE_CONTROL.zip` en `/kaggle/working`. Dentro deben revisarse `arm_results.json`, los archivos `actual_model_boundary.json` de cada brazo y los outputs de coordenadas y métricas. Un estado verde de Kaggle no basta: la celda conserva también fallos y el ZIP puede existir aunque un brazo haya fallado.

Para revisar el resumen tras descomprimir el ZIP de resultados:

```bash
python3 -m json.tool arm_results.json
```

`actual_model_boundary_pass` confirma la comprobación técnica de ese brazo. `prediction_complete` confirma que el driver encontró coordenadas y métricas válidas de su sample. Ninguna de esas dos etiquetas demuestra afinidad, unión a ratón o selectividad entre pH 6,5 y 7,4. La recuperación de pose requiere un análisis posterior contra la referencia.

## Archivos exactos y estado

Los YAML conservan PREPARED en el nombre porque ese es el nombre literal que espera el driver. En esta Version6 ya fueron usados para arrancar la ejecución. No renombres sus archivos sin cambiar también el driver.

El ZIP contiene el código verificado antes del arranque y no contiene resultados terminados. El estado observado a las 14:24 era Running.

https://www.kaggle.com/code/peperidruejo/notebook6198e9562d?scriptVersionId=355161133
