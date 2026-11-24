<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../assets/logo-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="../../assets/logo-light.svg">
  <img alt="brepmass logo" src="../../assets/logo-light.svg" width="380">
</picture>

**Biblioteca C99 sin dependencias, CLI y wrapper `ctypes` en Python para propiedades másicas, tensores de inercia, teselación NURBS/polígonos y STL/OBJ.**

[English](../../README.md) · [Русский](README.ru.md) · [中文](README.zh-CN.md) · [Deutsch](README.de.md) · [日本語](README.ja.md) · **Español**

[![C99](https://img.shields.io/badge/standard-C99%20(-Wall%20-Werror)-EA580C?style=flat-square)](../../CMakeLists.txt)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-libc%20%2B%20libm-F59E0B?style=flat-square)](../../CMakeLists.txt)
[![License: MIT](https://img.shields.io/badge/license-MIT-16A34A?style=flat-square)](../../LICENSE)
[![Support](https://img.shields.io/badge/support-BTC%20%7C%20YooMoney-F59E0B?style=flat-square)](../../DONATE.md)

</div>

---

## Arquitectura

```mermaid
flowchart LR
    P["Primitivas / NURBS<br/>Polígonos 2D con agujeros"] --> T["Teselador<br/>(tessellate.c)"]
    STL["STL binario + soldadura<br/>espacial de vértices"] --> M["Motor de masa<br/>(mass.c)"]
    T --> M
    M --> OUT["Volumen, área, centroide<br/>Tensores y ejes principales<br/>Cierre y Euler"]
    T --> EXP["STL (Bin/ASCII) / OBJ<br/>CLI y Python ctypes"]
```

## Características

- **C99 puro (`-std=c99 -Wall -Wextra -Werror`)** — solo requiere `libc` y `libm`.
- **Integrales exactas de masa e inercia** — flujo tetraédrico por teorema de la divergencia para sólidos cerrados, integrales de superficie (`dim=2`) para láminas abiertas, teorema de Steiner y diagonales de Jacobi para momentos y ejes principales.
- **Comprobación topológica** — emparejamiento de semiaristas orientadas (`is_closed`) y característica de Euler `V - E + F`.
- **Referencias analíticas** — caja, cilindro, tronco de cono, esfera y toro (`bm_analytic_*`).
- **Teselación estanca y E/S** — primitivas, triangulación 2D ear-clipping con agujeros, superficies NURBS racionales/no racionales, STL binario/ASCII (con soldadura por hash espacial) y Wavefront OBJ.

## Inicio rápido

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

```bash
brepmass make cylinder 0.5 2.0 -o cylinder.stl
brepmass info cylinder.stl --density 7850
```

## Apoyar el proyecto

Si `brepmass` te resulta útil, puedes apoyar el desarrollo con una donación voluntaria ([detalles](../../DONATE.md)):
- **Bitcoin (BTC)**: `bc1qqptjmpr684xglyqsp743ruufnwpsjv4dlzk9tx`
- **YooMoney**: [yoomoney.ru/to/4100118415983435](https://yoomoney.ru/to/4100118415983435) (`4100118415983435`)

## Documentación y licencia

- [Derivaciones matemáticas (EN)](../MATH.md) · [Contribuir](../../CONTRIBUTING.md) · [Historial de cambios](../../CHANGELOG.md)
- Licencia MIT — consulta [`LICENSE`](../../LICENSE).
