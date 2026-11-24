<div align="center">

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="../../assets/logo-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="../../assets/logo-light.svg">
  <img alt="brepmass logo" src="../../assets/logo-light.svg" width="380">
</picture>

**Лёгкая C99-библиотека без зависимостей, CLI и Python `ctypes`-обёртка для расчёта масс-инерционных характеристик B-rep/сеток, тесселяции полигонов и NURBS и работы с STL/OBJ.**

[English](../../README.md) · **Русский** · [中文](README.zh-CN.md) · [Deutsch](README.de.md) · [日本語](README.ja.md) · [Español](README.es.md)

[![C99](https://img.shields.io/badge/standard-C99%20(-Wall%20-Werror)-EA580C?style=flat-square)](../../CMakeLists.txt)
[![Zero Dependencies](https://img.shields.io/badge/dependencies-libc%20%2B%20libm-F59E0B?style=flat-square)](../../CMakeLists.txt)
[![License: MIT](https://img.shields.io/badge/license-MIT-16A34A?style=flat-square)](../../LICENSE)
[![Support](https://img.shields.io/badge/support-BTC%20%7C%20YooMoney-F59E0B?style=flat-square)](../../DONATE.md)

</div>

---

## Архитектура

```mermaid
flowchart LR
    P["Примитивы / NURBS<br/>2D-полигоны с отверстиями"] --> T["Тесселятор<br/>(tessellate.c)"]
    STL["Бинарный STL + пространственная<br/>сшивка вершин"] --> M["Вычислитель масс<br/>(mass.c)"]
    T --> M
    M --> OUT["Объём, площадь, ЦМ<br/>Тензоры инерции, оси Якоби<br/>Замкнутость и Эйлер"]
    T --> EXP["STL (Bin/ASCII) / OBJ<br/>CLI и Python ctypes"]
```

## Возможности

- **Чистый C99 (`-std=c99 -Wall -Wextra -Werror`)** — требуется только стандартная библиотека C (`libc` и `libm`).
- **Точные интегралы массы и инерции** — тетраэдральное потоковое интегрирование по теореме Гаусса–Остроградского для замкнутых тел, поверхностные интегралы (`dim=2`) для открытых оболочек, перенос тензора инерции к центру масс (теорема Штейнера) и циклический метод Якоби для главных моментов и осей.
- **Топологический контроль** — проверка парности направленных полурёбер (`is_closed`) и характеристика Эйлера `V - E + F`.
- **Аналитические эталоны** — `bm_analytic_box`, `bm_analytic_cylinder`, `bm_analytic_cone_frustum`, `bm_analytic_sphere`, `bm_analytic_torus`.
- **Водонепроницаемые тесселяторы** — базовые примитивы, триангуляция плоских 2D-полигонов с произвольными отверстиями методом ear-clipping (`bm_tessellate_polygon_2d`) и сэмплирование рациональных/нерациональных NURBS-поверхностей (`bm_tessellate_nurbs_patch`).
- **Файловый ввод-вывод и CLI** — бинарный/ASCII STL (со spatial-hash сшивкой вершин при чтении), Wavefront OBJ, утилита `brepmass` и [`python/brepmass.py`](../../python/brepmass.py).

## Быстрый старт

### Сборка через CMake

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

### CLI и Python

```bash
brepmass make cylinder 0.5 2.0 -o cylinder.stl
brepmass info cylinder.stl --density 7850
```

```python
from python.brepmass import tessellate_box, mesh_mass_props

verts, tris = tessellate_box(2.0, 3.0, 4.0)
props = mesh_mass_props(verts, tris, density=2.0)
print(props.volume, props.centroid, props.principal_moments)
```

## Поддержать проект

Если `brepmass` экономит вам время, проект можно добровольно поддержать ([подробнее](../../DONATE.md)):
- **Bitcoin (BTC)**: `bc1qqptjmpr684xglyqsp743ruufnwpsjv4dlzk9tx`
- **ЮMoney**: [yoomoney.ru/to/4100118415983435](https://yoomoney.ru/to/4100118415983435) (кошелёк `4100118415983435`)

## Документация и лицензия

- [Математические формулы и выводы (на русском)](MATH.ru.md) · [English](../MATH.md)
- [Участие в разработке](../../CONTRIBUTING.md) · [История изменений](../../CHANGELOG.md)
- Лицензия MIT — см. [`LICENSE`](../../LICENSE).
