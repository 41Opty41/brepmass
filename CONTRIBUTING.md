# Contributing to brepmass

Contributions, bug reports, and numerical accuracy tests are welcome.

## Code Standards

- **Strict C99**: All core sources (`src/*.c`, `cli/*.c`, `tests/*.c`) must compile cleanly with `-std=c99 -Wall -Wextra -Werror` (and MSVC `/W4 /WX`) using only `libc` and `libm`.
- **English only** in source code, comments, headers, CLI messages, and commit messages.
- **Public ABI**: Export new public symbols in `include/brepmass.h` with `BM_API`.

## Building & Testing

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
pytest -v tests/test_python_binding.py
```
