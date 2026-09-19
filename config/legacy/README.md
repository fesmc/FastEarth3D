# Legacy build configuration

Superseded by **configme**, which is the supported way to configure a build:

```bash
configme install FastEarth3D -m macbook -c gfortran
```

See the Install section of the top-level [README](../../README.md).

## What is here

- `config.py` — the pre-configme Makefile generator. It splices a compiler
  fragment into the `<COMPILER_CONFIGURATION>` placeholder of the template
  `config/Makefile` and writes the top-level `Makefile`. configme now performs
  exactly this splice, using its own shipped machine/compiler fragments.
- `macbook_gfortran` — the only fragment written for `config.py`; nothing else
  reads it.

Note that `config/Makefile`, `config/common.mk` and
`config/Makefile_fastearth.mk` are **not** legacy — they are the live template
and rule set that configme itself consumes.

## Still usable

`config.py` resolves paths relative to the working directory, not to its own
location, so it still works unchanged when run from the repository root:

```bash
python config/legacy/config.py config/legacy/macbook_gfortran
```

It is kept as the record of the pre-configme build path. The production
optimization flags it carried (`-O3 -mcpu=native -funroll-loops -ffast-math`,
see `doc/performance-assessment.md` §1) now live in the project-tier configme
fragment `.configme/machines/macbook.mk`, which is what the build actually
uses.
