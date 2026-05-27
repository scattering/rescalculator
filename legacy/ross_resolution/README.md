# Legacy Ross/RWE Resolution Reference

This directory preserves a legacy `res3` resolution implementation and related
multi-crystal resolution notes as reference material for `rescalculator`.

Provenance:

- Moved from `/Users/williamratcliff/tascopilot/ross resolution` on 2026-05-27.
- Imported at William Ratcliff's request so the material lives with the
  resolution calculator rather than in `tascopilot`.
- Original source comments and naming are preserved for attribution and
  convention tracking.

Use policy:

- Treat this as a legacy reference and numerical/convention cross-check.
- Do not make this the default runtime path.
- Prefer the maintained NumPy/Numba resolution implementation for production
  calculations.
- If a convention from this code is adopted, port the mathematical convention
  deliberately into the maintained implementation and add regression tests.
