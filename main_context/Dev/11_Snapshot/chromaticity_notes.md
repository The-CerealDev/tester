# MAD-X `DQ1`/`DQ2` And Operational Chromaticity

## Sources

- CERN MAD-X user guide PDF: <https://madx.web.cern.ch/madx/releases/last-rel/madxuguide.pdf>
- CERN MAD-X primer PDF: <https://madx.web.cern.ch/madx/doc/madx_primer.pdf>

The MAD-X user guide defines the `SUMM` table chromaticity fields as:

- `DQ1`: horizontal chromaticity, `dq1 = dQ1 / dpt`
- `DQ2`: vertical chromaticity, `dq2 = dQ2 / dpt`

The same guide states that MAD-X uses `PT` as the longitudinal variable for
dispersive and chromatic functions, rather than using `DELTAP` directly. It
also gives the approximate relation:

```text
PT ~= beta_rel * DELTAP
```

where:

```text
DELTAP = delta_p / p0
```

and `beta_rel` is the relativistic Lorentz beta of the beam at the selected
cycle time.

The MAD-X primer confirms that `DQ1` and `DQ2` are the global chromaticities
used in TWISS summary output and matching constraints, but the derivative
variable convention is specified in the user guide.

## Conversion

MAD-X reports:

```text
DQ_madx = dQ / dPT
```

Operational measurements are usually quoted as:

```text
chromaticity = dQ / d(delta_p / p0)
```

Using the MAD-X relation `PT ~= beta_rel * DELTAP`:

```text
dPT ~= beta_rel * dDELTAP
```

so:

```text
dQ / dDELTAP ~= beta_rel * dQ / dPT
```

Therefore:

```text
chromaticity_operational ~= beta_rel * DQ_madx
```

This is the convention used by the Dev/11 snapshot summary columns:

- `dqx`: horizontal operational chromaticity, approximately `beta_rel * MAD-X DQ1`
- `dqy`: vertical operational chromaticity, approximately `beta_rel * MAD-X DQ2`
- `madx_dqx_dpt`: raw MAD-X horizontal `DQ1 = dQ1/dPT`
- `madx_dqy_dpt`: raw MAD-X vertical `DQ2 = dQ2/dPT`

## Measurement Comparison

If a measurement fits tune as a function of momentum error, the fitted slope is:

```text
dQ / d(delta_p / p0)
```

and should be compared directly with the snapshot `dqx`/`dqy` columns.

If the measurement instead fits momentum error as a function of tune, the slope is:

```text
d(delta_p / p0) / dQ
```

and must be inverted before comparing with chromaticity.

At ISIS injection, the relativistic beam beta is about `0.366`. A typical raw
MAD-X value near `DQ1 ~= -13.5` therefore corresponds to:

```text
0.366 * -13.5 ~= -4.9
```

This is consistent with measured operational chromaticities of order `-4` and
normalised chromaticities of order `-1` when divided by tune.
