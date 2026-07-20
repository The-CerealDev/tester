"""
errors.py

Standalone MAD-X error-table utilities for the ISIS RCS optics GUI backend.

This layer owns pure file/DataFrame handling for MAD-X error tables. MAD-X
execution remains in MadxModel.
"""

from datetime import datetime
from pathlib import Path
import re

import numpy as np
import pandas as pd


MAD_X_ERROR_COLUMNS = ("name", "dx", "dy", "ds", "dphi", "dtheta", "dpsi")
MAD_X_ERROR_NUMERIC_COLUMNS = ("dx", "dy", "ds", "dphi", "dtheta", "dpsi")
MAD_X_SETERR_NUMERIC_COLUMNS = (
    tuple(
        item
        for order in range(21)
        for item in (f"k{order}l", f"k{order}sl")
    )
    + MAD_X_ERROR_NUMERIC_COLUMNS
    + (
        "mrex",
        "mrey",
        "mredx",
        "mredy",
        "arex",
        "arey",
        "mscalx",
        "mscaly",
        "rfm_freq",
        "rfm_harmon",
        "rfm_lag",
    )
    + tuple(
        item
        for order in range(21)
        for item in (f"p{order}l", f"p{order}sl")
    )
)
MAD_X_SETERR_COLUMNS = ("name",) + MAD_X_SETERR_NUMERIC_COLUMNS
MAD_X_ERROR_COLUMN_ALIASES = {
    "NAME": "name",
    "Name": "name",
    "element": "name",
    "Element": "name",
    "ELEMENT": "name",
    "DX": "dx",
    "DY": "dy",
    "DS": "ds",
    "DPHI": "dphi",
    "DTHETA": "dtheta",
    "DPSI": "dpsi",
}
MAD_X_ERROR_FAMILY_ALIASES = {
    "DIPOLE": "DIP",
    "DIP": "DIP",
    "QD": "QD",
    "QF": "QF",
    "QC": "QDS",
    "QDS": "QDS",
}


def _normalise_column_name(column):
    stripped = str(column).strip()
    return MAD_X_ERROR_COLUMN_ALIASES.get(stripped, stripped.lower())


def _strip_quotes(value):
    text = str(value).strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in ("'", '"'):
        return text[1:-1]
    return text


def normalise_error_table_columns(df, fill_missing=True, copy=True):
    """
    Return a DataFrame with normalised lower-case MAD-X error-table columns.

    Standard MAD-X columns such as NAME, DX and DTHETA are mapped to lower-case
    pandas columns. Extra columns are preserved after the standard columns.
    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    out = df.copy() if copy else df

    new_columns = [_normalise_column_name(column) for column in out.columns]
    if len(new_columns) != len(set(new_columns)):
        duplicates = sorted(
            {column for column in new_columns if new_columns.count(column) > 1}
        )
        raise ValueError(
            "Duplicate columns produced by normalisation: "
            f"{duplicates}"
        )

    out.columns = new_columns

    if "name" in out.columns:
        out["name"] = out["name"].map(_strip_quotes)

    if fill_missing:
        for column in MAD_X_ERROR_NUMERIC_COLUMNS:
            if column not in out.columns:
                out[column] = 0.0

    for column in MAD_X_ERROR_NUMERIC_COLUMNS:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype(float)

    ordered_columns = [
        column for column in MAD_X_ERROR_COLUMNS if column in out.columns
    ]
    extra_columns = [
        column for column in out.columns if column not in ordered_columns
    ]

    return out[ordered_columns + extra_columns]


def validate_error_table(
    df,
    required_columns=None,
    allow_extra_columns=True,
    require_non_empty=True,
    require_unique_names=False,
):
    """
    Validate and return a normalised MAD-X error-table DataFrame.

    The input DataFrame is normalised internally before validation. By default
    the table must contain name plus all standard MAD-X alignment columns.
    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")

    out = normalise_error_table_columns(df, fill_missing=False, copy=True)

    if required_columns is None:
        required_columns = MAD_X_ERROR_COLUMNS
    required_columns = tuple(_normalise_column_name(column) for column in required_columns)

    missing = [column for column in required_columns if column not in out.columns]
    if missing:
        raise ValueError(f"Missing required error-table columns: {missing}")

    if require_non_empty and out.empty:
        raise ValueError("Error table is empty.")

    if "name" in out.columns:
        if out["name"].isna().any():
            raise ValueError("Error table contains null or empty names.")

        names = out["name"].astype(str).str.strip()
        if (names == "").any() or (names.str.lower().isin({"nan", "none"})).any():
            raise ValueError("Error table contains null or empty names.")
        out["name"] = names

        if require_unique_names and out["name"].duplicated().any():
            duplicates = out.loc[out["name"].duplicated(), "name"].tolist()
            raise ValueError(f"Error table contains duplicate names: {duplicates}")

    for column in MAD_X_ERROR_NUMERIC_COLUMNS:
        if column in out.columns:
            out[column] = pd.to_numeric(out[column], errors="coerce").astype(float)
            if not np.isfinite(out[column].to_numpy()).all():
                raise ValueError(f"Column {column!r} contains non-finite values.")

    if not allow_extra_columns:
        extra_columns = [column for column in out.columns if column not in MAD_X_ERROR_COLUMNS]
        if extra_columns:
            raise ValueError(f"Unexpected error-table columns: {extra_columns}")

    return out


def _split_tfs_row(line):
    return re.findall(r'"[^"]*"|\'[^\']*\'|\S+', line)


def read_error_table(path, validate=True):
    """
    Read a MAD-X/TFS-style error table into a normalised pandas DataFrame.

    The parser supports the files written by write_error_table() and broader
    MAD-X EFIELD tables that contain the standard alignment columns.
    """

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing error table: {path}")

    lines = path.read_text().splitlines()
    header = None
    data_start = None

    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("*"):
            header = stripped.lstrip("*").split()
        elif stripped.startswith("$") and header is not None:
            data_start = index + 1
            break

    if header is None or data_start is None:
        raise ValueError(f"Could not find MAD-X table header/type lines in {path}")

    rows = []
    for line in lines[data_start:]:
        stripped = line.strip()
        if (
            not stripped
            or stripped.startswith("@")
            or stripped.startswith("!")
            or stripped.startswith("*")
            or stripped.startswith("$")
        ):
            continue

        parts = _split_tfs_row(stripped)
        if len(parts) != len(header):
            raise ValueError(
                f"Could not parse row with {len(parts)} fields; "
                f"expected {len(header)} in {path}: {stripped}"
            )
        rows.append([_strip_quotes(part) for part in parts])

    df = pd.DataFrame(rows, columns=header)
    df = normalise_error_table_columns(df, fill_missing=True, copy=False)

    if validate:
        df = validate_error_table(df)

    return df


def write_error_table(
    df,
    path,
    table_name="ERROR_TABLE",
    validate=True,
    include_metadata=True,
):
    """
    Write a MAD-X-readable error table suitable for readtable/seterr.

    The written file uses the fuller MAD-X EFIELD column set seen in the ISIS
    reference tables, because SETERR expects columns such as K0L as well as the
    alignment columns. Missing non-alignment fields are written as zero. Reading
    the file back through read_error_table() returns lower-case pandas columns.
    """

    if validate:
        out = normalise_error_table_columns(df, fill_missing=True, copy=True)
        out = validate_error_table(out)
    else:
        out = normalise_error_table_columns(df, fill_missing=True, copy=True)

    write_columns = {"name": out["name"].astype(str)}
    for column in MAD_X_SETERR_NUMERIC_COLUMNS:
        if column in out.columns:
            write_columns[column] = pd.to_numeric(
                out[column],
                errors="coerce",
            ).astype(float)
        else:
            write_columns[column] = 0.0

    write_out = pd.DataFrame(
        write_columns,
        columns=MAD_X_SETERR_COLUMNS,
    )

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with path.open("w") as file_handle:
        if include_metadata:
            file_handle.write(f'@ NAME %08s "{table_name}"\n')
            file_handle.write('@ TYPE %08s "EFIELD"\n')
            file_handle.write('@ ORIGIN %08s "errors.py"\n')

        upper_columns = [column.upper() for column in MAD_X_SETERR_COLUMNS]
        file_handle.write("* " + " ".join(upper_columns) + "\n")
        file_handle.write("$ %s " + " ".join(["%le"] * len(MAD_X_SETERR_NUMERIC_COLUMNS)) + "\n")

        for _, row in write_out.iterrows():
            values = [f'"{row["name"]}"']
            values.extend(
                f"{float(row[column]):.16e}"
                for column in MAD_X_SETERR_NUMERIC_COLUMNS
            )
            file_handle.write(" ".join(values) + "\n")

    return str(path)


def _load_error_table_like(table, fill_missing=True):
    if isinstance(table, (str, Path)):
        return read_error_table(table)
    if isinstance(table, pd.DataFrame):
        return normalise_error_table_columns(
            table,
            fill_missing=fill_missing,
            copy=True,
        )
    raise TypeError("tables must contain pandas DataFrames or file paths.")


def combine_error_tables(tables, mode="sum", fill_missing=True, validate=True):
    """
    Combine multiple error tables into one DataFrame.

    Only mode="sum" is currently supported. Element order follows first
    appearance across the input tables.
    """

    tables = list(tables)
    if not tables:
        raise ValueError("At least one error table is required.")

    if mode != "sum":
        raise ValueError(f"Unsupported combine mode: {mode!r}")

    order = []
    combined = {}

    for table in tables:
        df = _load_error_table_like(table, fill_missing=fill_missing)
        if validate:
            df = validate_error_table(df, require_unique_names=False)

        for _, row in df.iterrows():
            name = str(row["name"])
            if name not in combined:
                order.append(name)
                combined[name] = {column: 0.0 for column in MAD_X_ERROR_NUMERIC_COLUMNS}

            for column in MAD_X_ERROR_NUMERIC_COLUMNS:
                value = row[column] if column in row.index else 0.0
                if pd.isna(value):
                    value = 0.0
                combined[name][column] += float(value)

    rows = [
        {"name": name, **combined[name]}
        for name in order
    ]
    out = pd.DataFrame(rows, columns=MAD_X_ERROR_COLUMNS)

    if validate:
        out = validate_error_table(out)

    return out


def filter_error_table(
    df,
    names=None,
    patterns=None,
    element_types=None,
    case=False,
    regex=True,
    invert=False,
):
    """
    Filter an error table by explicit element names, name patterns, or type.
    """

    out = normalise_error_table_columns(df, fill_missing=True, copy=True)
    mask = pd.Series(False, index=out.index)
    any_filter = False

    name_series = out["name"].astype(str)
    comparable_names = name_series if case else name_series.str.lower()

    if names is not None:
        any_filter = True
        names = {str(name) for name in names}
        if not case:
            names = {name.lower() for name in names}
        mask |= comparable_names.isin(names)

    if patterns is not None:
        any_filter = True
        for pattern in patterns:
            pattern = str(pattern)
            if regex:
                mask |= name_series.str.contains(pattern, case=case, regex=True, na=False)
            else:
                comparable_pattern = pattern if case else pattern.lower()
                mask |= comparable_names.str.contains(
                    re.escape(comparable_pattern),
                    regex=True,
                    na=False,
                )

    if element_types is not None:
        any_filter = True
        type_column = None
        for candidate in ("element_type", "keyword"):
            if candidate in out.columns:
                type_column = candidate
                break

        if type_column is None:
            raise ValueError(
                "element_types requires an 'element_type' or 'keyword' column."
            )

        types = {str(element_type) for element_type in element_types}
        type_series = out[type_column].astype(str)
        if not case:
            types = {element_type.lower() for element_type in types}
            type_series = type_series.str.lower()
        mask |= type_series.isin(types)

    if not any_filter:
        return out

    if invert:
        mask = ~mask

    return out.loc[mask].reset_index(drop=True)


def flip_error_columns(df, columns, copy=True):
    """
    Flip signs of selected MAD-X error columns.
    """

    out = normalise_error_table_columns(df, fill_missing=True, copy=copy)
    requested_columns = [_normalise_column_name(column) for column in columns]

    invalid = [
        column for column in requested_columns
        if column not in MAD_X_ERROR_NUMERIC_COLUMNS
    ]
    if invalid:
        raise ValueError(
            "Can only flip numeric MAD-X error columns: "
            f"{invalid}"
        )

    for column in requested_columns:
        out[column] = -1.0 * out[column].astype(float)

    return out


def _normalise_error_family(family):
    key = str(family).strip().upper()
    if key not in MAD_X_ERROR_FAMILY_ALIASES:
        raise ValueError(
            "Unsupported magnet family "
            f"{family!r}. Supported families: {sorted(MAD_X_ERROR_FAMILY_ALIASES)}"
        )
    return MAD_X_ERROR_FAMILY_ALIASES[key]


def _normalise_error_names(names, case=False):
    values = {str(name).strip().strip("\"'") for name in names}
    if case:
        return values
    return {name.lower() for name in values}


def _error_table_zero_mask(
    df,
    names=None,
    families=None,
    patterns=None,
    case=False,
    regex=True,
):
    out = normalise_error_table_columns(df, fill_missing=True, copy=True)
    name_series = out["name"].astype(str)
    comparable_names = name_series if case else name_series.str.lower()
    mask = pd.Series(False, index=out.index)
    any_selector = False

    if names is not None:
        any_selector = True
        requested_names = _normalise_error_names(names, case=case)
        mask |= comparable_names.isin(requested_names)

    if families is not None:
        any_selector = True
        requested_families = {_normalise_error_family(family) for family in families}
        upper_names = name_series.str.upper()
        family_mask = pd.Series(False, index=out.index)
        for family in requested_families:
            family_mask |= upper_names.str.endswith(f"_{family}")
        mask |= family_mask

    if patterns is not None:
        any_selector = True
        for pattern in patterns:
            pattern = str(pattern)
            if regex:
                mask |= name_series.str.contains(pattern, case=case, regex=True, na=False)
            else:
                comparable_pattern = pattern if case else pattern.lower()
                mask |= comparable_names.str.contains(
                    re.escape(comparable_pattern),
                    regex=True,
                    na=False,
                )

    if not any_selector:
        raise ValueError("At least one of names, families or patterns is required.")

    return out, mask


def zero_error_table_magnets(
    df,
    names=None,
    families=None,
    patterns=None,
    columns=None,
    case=False,
    regex=True,
    require_match=True,
    return_zeroed_names=False,
):
    """
    Return an error table with selected magnet misalignment columns set to zero.

    Family names follow the GUI/survey language. In particular, ``QC`` maps to
    the MAD-X table rows named ``SPn_QDS``.
    """

    out, mask = _error_table_zero_mask(
        df,
        names=names,
        families=families,
        patterns=patterns,
        case=case,
        regex=regex,
    )
    if require_match and not mask.any():
        raise ValueError("No error-table rows matched the requested magnets.")

    if columns is None:
        zero_columns = MAD_X_ERROR_NUMERIC_COLUMNS
    else:
        zero_columns = tuple(_normalise_column_name(column) for column in columns)
        invalid = [column for column in zero_columns if column not in MAD_X_ERROR_NUMERIC_COLUMNS]
        if invalid:
            raise ValueError(
                "Can only zero MAD-X misalignment columns: "
                f"{invalid}"
            )

    zeroed_names = out.loc[mask, "name"].astype(str).tolist()
    for column in zero_columns:
        out.loc[mask, column] = 0.0

    if return_zeroed_names:
        return out, zeroed_names
    return out


def write_zeroed_error_table_copy(
    source_path,
    output_dir=None,
    names=None,
    families=None,
    patterns=None,
    columns=None,
    timestamp=None,
    suffix="zeroed",
    table_name=None,
    case=False,
    regex=True,
    require_match=True,
    return_table=False,
):
    """
    Write a timestamped copy of an error table with selected magnets zeroed.

    Returns a dictionary containing the written path and the row names that were
    zeroed. Set ``return_table=True`` to include the edited DataFrame.
    """

    source_path = Path(source_path)
    if not source_path.is_file():
        raise FileNotFoundError(f"Missing source error table: {source_path}")

    original = read_error_table(source_path)
    edited, zeroed_names = zero_error_table_magnets(
        original,
        names=names,
        families=families,
        patterns=patterns,
        columns=columns,
        case=case,
        regex=regex,
        require_match=require_match,
        return_zeroed_names=True,
    )

    timestamp = timestamp or datetime.now().strftime("%Y%m%d_%H%M%S")
    suffix = str(suffix).strip() or "zeroed"
    output_dir = source_path.parent if output_dir is None else Path(output_dir)
    output_path = output_dir / f"{source_path.stem}_{suffix}_{timestamp}.tfs"
    written_path = write_error_table(
        edited,
        output_path,
        table_name=table_name or f"{source_path.stem}_{suffix}",
    )

    result = {
        "path": written_path,
        "zeroed_names": zeroed_names,
        "n_zeroed": len(zeroed_names),
    }
    if return_table:
        result["table"] = edited
    return result


def map_error_names_to_lattice(
    df,
    twiss_df,
    name_col="name",
    lattice_name_col="name",
    case=False,
):
    """
    Check error-table names against a lattice/TWISS DataFrame.

    Matching is exact after optional case normalisation. Unmatched names are
    retained.
    """

    out = normalise_error_table_columns(df, fill_missing=True, copy=True)

    if not isinstance(twiss_df, pd.DataFrame):
        raise TypeError("twiss_df must be a pandas DataFrame.")

    name_col = _normalise_column_name(name_col)
    if name_col not in out.columns:
        raise ValueError(f"Error table is missing name column {name_col!r}.")

    if lattice_name_col not in twiss_df.columns:
        raise ValueError(
            f"TWISS DataFrame is missing lattice name column {lattice_name_col!r}."
        )

    lattice_names = twiss_df[lattice_name_col].astype(str)
    lookup = {}
    for lattice_name in lattice_names:
        key = lattice_name if case else lattice_name.lower()
        lookup.setdefault(key, lattice_name)

    matched_names = []
    in_lattice = []
    for name in out[name_col].astype(str):
        key = name if case else name.lower()
        matched_name = lookup.get(key)
        matched_names.append(matched_name)
        in_lattice.append(matched_name is not None)

    mapped = out.copy()
    mapped["matched_name"] = matched_names
    mapped["in_lattice"] = in_lattice

    leading_columns = ["name", "matched_name", "in_lattice"]
    remaining_columns = [
        column for column in mapped.columns if column not in leading_columns
    ]

    return mapped[leading_columns + remaining_columns]


def summarise_error_table(df):
    """
    Return a compact summary dictionary for an error table.
    """

    out = validate_error_table(df, require_non_empty=False)

    if out.empty:
        nonzero_rows = 0
    else:
        numeric = out[list(MAD_X_ERROR_NUMERIC_COLUMNS)].abs()
        nonzero_rows = int((numeric > 0.0).any(axis=1).sum())

    summary = {
        "n_rows": int(len(out)),
        "n_unique_names": int(out["name"].nunique()) if "name" in out.columns else 0,
        "n_nonzero_rows": nonzero_rows,
        "columns": list(out.columns),
    }

    for column in MAD_X_ERROR_NUMERIC_COLUMNS:
        if column in out.columns and not out.empty:
            summary[f"max_abs_{column}"] = float(out[column].abs().max())
        else:
            summary[f"max_abs_{column}"] = 0.0

    return summary


def apply_error_table(madx_model, path, table_name="error_table"):
    """
    Thin wrapper around MadxModel.apply_error_table(...).

    MAD-X execution remains owned by MadxModel.
    """

    return madx_model.apply_error_table(
        error_table_path=path,
        table_name=table_name,
    )
