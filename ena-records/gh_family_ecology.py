#!/usr/bin/env python3
# ------------------------------------------------------------------------------
# Script: gh_family_ecology.py
# Author: Δημήτρης
# Date: 2025-08-19
#
# Description:
#   Summarize ecological patterns for GH chitinase families from ENA TSVs.
#   For each ghXX.tsv it computes:
#     - total records
#     - distinct species
#     - environmental vs non-environmental (%)
#     - host info presence (%), top hosts
#     - unique countries, top countries
#     - top isolation sources
#     - altitude stats (n, min, max, mean)
#     - collection year range (min, max)
#   It prints a readable per-family report and also writes CSV summaries.
#
# Input:
#   TSV files named gh*.tsv in the current directory (or run with --folder).
#
# Output:
#   - ecology_summary.csv (one row per family with key metrics)
#   - top_hosts.csv (family, host, count)
#   - top_sources.csv (family, isolation_source, count)
#   - top_countries.csv (family, country_norm, count)
#
# Usage:
#   python3 gh_family_ecology.py
#   python3 gh_family_ecology.py --folder ena-records
# ------------------------------------------------------------------------------
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

def normalize_country(s: str) -> str:
    """Heuristic country normalization: take the first token before ';' or ',' or ':'."""
    if not isinstance(s, str):
        return np.nan
    s = s.strip()
    for sep in [';', ',', ':', '|', '/']:
        if sep in s:
            s = s.split(sep)[0]
    return s.strip()

def coerce_bool_series(x: pd.Series) -> pd.Series:
    """Coerce mixed boolean-like series to True/False."""
    return x.astype(str).str.strip().str.lower().isin(['true', '1', 'yes'])

def safe_numeric(x: pd.Series) -> pd.Series:
    """Coerce to numeric; invalid values become NaN."""
    return pd.to_numeric(x, errors='coerce')

def extract_years(df: pd.DataFrame) -> pd.Series:
    """Extract year from collection_date, collection_date_start, collection_date_end."""
    cand_cols = [c for c in ['collection_date', 'collection_date_start', 'collection_date_end'] if c in df.columns]
    if not cand_cols:
        return pd.Series(dtype='Int64')
    years = pd.Series(dtype='Int64')
    for c in cand_cols:
        y = pd.to_datetime(df[c], errors='coerce', utc=True).dt.year
        years = y.fillna(years)
    return years.astype('Int64')

def top_counts(series: pd.Series, k: int = 10) -> pd.DataFrame:
    """Return top-k value counts as DataFrame with columns [value, count]."""
    if series is None:
        return pd.DataFrame(columns=['value', 'count'])
    s = series.dropna().astype(str).str.strip()
    vc = s.value_counts().head(k)
    return vc.rename_axis('value').reset_index(name='count')

def summarize_family(tsv_path: Path) -> dict:
    fam = tsv_path.stem.upper()  # e.g. GH18
    df = pd.read_csv(tsv_path, sep="\t")

    out = {'family': fam}
    n = len(df)
    out['total_records'] = n

    # species
    if 'scientific_name' in df.columns:
        out['distinct_species'] = int(df['scientific_name'].nunique())
    else:
        out['distinct_species'] = np.nan

    # environmental
    env_pct = np.nan
    env_n = np.nan
    if 'environmental_sample' in df.columns:
        env = coerce_bool_series(df['environmental_sample'])
        env_n = int(env.sum())
        env_pct = float(env_n / n * 100.0) if n > 0 else np.nan
    out['environmental_n'] = env_n
    out['environmental_pct'] = env_pct

    # host
    host_n = np.nan
    if 'host' in df.columns:
        host_n = int(df['host'].notna().sum())
    out['with_host_n'] = host_n
    out['with_host_pct'] = float(host_n / n * 100.0) if n and not np.isnan(host_n) else np.nan

    # countries
    unique_countries = np.nan
    country_norm_series = None
    if 'countrylocation' in df.columns:
        country_norm_series = df['countrylocation'].map(normalize_country)
        unique_countries = int(country_norm_series.dropna().nunique())
    out['unique_countries'] = unique_countries

    # isolation sources (raw strings)
    # We don't try to ontologize; just surface the most frequent labels.
    # altitude
    alt_n = alt_min = alt_max = alt_mean = np.nan
    if 'altitude' in df.columns:
        alt = safe_numeric(df['altitude'])
        alt_n = int(alt.notna().sum())
        if alt_n > 0:
            alt_min = float(np.nanmin(alt))
            alt_max = float(np.nanmax(alt))
            alt_mean = float(np.nanmean(alt))
    out['altitude_n'] = alt_n
    out['altitude_min'] = alt_min
    out['altitude_max'] = alt_max
    out['altitude_mean'] = alt_mean

    # years
    years = extract_years(df)
    out['year_min'] = int(years.min()) if years.notna().any() else np.nan
    out['year_max'] = int(years.max()) if years.notna().any() else np.nan

    # top tables
    top_hosts_df = top_counts(df['host'] if 'host' in df.columns else None)
    top_sources_df = top_counts(df['isolation_source'] if 'isolation_source' in df.columns else None)
    top_countries_df = top_counts(country_norm_series if country_norm_series is not None else None)

    # print readable report
    print(f"\n=== Family {fam} ===")
    print(f"Total records: {n}")
    print(f"Distinct species: {out['distinct_species']}")
    if not np.isnan(out['environmental_n']):
        print(f"Environmental: {out['environmental_n']} ({out['environmental_pct']:.1f}%)")
        print(f"Non-environmental: {n - out['environmental_n']} ({(100.0 - out['environmental_pct']):.1f}%)")
    if not np.isnan(out['with_host_n']):
        print(f"With host info: {out['with_host_n']} ({out['with_host_pct']:.1f}%)")
    if not np.isnan(out['unique_countries']):
        print(f"Unique countries: {out['unique_countries']}")
    if alt_n and not np.isnan(alt_n) and alt_n > 0:
        print(f"Altitude (n={alt_n}): min={alt_min:.1f}, max={alt_max:.1f}, mean={alt_mean:.1f}")
    if years.notna().any():
        print(f"Collection years: {out['year_min']}–{out['year_max']}")

    if not top_sources_df.empty:
        print("\nTop isolation sources:")
        for _, r in top_sources_df.iterrows():
            print(f"  {r['value']}: {r['count']}")
    if not top_hosts_df.empty:
        print("\nTop hosts:")
        for _, r in top_hosts_df.iterrows():
            print(f"  {r['value']}: {r['count']}")
    if not top_countries_df.empty:
        print("\nTop countries:")
        for _, r in top_countries_df.iterrows():
            print(f"  {r['value']}: {r['count']}")

    # attach family to top tables for CSV export
    for tdf in (top_hosts_df, top_sources_df, top_countries_df):
        if not tdf.empty:
            tdf.insert(0, 'family', fam)

    return out, top_hosts_df, top_sources_df, top_countries_df

def main():
    ap = argparse.ArgumentParser(description="Ecological summary for GH chitinase families from ENA TSVs.")
    ap.add_argument('--folder', type=str, default='.', help="Folder containing gh*.tsv (default: current directory)")
    ap.add_argument('--topk', type=int, default=10, help="Top-K items for hosts/sources/countries (default 10)")
    args = ap.parse_args()

    folder = Path(args.folder)
    tsvs = sorted(folder.glob('gh*.tsv'))
    if not tsvs:
        print(f"No TSV files found in {folder.resolve()} matching gh*.tsv")
        return

    summary_rows = []
    all_hosts = []
    all_sources = []
    all_countries = []

    for tsv in tsvs:
        out, top_hosts, top_sources, top_countries = summarize_family(tsv)
        summary_rows.append(out)
        if not top_hosts.empty:
            all_hosts.append(top_hosts)
        if not top_sources.empty:
            all_sources.append(top_sources)
        if not top_countries.empty:
            all_countries.append(top_countries)

    # Write CSV summaries
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv('ecology_summary.csv', index=False)

    if all_hosts:
        pd.concat(all_hosts, ignore_index=True).to_csv('top_hosts.csv', index=False)
    if all_sources:
        pd.concat(all_sources, ignore_index=True).to_csv('top_sources.csv', index=False)
    if all_countries:
        pd.concat(all_countries, ignore_index=True).to_csv('top_countries.csv', index=False)

    print("\n[OK] Wrote: ecology_summary.csv"
          f"{', top_hosts.csv' if all_hosts else ''}"
          f"{', top_sources.csv' if all_sources else ''}"
          f"{', top_countries.csv' if all_countries else ''}")

if __name__ == '__main__':
    main()
