#!/usr/bin/env python3
# =====================================================================
# Compute per-protein physiochemical properties for Cave vs IMG/VR phage
# proteins, from the amino-acid sequences (FASTA), with Biopython ProtParam.
# This is the upstream step for FigS10_physicochemical_features.R:
# it produces the feature table that the plotting script reads.
#
# INPUTS (amino-acid sequences, one FASTA per protein type per group):
#   Cave  : <SRC>/sulfurcave_environment.tar.gz   (12 sulfurcave_<type>.faa)
#   IMG/VR: <SRC>/all_environment/imgvr_<type>.faa (12 files, 43K-1.9M seqs)
# Cave = ALL proteins; IMG/VR = seeded reservoir subsample of 3000 per type
# (the full IMG/VR set is far too large to compute/plot whole). seed = 123.
#
# OUTPUT (full per-protein feature table):
#   <repo>/viromics/cave_imgvr_features_raw.csv
#   columns: protein_id (FASTA header token), origin, protein, then all
#   features + the 20 AA percents.
#   read directly by plotting/FigS10_physicochemical_features.R
#
# FEATURES (Bio.SeqUtils.ProtParam.ProteinAnalysis): sequence length,
#   molecular weight, isoelectric point, aromaticity (F+W+Y fraction),
#   instability index, GRAVY hydrophobicity, mean flexibility (Vihinen
#   sliding-window index), net charge at pH 1.0, helix/turn/sheet fractions,
#   and per-amino-acid composition.
#
# NB on the *_Percent scale: values are stored as amino_acids_percent * 100.
#   In the Biopython version used here amino_acids_percent is already on a
#   0-100 scale, so the stored columns are (percent * 100) = (fraction * 10000).
#   Downstream: divide a *_Percent column by 100 for percent, by 10000 for a
#   0-1 fraction. (The plotting script relies on this convention.)
# =====================================================================
import os, glob, random, re, tarfile, tempfile, warnings
import numpy as np, pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
warnings.filterwarnings("ignore")

# --- paths (see README, "Running the scripts") ------------------------
REPO      = os.environ.get("SULFURCAVE_REPO", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
EXTERNAL  = os.environ.get("SULFURCAVE_EXTERNAL", os.path.join(REPO, "external_data"))
SRC       = os.path.join(EXTERNAL, "2026_analysis/SulfurCave_imgvr_proteins/all")
CAVE_TAR  = os.path.join(SRC, "sulfurcave_environment.tar.gz")
IMGVR_DIR = os.path.join(SRC, "all_environment")
OUT_CSV   = os.path.join(REPO, "viromics/cave_imgvr_features_raw.csv")

N_IMGVR = 3000
random.seed(123)
STD      = set("ACDEFGHIKLMNPQRSTVWY")
AA_ORDER = list("AVLIMFWYRKHDESTNQGPC")

def clean(seq):
    return "".join(c for c in seq.upper() if c in STD)

def iter_fasta(path):
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None: yield name, "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
        if name is not None: yield name, "".join(buf)

def reservoir(path, k):
    """Seeded reservoir subsample of k (name, seq) pairs from a (large) FASTA."""
    res = []
    for i, (nm, sq) in enumerate(iter_fasta(path)):
        if i < k:
            res.append((nm, sq))
        else:
            j = random.randint(0, i)
            if j < k: res[j] = (nm, sq)
    return res

def features(seq):
    s = clean(seq)
    if len(s) < 10:
        return None
    pa = ProteinAnalysis(s)
    try:
        flex = pa.flexibility(); fm = float(np.mean(flex)) if len(flex) else np.nan
        h, t, sh = pa.secondary_structure_fraction(); aap = pa.amino_acids_percent
        d = {"ProteinLength": len(s), "ProteinWeight": pa.molecular_weight(),
             "IsoelectricPoint": pa.isoelectric_point(), "Aromaticity": pa.aromaticity(),
             "InstabilityIndex": pa.instability_index(), "Hydrophobicity": pa.gravy(),
             "FlexibilityMean": fm, "ChargeAtPh1": pa.charge_at_pH(1.0),
             "HelixFraction": h, "TurnFraction": t, "SheetFraction": sh}
        for aa in AA_ORDER:
            d[f"{aa}_Percent"] = aap.get(aa, 0.0) * 100.0   # see scale note in header
        return d
    except Exception:
        return None

def collect(path, origin, protein_type, pairs=None):
    rows = []
    it = pairs if pairs is not None else iter_fasta(path)
    for nm, sq in it:
        f = features(sq)
        if f:
            # protein_id = first whitespace-delimited token of the FASTA header
            f.update(protein_id=nm.split()[0], origin=origin, protein=protein_type)
            rows.append(f)
    return rows

rows = []

# --- Cave: all proteins (extract the tarball to a temp dir) -----------
with tempfile.TemporaryDirectory() as tmp:
    with tarfile.open(CAVE_TAR) as tar:
        tar.extractall(tmp)
    for f in sorted(glob.glob(os.path.join(tmp, "**", "sulfurcave_*.faa"), recursive=True)):
        pt = re.sub(r"^sulfurcave_|\.faa$", "", os.path.basename(f))
        rows += collect(f, "sulfurcave", pt)

# --- IMG/VR: seeded 3000-per-type reservoir subsample ----------------
for f in sorted(glob.glob(os.path.join(IMGVR_DIR, "imgvr_*.faa"))):
    pt = re.sub(r"^imgvr_|\.faa$", "", os.path.basename(f))
    rows += collect(None, "imgvr", pt, pairs=reservoir(f, N_IMGVR))

df = pd.DataFrame(rows)
lead = ["protein_id", "origin", "protein"]   # id first, then origin, protein
df = df[lead + [c for c in df.columns if c not in lead]]
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
df.to_csv(OUT_CSV, index=False)

print("shape:", df.shape)
print(df.groupby("origin")["protein"].value_counts().unstack(0).fillna(0).astype(int))
