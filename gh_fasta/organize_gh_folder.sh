#!/usr/bin/env bash
set -euo pipefail

# 1) φτιάξε φακέλους
mkdir -p data matrices plots scripts

# 2) data: fasta + raw tsv (όχι τα _matrix.tsv)
mv -n -- *.fasta data/ 2>/dev/null || true
mv -n -- gh18.tsv gh19.tsv gh23.tsv gh48.tsv gh75.tsv gh80.tsv data/ 2>/dev/null || true

# 3) matrices: μόνο τα distance matrices
mv -n -- *_matrix.tsv matrices/ 2>/dev/null || true

# 4) plots: όλα τα παραγόμενα pdf
mv -n -- *_dendrogram.pdf *_nj.pdf plots/ 2>/dev/null || true
# Αν υπάρχει Rplots.pdf, μετονομάσ’ το για να μη σε μπερδεύει
if [ -f Rplots.pdf ]; then
  mv -n Rplots.pdf plots/legacy_Rplots.pdf
fi

# 5) scripts: ό,τι .R και .sh (εκτός του παρόντος)
for s in *.R *.sh; do
  [ "$s" != "organize_gh_folder.sh" ] && [ -f "$s" ] && mv -n -- "$s" scripts/
done

# 6) .gitignore (αν δεν υπάρχει)
if [ ! -f .gitignore ]; then
  cat > .gitignore <<'EOF'
# R junk
.Rhistory
.RData
.Rproj.user/

# PDFs and plots (tracked in /plots only)
Rplots.pdf

# Temporary files
*.tmp
*.swp
EOF
fi

echo "✅ Done. New layout:"
find . -maxdepth 2 -type d -print | sed 's|^\./||'
