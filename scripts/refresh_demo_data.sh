#!/usr/bin/env bash
# Retélécharge les génomes de démonstration depuis le NCBI et compare les sommes SHA-256
# à celles versionnées dans data/demo/SHA256SUMS.
#
# Usage : NCBI_EMAIL=vous@exemple.org bash scripts/refresh_demo_data.sh
#
# Une différence de somme n'est pas une erreur en soi : le NCBI peut mettre à jour
# l'annotation d'un enregistrement sans changer son numéro de version. Les résultats
# documentés correspondent aux fichiers versionnés dans le dépôt.
set -euo pipefail

cd "$(dirname "$0")/.."
ACCESSIONS=(NC_001416.1 NC_000908.2)
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

orfeval fetch "${ACCESSIONS[@]}" --outdir "$TMP_DIR" --force

cd "$TMP_DIR"
if sha256sum --check --quiet "$OLDPWD/data/demo/SHA256SUMS"; then
    echo "Les fichiers du NCBI sont identiques à ceux du dépôt."
else
    echo "Différence détectée : l'enregistrement NCBI a changé depuis la version du dépôt." >&2
    echo "Fichiers téléchargés conservés dans data/raw/ pour comparaison." >&2
    mkdir -p "$OLDPWD/data/raw"
    cp ./*.gb.gz "$OLDPWD/data/raw/"
    exit 1
fi
