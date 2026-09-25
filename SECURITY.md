# Politique de sécurité

## Versions prises en charge

| Version | Correctifs de sécurité |
|---|---|
| 0.1.x | oui |

## Signaler une vulnérabilité

Merci de **ne pas ouvrir d'issue publique**. Utilisez le signalement privé de GitHub :
onglet **Security**, puis **Report a vulnerability**. Une première réponse est donnée sous
7 jours.

## Périmètre

orfeval est un outil d'analyse local. Points d'attention :

- **Fichiers d'entrée** : les fichiers GenBank/FASTA sont lus par Biopython. N'analysez pas
  de fichiers d'origine inconnue sur une machine sensible.
- **Réseau** : seule la commande `fetch` accède au réseau (E-utilities du NCBI, en HTTPS).
  Une clé d'API NCBI éventuelle se passe par la variable d'environnement `NCBI_API_KEY` et
  n'est jamais écrite dans les résultats.
- **Dashboard** : Streamlit n'a pas d'authentification et accepte des fichiers importés. Il
  est prévu pour un usage local : ne l'exposez pas sur Internet sans protection (proxy
  authentifié, réseau privé).
- **Docker** : l'image s'exécute avec un utilisateur non root.
