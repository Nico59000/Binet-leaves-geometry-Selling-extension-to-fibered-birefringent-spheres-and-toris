# Programme KEP — Handoff de reprise complète après KEP-R24

## Point actif

\[
\texttt{KEP-MB}\Sigma\texttt{-4}
=
\texttt{OPEN\_786\_TRANSVERSE\_CLASSES\_PLUS\_}
\texttt{17\_ENDPOINT\_ROUTES\_5068\_SUBLEAVES}.
\]

État calculatoire à reprendre :

- 786 classes transversales ;
- 6413 boîtes transversales sérialisées ;
- 16 routes endpoint split et une route frontière ;
- 5068 sous-feuilles endpoint ;
- 4052 représentants géométriques exacts S4.

## Fichiers prioritaires

1. `stages/kep_r24_formal_bundle.zip` ;
2. `certificates/transverse/persistent_frontier_6413.json` ;
3. `certificates/endpoints/endpoint_open_subleaves_5068.json` ;
4. `certificates/audit/final_5068_exact_S4_orbit_ledger.json` ;
5. `certificates/decision/KEP_MB_SIGMA_4_decision.json` ;
6. les sources du dossier `src/`.

## Garde-fous

- la couche Fibonacci–Binet est seulement un ordonnanceur 3+2, sans transfert analytique ;
- ne jamais promouvoir un résultat Krawczyk `SPLIT` ;
- conserver séparément la route frontière ;
- KEP-MC et la fermeture globale restent ouvertes ;
- fermer KEP-MB Sigma-4 seulement après extinction de toutes les files ou contre-certificat explicite.

## Intégrité

Les bundles byte-exacts KEP-R4 à KEP-R24, les sources principales et les
références Binet utilisées à R24 sont inclus. KEP-R0 à KEP-R3 restent
reconstruits par handoff.
