# Rechnung der Rohdaten-Ableitung in der Cloud (AWS, Hetzner, Vultr; Stand 19.10.2026)

Ziel: die Teilnehmer ab Rang 1001 der festen Zufallsreihenfolge (`../analytic_random_order.tsv`) in der Cloud rechnen, waehrend der Mac die Raenge 1 bis 1000 abarbeitet.

Ablauf
1. `make_cloud_variants.py` erzeugt `ggir_one_cloud.R` und `swan_one_cloud.sh` aus den lokal benutzten Skripten (nur Pfade/Hilfsprogramme, Diff wird ausgegeben).
2. Docker-Image aus `Dockerfile` (R 4.6.1, GGIR 3.3-6 und alle R-Pakete wie auf dem Mac, Python 3.8 mit SWaN_accel 1.18 und gleichen Versionen).
3. Plattformvergleich: `make_crosscheck_list.py` (20 lokal fertige Teilnehmer), in der Cloud neu rechnen, `cross_check.R` vergleicht Epoche fuer Epoche. Erst bei Uebereinstimmung geht der Hauptlauf los.
4. Hauptlauf: `make_list.py <von> <bis>` erzeugt die Liste, `queue_worker.sh <liste> <slots>` rechnet; jeder Slot laedt sein Archiv von ftp.cdc.gov, rechnet GGIR und SWaN parallel, behaelt nur die Ergebnisse, loescht das Archiv.
5. Ergebnisse (~10 MB je Teilnehmer) nach S3, von dort auf die LaCie; die `done`-Marker werden zuletzt gesetzt, dann sieht der lokale Dispatcher die Teilnehmer als fertig.
6. Instanz beendet sich selbst (Shutdown-Verhalten terminate, Zeitlimit).
