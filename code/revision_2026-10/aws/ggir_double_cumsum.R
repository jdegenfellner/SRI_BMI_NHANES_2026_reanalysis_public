# Emuliert in GGIRs Autokalibrierung die Rechengenauigkeit von R auf Apple-arm64 (long double = double):
# g.calibrate bildet Fenstermittel ueber cumsum(); R rechnet cumsum auf x86 mit 80-Bit-Akkumulator, auf Apple-arm64 mit 64 Bit.
# Bei flachen Signalstuecken entscheidet das, welche 10-s-Fenster GGIR als Duplikate entfernt. Hier wird cumsum in g.calibrate in doppelter Genauigkeit akkumuliert.
suppressMessages(library(Rcpp))
Rcpp::cppFunction("NumericVector cumsum_dbl(NumericVector x) { int n = x.size(); NumericVector y(n); double s = 0; for (int i = 0; i < n; i++) { s += x[i]; y[i] = s; } return y; }")
.gc <- GGIR::g.calibrate
.e <- new.env(parent = environment(.gc)); .e$cumsum <- function(x) cumsum_dbl(x)
environment(.gc) <- .e
utils::assignInNamespace("g.calibrate", .gc, "GGIR")
