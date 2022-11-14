# Phylogenetic Tree Analyzer
a seminar program by Adél Vancl, as part of _NPRG030/Programming I_ at Charles University, Faculty of Mathematics and Physics

## Program Specification (in Czech)
- Program dle definovaných _množin taxonů_ (systematicky pojmenovaných _listů_ stromu) v každém stromu **nalezne nejvyšší ohodnocení hrany** (_bootstrap_), které určuje pravděpodobnost správného rozdělení příbuznosti listů ve dvou **podstromech oddělených touto hranou**, kde jeden ze stromů obsahuje kromě zadaného _seed taxonu_ pouze _taxony_ popsané _kvantifikovanými množinovými operacemi_, a nalezenou hodnotu systematicky uloží do výstupního _CSV souboru_.

- **CLI programu** očekává jako argumenty soubor(y) fylogenetických stromů v závorkovém formátu (_newick_, _nexus_), název seed taxonu, jména množin a výčet (s _wildcardy_) názvů taxonů do množiny spadajících, zkoumané podmínky pro podstromy definované kvantifikovanými množinovými operacemi a výstupní soubor.

- Program pro každý soubor fylogenetického stromu vygeneruje abstraktní strukturu stromu s ohodnocenými hranami a názvy taxonů v listech. Zpracuje zadané výčty množin taxonů a pro každou ze zadaných podmínek **hledá hranu s nejvyšší hodnotou**, která odděluje podstrom obsahující seed taxon a kde všechny ostatní listy spňují řešenou podmínku. Nejvyšší hodnotu ukládá do CSV souboru do příslušného sloupce.

- Podmínky jsou zadány kvantifikovanými množinovými operacemi, kde lze **kromě obvyklých operací definovat i minimální nebo maximální četnost výskytu** prvku množiny.

- Ve výsledném CSV souboru každý **řádek odpovídá jednomu ze souborů stromu** a sloupce odpovídají argumenty zadaným podmínkám.

_This specification has been approved by Mgr. Rudolf Rosa, Ph.D. (seminar program supervisor) and Mgr. Anna M. G. Novák Vanclová (bioinformatics consultant and project owner)_
