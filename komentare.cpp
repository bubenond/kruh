int main() {
    // 1. Nastavení Modelu a Sítě

    // Inicializace fyzikálních konstant (tepelná vodivost, okrajové h) a geometrie (Lx, Nx, dx).
    // Vytvoření indexovací funkce pro převod 3D souřadnic (i, j, k) na 1D index.
    // Vytvoření a inicializace pole T (Teplotní pole) s počátečním odhadem.
    // Vytvoření a nastavení pole Q (objemový zdroj tepla - Jouleovo teplo) v definované oblasti.


    // 2. Iterativní Řešení Stacionárního Stavu (Gauss–Seidel)
    // Nastavení numerické tolerance a maximálního počtu iterací

    // HLAVNÍ ITERAČNÍ SMYČKA: Dokud není dosaženo konvergence nebo maxIterace
    // Sledujeme maximální změnu teploty v iteraci

    // Aktualizace Teplot Vnitřních Uzlů:
    // PRO VŠECHNY VNITŘNÍ UZLY (i, j, k):
    // Vypočti NOVOU teplotu T_new pomocí diskretizované Poissonovy rovnice (FDM)
    // T_new = Fce(Šest sousedů, Vliv zdroje Q).
            
    // Aktualizuj T s T_new (Gauss-Seidel) a sleduj maxChange.

    // Aplikace a Aktualizace Okrajových Podmínek (Robinovy)
    // PRO VŠECHNY OKRAJOVÉ UZLY (x=0, x=Lx, y=0 atd.):
    // T[okraj] = VÝPOČET Z ROBINOVY BC (závislé na T_vnitřní_soused, konvekci a T_okolí).
            


    // 3. Výstup Výsledků 
    // Zjisti a vytiskni minimální a maximální teplotu v nalezeném ustáleném stavu.
    // Uložení dat:
    // výstupní soubory (CSV)
    // vizualizace výsledů pomocí grafů a obrázků
    
    return 0;
}