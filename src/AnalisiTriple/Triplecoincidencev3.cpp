#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <cstring>

#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TAxis.h>
#include <TTree.h>

// #include "FlagClass.h"  // Commentato se non serve
// #include "Hist.h"       // Commentato se non serve

using namespace std;

// Versione con controllo più elegante usando una funzione di supporto
bool isEventoValido(unsigned int ev) {
    return (ev == 1 || ev == 3 || ev == 12 || ev == 4 || ev == 48 || ev == 16);
}

void ratetriple(const char* filename = "Fiforead_example.txt") {
    
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "Errore: impossibile aprire il file " << filename << endl;
        return;
    }

    // Crea i vettori
    vector<unsigned int> V_Ev;
    vector<unsigned int> V_Ts;
    
    // Parametri
    unsigned int cl = 2147483648;  // 2^31
    bool primoResetTrovato = false;
    
    // Leggi il file riga per riga
    unsigned int first, second;
    while (file >> first >> second) {
        
        // Controlla se questo è il primo reset
        if (!primoResetTrovato && second == cl) {
            primoResetTrovato = true;
            cout << "Primo reset trovato alla riga con Ts = " << second << endl;
            // Non aggiungiamo questo reset? Dipende se vuoi includerlo o no
            // Se vuoi includere anche il reset, decommenta le prossime 2 righe:
            V_Ev.push_back(first);
            V_Ts.push_back(second);
            //continue;   Skip this line? Dipende se vuoi includere il reset
        }
        
        // Se abbiamo già trovato il primo reset, aggiungi tutti gli eventi
        if (primoResetTrovato) {
            V_Ev.push_back(first);
            V_Ts.push_back(second);
        }
    }

    file.close();

    if (!primoResetTrovato) {
        cout << "Attenzione: nessun reset trovato nel file!" << endl;
        return;
    }

    cout << "File: " << filename << endl;
    cout << "Lette " << V_Ev.size() << " righe (dopo il primo reset)" << endl;

    //unsigned int cl = 2147483648;  // 2^31
    unsigned int Rs = 0, Rt = 0;
    int n = 1;
    int tsh = 0;
    int Rth = 0;
    vector<double> RH;
    vector<double> TSH;
    int Rb = 0;

    for(int j = 0; j < (int)V_Ts.size() - 2; j++){ 
        if(V_Ts[j] != cl){
            if(isEventoValido(V_Ev[j])){
                if(V_Ts[j+1] <= V_Ts[j] + 208 && isEventoValido(V_Ev[j+1]) && V_Ev[j+1] != V_Ev[j]){
                    if(V_Ts[j+2] <= V_Ts[j] + 208 && isEventoValido(V_Ev[j+2]) && V_Ev[j+2] != V_Ev[j]){ // controlla per verificare che ci sia una tripla poco dopo
                        Rt++;
                        j += 2;  // +2 perchè il loop farà un altro +1
                    }
                }
            }
        }
        else{
            Rs++;
            if(Rs >= n*680){
                n++;
                Rth = (Rt - Rb);
                tsh++;
                TSH.push_back(tsh);
                RH.push_back((Rth)/3.6);
                Rb = Rt;
            }
            cl++;
        }
    }

    long double Ts = Rs * 5.3687091;  // secondi?
    cout << "Rs: " << Rs << ", Rt: " << Rt << endl;
    cout << "Ts: " << Ts << " s, Th: " << Ts/3600 << " h" << endl;
    long double RT = ((Rt)/Ts)*1000;
    cout << "RATE COINCIDENZE TRIPLE: "<<RT << " mHz" << endl;
    
    if (!TSH.empty() && !RH.empty()) {
        TCanvas *c1 = new TCanvas("c1", "Rate triple", 800, 600);
        c1->SetGrid();
        
        TGraph *gr = new TGraph(RH.size(), TSH.data(), RH.data());
        gr->SetTitle("Rate coincidenze triple nel tempo;Ore (h);Rate (mHz)");
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kBlue);
        gr->SetLineColor(kRed);
        gr->Draw("APL");
        
        c1->SaveAs("rate_triple.png");
        c1->Draw();
    } else {
        cout << "Nessun dato da plottare!" << endl;
    }
}