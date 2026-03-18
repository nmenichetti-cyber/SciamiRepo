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
    return (ev == 1 || ev == 16 || ev == 256 || ev == 3 || ev == 5 || ev == 9 || ev == 48 || ev == 80 || ev == 144 || ev == 9 || ev == 768 || ev == 1280 || ev == 2304);
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
            continue;   //Skip this line? Dipende se vuoi includere il reset
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
    int BUF[1] = {0};
    int Rb = 0;  
    int rt = 0;
    vector<double> RH_err;
    vector<double> TSH_err(TSH.size(), 0.0);
    int triplasecondopmt = 0;
    int difference = 0;

    for(int j = 0; j < (int)V_Ts.size(); j++){ 
        if(V_Ts[j] != cl){
            if(isEventoValido(V_Ev[j])){
                for(int i = j+1; i < (int)V_Ts.size(); i++){
                    if(V_Ts[i] == cl){
                        rt=0;
                        break;
                    }
                    if(V_Ts[i] <= V_Ts[j] + 42 && isEventoValido(V_Ev[i]) && V_Ev[i] != V_Ev[j] && V_Ev[i] != BUF[0]){
                            rt++;
                            BUF[0] = V_Ev[i];
                            if(rt == 2){
                                Rt++;
                                rt = 0;
                                if( BUF[0] == 16 || BUF[0] == 48 || BUF[0] == 80 || BUF[0] == 144){
                                triplasecondopmt++;            
                                }
                                BUF[0]=0;
                                j=i;
                                break;
                            }
                    }
                    else if(V_Ts[i] > V_Ts[j] + 42) {
                        rt = 0;
                        break;
                    }
                    else{
                        continue;
                    };
                }
            }
        }
        else{
            Rs++;
            if(Rs >= n*671){
                n++;
                Rth = (Rt - Rb);
                tsh++;
                TSH.push_back(tsh);
                RH.push_back((Rth));
                RH_err.push_back(sqrt(Rth));
                Rb = Rt;
            }
            cl++;
        }
    }

    difference = Rt -triplasecondopmt;

    long double Ts = Rs * 5.3687091;  // secondi?
    cout << "Rs: " << Rs << ", Rt: " << Rt << endl;
    cout << "Ts: " << Ts << " s, Th: " << Ts/3600 << " h" << endl;
    long double RT = ((Rt)/Ts)*1000;
    cout << "RATE COINCIDENZE TRIPLE: "<<RT << " conteggi all'ora" << endl;
    cout << "Triple secondo PMT accidentali: "<< difference << endl;
    
    if (!TSH.empty() && !RH.empty()) {
        TCanvas *c1 = new TCanvas("c1", "Rate triple", 800, 600);
        c1->SetGrid();
        
        TGraphErrors *gr = new TGraphErrors(RH.size(), TSH.data(), RH.data(), TSH_err.data(), RH_err.data());
        gr->SetTitle("Rate coincidenze triple nel tempo;Ore (h);Rate (Counts/h)");
        gr->SetMarkerStyle(20);
        gr->SetMarkerColor(kBlue);
        gr->SetLineColor(kRed);
        gr->Draw("AP");
        
        c1->SaveAs("rate_triple.png");
        c1->Draw();
    } else {
        cout << "Nessun dato da plottare!" << endl;
    }
}

// da riportare incertezza e cambiare da mHz a conteggi orari