#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <tuple>

#include <TGraph.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TAxis.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TLegend.h>

#include "Rate.h"

//Definizione del metodo fill per la classe container

void Container::Fill(double y_val, double yerr_val){
    y.push_back(y_val);
    y_err.push_back(yerr_val);
}

//Blocco per fare il grafico +  fit + larghezza banda

std::tuple<TGraphErrors*, float> rateGraph (const std::vector<double> &x,
                         const std::vector<double> &y,
                         const std::vector<double> &x_err,
                         const std::vector<double> &y_err){

    //Grafico dei rate

    TGraphErrors *g = new TGraphErrors(x.size(), x.data(), y.data(), x_err.data(), y_err.data());

    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.1);
    g->SetTitle("Rate vs Tempo");

    g->GetXaxis()->SetTitle("Tempo [h]");
    g->GetYaxis()->SetTitle("Rate [Hz m^{-2} s^{-1}]");

    double xmax = *std::max_element(x.begin(), x.end());

    //Fit con una costante

    TF1 *f = new TF1("f", "[0]", 0, xmax);
    f->SetParameter(0, 5);
    g->Fit(f, "Q");

    //Calcolo variazione percentuale

    float p = f->GetParameter(0);
    std::vector<float> diff;
    for (long unsigned int h = 0; h < y.size(); h++){

        diff.push_back(std::abs(100*(y[h]-p)/p));

    }

    float max = *std::max_element(diff.begin(), diff.end());

    //Tupla di output

    std::tuple<TGraphErrors*, float> res;

    res = make_tuple(g, max);

    return res;
}

//Blocco per il calcolo dei rate e output del programma

void rateMain(TTree *t, int n_ch, double Dt){

    //Inizializzazione del TTree

    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    //Numero di eventi

    ULong64_t N = t->GetEntries();

    //Variabili di timing : tempo passato e tempo totale, con vettore dei tempi e errori

    double el = 0.;
    double tot = 0.;

    std::vector<double> timings;
    std::vector<double> timings_err;

    //Variabili di conteggio e vettori di contenimento dei rate

    std::vector<double> counter(n_ch, 0.);
    std::vector<Container> data(n_ch);

    timings.reserve(N / 10);
    timings_err.reserve(N / 10);

    //Indice del loop e step (serve dopo è il numero di telescopi)

    ULong64_t k = 0;
    int step = n_ch / 3;

    //Loop sugli eventi

    while (k < N) {

        t->GetEntry(k);

        //Caso di parola di riscrittura, aumento le variabili temporali

        if (ch == 2147483648) {
            el += 5;
            tot += 5./3600;
        }

        //Caso diifferente : loop sui bit e aumento i conteggi

        if (ch != 2147483648){

            ULong64_t mask = ch;

            while (mask){
                int j = __builtin_ctzll(mask);
                if (j < n_ch) counter[j] += 1;
                mask &= (mask - 1);
            }
        }

        //Se raggiungo la finestra desiderata, riempio i vettori

        if (el >= Dt){

            el = 0.;

            timings.push_back(tot);
            timings_err.push_back(5./3600);

            for (int l = 0; l < n_ch; l++){
                double rate = counter[l]/Dt;
                data[l].Fill(rate, std::sqrt(counter[l])/Dt);
                counter[l] = 0.;
            }
        }

        k++;
    }

    //Matrice accettanze

    std::vector<std::vector<double>> A = {
        {0.97, 1, 0.92},
        {0.36, 1, 0.42},
        {0.61, 1, 0.61}
    };

    //Supperfici e angoli solidi sottesi

    std::vector<float> S = {0.1218, 0.08, 0.2};
    std::vector<float> O = {0.776, 0.171, 0.495};

    //Correzione dell'efficienza (da fare un plot in seguito)

    for (int m = 0; m < n_ch; m += 4){

        auto& T   = data[m].y;
        auto& D12 = data[m+1].y;
        auto& D23 = data[m+2].y;
        auto& D31 = data[m+3].y;

        const std::vector<double>& A_v = A[(m / 4) % 3];

        for (size_t e = 0; e < T.size(); e++){

            if (D12[e] > 0 && D23[e] > 0 && D31[e] > 0){

                double eff_T  = (T[e]/(D12[e]*A_v[0])) * (T[e]/(D23[e]*A_v[1])) * (T[e]/(D31[e]*A_v[2]));
                double eff_12 = (T[e]/(D23[e]*A_v[0])) * (T[e]/(D31[e]*A_v[1]));
                double eff_23 = (T[e]/(D12[e]*A_v[2])) * (T[e]/(D31[e]*A_v[1]));
                double eff_31 = (T[e]/(D12[e]*A_v[2])) * (T[e]/(D23[e]*A_v[0]));

                data[m].y[e]   /= eff_T * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m].y_err[e] /= eff_T * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m+1].y[e] /= eff_12 * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m+1].y_err[e] /= eff_12 * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m+2].y[e] /= eff_23 * O [(m / 4) % 3] * S[(m / 4) % 3] ;
                data[m+2].y_err[e] /= eff_23 * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m+3].y[e] /= eff_31 * O [(m / 4) % 3] * S[(m / 4) % 3];
                data[m+3].y_err[e] /= eff_31 * O [(m / 4) % 3] * S[(m / 4) % 3];

                }
            }
        }

    //Riempimento dei file e output

    int n_files = n_ch / step;
    std::vector<std::string> labels = {"1&2&3", "1&2", "2&3", "1&3"};

    for (int n = 0; n < n_files; n++) {

        TCanvas* c = new TCanvas(Form("c_setup_%d", n), "Rate vs Tempo", 900, 700);
        TLegend* legend = new TLegend(0.65, 0.7, 0.88, 0.88);
        legend->SetHeader("Tipo di evento", "C");

        // Calcolo dinamico di ymin/ymax su tutti i canali del setup
        double ymin = 1e9, ymax = -1e9;
        for (int ch_index = n*step; ch_index < (n+1)*step; ch_index++) {
            for (size_t i = 0; i < data[ch_index].y.size(); i++) {
                if (data[ch_index].y[i] > 0) {
                    ymin = std::min(ymin, data[ch_index].y[i] - data[ch_index].y_err[i]);
                    ymax = std::max(ymax, data[ch_index].y[i] + data[ch_index].y_err[i]);
                }
            }
        }

        ymin = std::max(0.0, ymin * 0.9);
        ymax = ymax * 1.1;

        bool first = true;

        for (int ch_index = n*step; ch_index < (n+1)*step; ch_index++) {

            std::tuple<TGraphErrors*, float> res = rateGraph(timings, data[ch_index].y, timings_err, data[ch_index].y_err);

            TGraphErrors* g = std::get<0>(res);

            g->SetMarkerColor((ch_index % step) + 1);
            g->SetLineColor((ch_index % step) + 1);

            g->SetMinimum(ymin);
            g->SetMaximum(ymax);

            if (first) {
                g->Draw("AP");
                first = false;
            } else {
                g->Draw("P SAME");
            }

            cout << std::get<1>(res) <<endl;

            // Fit del grafico
            TF1* f = g->GetFunction("f");
            if (f) {
                double p = f->GetParameter(0);
                double ep = f->GetParError(0);
                legend->AddEntry(g, Form("%s : %.2f #pm %.2f Hz m^{-2} str^{-1} ", labels[ch_index % step].c_str(), p, ep), "p");
            } else {
                // Se non c'è fit, aggiungi solo il canale
                legend->AddEntry(g, labels[ch_index % step].c_str(), "p");
            }
        }

        //Salvataggio e fine programma

        legend->Draw();
        c->Update();

        std::string name = "FileRoot/Rates_setup" + std::to_string(n+1) + ".root";
        TFile f(name.c_str(), "RECREATE");
        c->Write();
        f.Close();

        delete c;
    }
}






