#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <cstring>

#include <TGraph.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TAxis.h>
#include <TGraphErrors.h>

#include "FlagClass.h"
#include "Rate.h"

void rateGraph(TTree *t, double Dt){

    //Allestimento del Tree

    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    //Numero di eventi

    Long64_t N = t->GetEntries();

    //Variabili temporali

    double elapsed_time = 0.;
    double total_time = 0.;

    //Variabili di conteggio

    std::array<int,3> triples = {0,0,0};

    std::array<int,3> doubles = {0,0,0};

    //Indice del loop while

    int k = 0;

    //Strutture e vettori per i plot vari

    TelData rate;

    TelData err_rate;

    TelData eff;

    TelData err_eff;

    std::vector<double> timings;

    std::vector<double> timings_err;

    //Elaborazione e  riempimento

    while (elapsed_time < Dt && k < N){

        t->GetEntry(k);

        Flag Channel;

        if (Channel.IsResetKW(ch)) {

            elapsed_time += 5;
            total_time += 5; 

            if (elapsed_time >= Dt){

                elapsed_time *= 0;

                timings.push_back(total_time);

                timings_err.push_back(5./3600);

                for (int i = 0; i < 3; i++){

                if (doubles[i] != 0){

                    rate[i].push_back(triples[i]/Dt);

                    err_rate[i].push_back(sqrt(triples[i])/Dt);

                    double r = (double)triples[i]/(double)doubles[i];

                    eff[i].push_back(r);

                    err_eff[i].push_back(sqrt(r * (1-r)/(double)doubles[i]));


                }

                    triples[i] *= 0;
                    doubles[i] *= 0;
                }

            }

        }

        if (Channel.IsTriple08(ch)) {triples[0] += 1;}

        if (Channel.IsTriple06(ch)) {triples[1] += 1;}

        if (Channel.IsTriple04(ch)) {triples[2] += 1;}

        if (Channel.IsDouble08(ch)) {doubles[0] += 1;}

        if (Channel.IsDouble06(ch)) {doubles[1] += 1;}

        if (Channel.IsDouble04(ch)) {doubles[2] += 1;}

        k += 1;

    }

    

    //Grafici per i rate in funzione del tempo/efficienze

    for (size_t j = 0; j < rate.size(); j++) {
        // Grafico Rate vs tempo

        std::cout << "rate size = " << rate[j].size() << std::endl;
        std::cout << "band = " << calculate_band(rate[j]) << std::endl;

        TGraphErrors* g = new TGraphErrors(timings.size(), timings.data(), rate[j].data(), timings_err.data(), err_rate[j].data());
        g->SetMarkerStyle(21);
        g->SetMarkerColor(kRed);
        g->SetLineColor(kBlack);
        g->SetTitle("Rate vs tempo");
        g->GetXaxis()->SetTitle("Tempo [h]");
        g->GetYaxis()->SetTitle("Rate [Hz]");

        // Canvas unico per ogni grafico
        TCanvas* c = new TCanvas(Form("c_%zu", j), Form("Canvas_%zu", j), 800, 600);
        g->Draw("AP");
        c->Update();

        // Grafico Efficienza vs tempo 
        TGraphErrors* g_e = new TGraphErrors(timings.size(), timings.data(), eff[j].data(), timings_err.data(), err_eff[j].data());
        g_e->SetMarkerStyle(21);
        g_e->SetMarkerColor(kBlue); // cambia colore per distinguere
        g_e->SetLineColor(kBlack);
        g_e->SetTitle("Efficienza vs tempo");
        g_e->GetXaxis()->SetTitle("Tempo [h]");
        g_e->GetYaxis()->SetTitle("Efficienza [%]");

        TCanvas* c_e = new TCanvas(Form("c_e_%zu", j), Form("Canvas_e_%zu", j), 800, 600);
        g_e->Draw("AP");
        c_e->Update();

    }

}

double calculate_band(const std::vector<double>& v) {

    ULong64_t N = v.size();

    double sum = 0;

    for (ULong64_t k = 0; k < N; k++){

        sum += v[k];

    }

    double avg = sum /N;

    std::vector<double> diff;

    for (ULong64_t l = 0; l < N; l++){

        diff.push_back(std::abs(v[l]-avg));

    }

    auto bandwidth = std::max_element(diff.begin(), diff.end());

    return *bandwidth;

}