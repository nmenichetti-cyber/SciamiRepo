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

#include "FlagClass.h"
#include "Hist.h"

void histGraph(TTree *t){

    //Allestimento del Tree

    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    //Numero di eventi

    Long64_t N = t->GetEntries();

    //Calibrazione del FIFO

    double c = 5e-9;

    //Trio di istogrammi

    TH1F *h_08 = new TH1F("Setup 08", "Tempi fra eventi successivi, setup08", 100, 0, 4e-6);

    TH1F *h_06 = new TH1F("Setup 06", "Tempi fra eventi successivi, setup06", 100, 0, 4e-6);

    TH1F *h_04 = new TH1F("Setup 04", "Tempi fra eventi successivi, setup04", 100, 0, 4e-6);

    //Loop sugli eventi

    for (Long64_t k = 0; k < N; k++){

        t->GetEntry(k);

        double t_1 = time;

        Flag Channel;

        BoolTrio FirstMask(Channel.IsTriple08(ch), Channel.IsTriple06(ch), Channel.IsTriple04(ch));

        if (Channel.IsResetKW(ch)) {continue;}
        if (!FirstMask.Or()) {continue;} 

        if (FirstMask.Or()) {

            Long64_t h = k+1;

            BoolTrio SecondMask;

            while (h < N) {

                t->GetEntry(h);
                Flag Temp_Ch;

                SecondMask.update(Temp_Ch.IsTriple08(ch),Temp_Ch.IsTriple06(ch),Temp_Ch.IsTriple04(ch));

                if (FirstMask.Sum(SecondMask)) break;

                h++;
            }

            if (h >= N) break;

            double t_2 = time;

            if (FirstMask.check_first(SecondMask)) {h_08->Fill(c*(t_2-t_1));}

            if (FirstMask.check_second(SecondMask)) {h_06->Fill(c*(t_2-t_1));}

            if (FirstMask.check_third(SecondMask)) {h_04->Fill(c*(t_2-t_1));}

            k = h-1;

            }
    }

    TCanvas *C_08 = new TCanvas();
    h_08->Draw();
    C_08->Update();

    TCanvas *C_06 = new TCanvas();
    h_06->Draw();
    C_06->Update();

    TCanvas *C_04 = new TCanvas();
    h_04->Draw();
    C_04->Update();

}