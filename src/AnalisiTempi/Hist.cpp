#include <cmath>
#include <vector>
#include <iostream>
#include <string>

#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TLine.h>

#include "Hist.h"

void histFit(TH1F *h, const std::string &fname) {
    int nBins = h->GetNbinsX();
    double Max_est = h->GetMaximum();
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();
    double avg = h->GetMean();

    TF1 *f = new TF1("f", "[0] * exp(-[1]*x) + [2]", xmin, xmax);
    f->SetParameters(Max_est, 1/avg, 0);
    f->SetNpx(1000);
    h->Fit(f, "ILS Q"); 

    // Calcolo chi²
    double chi2 = f->GetChisquare();
    int ndf = f->GetNDF();

    //calcolo residui

    TGraph *gResidui = new TGraph();
    int point = 0;
    for (int i = 1; i <= nBins; i++) {
        double err = h->GetBinError(i);
        if (err > 0) {
            double x = h->GetXaxis()->GetBinCenter(i);
            double res = (h->GetBinContent(i) - f->Eval(x)) / err;
            gResidui->SetPoint(point++, x, res);
        }
    }

    // Disegno
    TCanvas *C = new TCanvas("c", "c", 800, 900);

    // Dividi in due pad verticali (asimmetrici)
    TPad *pad1 = new TPad("pad1", "fit",      0, 0.4, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "residui",  0, 0.0, 1, 0.3);

    pad1->SetBottomMargin(0.02);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.3);

    pad1->Draw();
    pad2->Draw();

    // Disegna il fit nel pad grande
    pad1->cd();
    h->SetStats(0);
    h->GetXaxis()->SetTitle("Differenze temporali [s]");
    h->GetYaxis()->SetTitle("Occorrenze [pure]");
    h->Draw();

    // Creazione della casella
    double x1 = 0.65, y1 = 0.55;
    double x2 = 0.9,  y2 = 0.85;
    TPaveText *pt = new TPaveText(x1, y1, x2, y2, "NDC");
    pt->SetFillColorAlpha(kWhite, 0.8); 
    pt->SetLineColor(kBlack);
    pt->SetTextFont(42);
    pt->SetTextSize(0.03);
    pt->SetTextAlign(12); 

    pt->AddText(Form("Parametri di fit:"));
    pt->AddText(Form("Parameter A : %.2f ± %.2f", f->GetParameter(0), f->GetParError(0)));
    pt->AddText(Form("Parameter R : %.2f ± %.2f Hz", f->GetParameter(1), f->GetParError(1)));
    pt->AddText(Form("Parameter B : %.2f ± %.2f", f->GetParameter(2), f->GetParError(2)));
    pt->AddText(Form("Chi2/NDF = %.2f", f->GetChisquare()/f->GetNDF()));

    pt->Draw(); 
    f->Draw("same");

    // Disegna i residui nel pad piccolo
    pad2->cd();
    gResidui->SetMarkerStyle(20);
    gResidui->SetMarkerSize(0.5);
    gResidui->GetXaxis()->SetTitle("Differenze temporali [s]");
    gResidui->GetYaxis()->SetTitle("Residui normalizzati [pure]");
    gResidui->Draw("AP");

    TLine *line = new TLine(xmin, 0, xmax, 0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->Draw();
    C->Update();

    // Salvataggio
    TFile file(fname.c_str(), "RECREATE");
    C->Write();
    file.Close();
    delete C;
}

void histMain(TTree* t, int n_ch) {

    // Variabili originali
    ULong64_t ch;
    ULong64_t time;

    t->SetBranchAddress("Channels", &ch);
    t->SetBranchAddress("Time", &time);

    ULong64_t N = t->GetEntries();

    std::vector<TH1F*> hist(n_ch);
    for (int i = 0; i < n_ch; i++)
        hist[i] = new TH1F(Form("h_ch%d", i), Form("Differenze temporali ch %d;time [s]", i), 100, 0, 4);

    // Riempimento istogrammi
    for (ULong64_t k = 0; k < N; k++){
        t->GetEntry(k);

        if (ch == 2147483648) continue;

        for (int bit = 0; bit < n_ch; bit++){
            if ((ch >> bit) & 1){
                ULong64_t start_time = time;
                ULong64_t stop_time = 0;

                for (ULong64_t h = k+1; h < N; h++){
                    t->GetEntry(h);
                    if ((ch >> bit) & 1){
                        stop_time = time;
                        break;
                    }
                }

                if (stop_time > start_time){
                    hist[bit]->Fill((stop_time - start_time) * 5e-9);
                }
            }
        }
    }

    for (int j = 0; j < n_ch; j++)
        histFit(hist[j], Form("FileRoot/h_ch%d.root", j));
}