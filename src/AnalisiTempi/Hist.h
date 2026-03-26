#ifndef HIST_H
#define HIST_H

#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <cstring>

#include <TGraph.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TTree.h>
#include <TAxis.h>



void histFit(TH1F *h, const std::string &fname);

void histMain(TTree *t, int n_ch);


#endif