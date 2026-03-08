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

#include "FlagClass.h"

struct BoolTrio {

    bool a, b, c;

    BoolTrio(bool x = false, bool y = false, bool z = false)
    : a(x), b(y), c(z) {}

    bool check_first (const BoolTrio&other) const {return a && other.a;}

    bool check_second (const BoolTrio&other) const {return b && other.b;}

    bool check_third (const BoolTrio&other) const {return c && other.c;}

    bool Or() const {return a  || b || c;}

    bool Sum(const BoolTrio&other) {return a && other.a || b && other.b || c && other.c;}

    void update(bool x, bool y, bool z) {
        a = x;
        b = y;
        c = z;
    }

};

void histFit(TH1F *h);

void histGraph(TTree *t);


#endif