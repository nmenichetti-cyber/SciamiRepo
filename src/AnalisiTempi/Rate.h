#ifndef RATE_H
#define RATE_H

#include <vector>
#include <array>
#include <stdexcept>

#include <TTree.h>

struct TelData{
    std::vector<double> tel_08;
    std::vector<double> tel_06;
    std::vector<double> tel_04;

    std::vector<double>& operator[](size_t idx) {
        switch(idx) {
            case 0: return tel_08;
            case 1: return tel_06;
            case 2: return tel_04;
            default: throw std::out_of_range("Indice fuori dai limiti");
        }
    }

    size_t size() const { return 3; }
};


void rateGraph(TTree *t, double Dt);

#endif