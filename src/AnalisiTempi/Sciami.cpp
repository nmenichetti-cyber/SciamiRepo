#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <array>
#include <cstring>

#include <TGraph.h>
#include <TCanvas.h>

#include "FlagClass.h"
#include "Rate.h"
#include "Hist.h"

//Crea il file .root con i dati del FIFO

void FileToTree (const char * fname, const char * option = "RECREATE") {

    //Apertura del file con test

    std::ifstream data(fname);

    if (!data) {cout << "Errore di apertura \n" << endl;}

    //Inizializzazione del File di output e del Tree

    TFile OutFile("Data.root", option);

    TTree Tree("T", "Contiene i dati del DE10-NANO");

    //Variabili, canali scattati e tempi in digit

    ULong64_t ch = 0;

    ULong64_t time = 0;

    //Creazione dei branch

    Tree.Branch("Channels", &ch, "Channel/l");

    Tree.Branch("Time", &time, "Time/l");

    //Riempimento del TTree

    while (data >> ch >> time){

        Tree.Fill() ;

        //Per controllo puoi printare le entries del file, ma è sconsigliato per file lunghi

        /*

        std::cout << ch << "," << time <<endl;

        */

    }

    Tree.Write();
    OutFile.Close();
}

//Fa partire un'analisi a scelta; le keyword sono "rate" e "hist"

void Sciami(const char * path, const char* kw){

    //Apertura file e ottenimento del TTree

    TFile *f = TFile::Open(path);
    TTree *t = (TTree*)f->Get("T");

    //Stringa di scelta
    if (strcmp(kw, "rate") == 0) {rateGraph(t, 100);}

    if (strcmp(kw, "hist") == 0) {histGraph(t);}

}