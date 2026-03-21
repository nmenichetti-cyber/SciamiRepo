///////////////////////////////////////////////////////////////////////////////
//  DRS4Browser_v3.cpp
//  ====================
//  Browser interattivo per la visualizzazione delle forme d'onda campionate
//  dalla DRS4 Evaluation Board, a partire dai file XML prodotti da DRSOsc.
//
//  NOVITA' v3: TIMING DUALE (CFD + soglia fissa con correzione walk)
//  -----------------------------------------------------------------
//  Implementa due strategie di timing:
//    (A) Segnali NON clippati → CFD al 35% (immune al time walk)
//    (B) Segnali CLIPPATI     → Soglia fissa + correzione walk tramite
//                                lo slew rate massimo (derivata) del fronte
//  Il campo t_timing contiene SEMPRE il miglior tempo disponibile.
//  Il campo timing_method indica il metodo usato (0=nessuno, 1=CFD, 2=walk).
//
//  STRUTTURA DEL PROGRAMMA (in ordine di esecuzione):
//  ---------------------------------------------------
//  1. DEFINIZIONI: costanti, strutture dati (ChannelData, EventData)
//  2. ComputeSmoothedDerivative(): derivata numerica con filtro a media mobile
//  3. AnalyzeChannel():  estrae baseline, ampiezza, timing duale, integrale
//  4. ParseXML():        legge il file XML con macchina a stati finiti
//  5. DrawEvent():       disegna un evento nel canvas multi-pad
//  6. Funzioni globali:  Next(), Prev(), GoTo(), Overlay*(), Save(), Summary()
//  7. DRS4Browser_v3():  entry point — lancia il browser
//  8. main():            versione compilata (opzionale)
//
//  UTILIZZO:
//    root -l 'DRS4Browser_v3.cpp("file.xml")'
//    root -l 'DRS4Browser_v3.cpp("file.xml", "CHN1=PMT top,CHN2=PMT mid,CHN3=PMT bot,CHN4=NIM trigger")'
//
//  Autori: Luca & Claude — Laboratorio Sciami Estesi, Marzo 2026
///////////////////////////////////////////////////////////////////////////////


// ============================= INCLUDE ====================================
// Le librerie ROOT (Txxx.h) forniscono classi per la grafica e l'interfaccia.
// Le librerie C++ standard gestiscono I/O e strutture dati.

#include <TCanvas.h>         // TCanvas: finestra grafica ROOT
#include <TGraph.h>          // TGraph: grafico XY a punti/linea
#include <TLine.h>           // TLine: linea dritta tra due punti
#include <TLatex.h>          // TLatex: testo con formattazione LaTeX
#include <TLegend.h>         // TLegend: legenda del grafico
#include <TStyle.h>          // gStyle: stile globale dei grafici
#include <TPad.h>            // TPad: sotto-area del canvas
#include <TSystem.h>         // gSystem: interfaccia al sistema operativo
#include <TROOT.h>           // gROOT: oggetto globale ROOT
#include <TMath.h>           // TMath: funzioni matematiche ROOT
#include <TPaveText.h>       // TPaveText: box di testo con bordo
#include <TApplication.h>    // TApplication: gestione event loop (versione compilata)
#include "TAxis.h"              // TAxis: gestione assi dei grafici  
#include <TEllipse.h>        // TEllipse: ellissi e cerchi (per skymap)
#include <TH1D.h>            // TH1D: istogramma 1D (per distribuzioni angolari)
#include <TF1.h>             // TF1: funzione 1D (per fit)



#include <fstream>           // std::ifstream: lettura file
#include <iostream>          // std::cout, std::cerr: output a terminale
#include <sstream>           // std::istringstream: parsing stringhe
#include <string>            // std::string: stringhe C++
#include <vector>            // std::vector: array dinamico (per lista eventi)
#include <map>               // std::map: dizionario (per etichette canali)
#include <cstring>           // strncpy, memset, strlen: operazioni su stringhe C
#include <cstdlib>           // atoi: conversione stringa→intero
#include <cmath>             // sqrt, fabs: funzioni matematiche
#include <algorithm>         // std::min_element, std::max_element
#include <utility>           // std::pair, std::make_pair


// ========================= COSTANTI =======================================
// Parametri fissi del hardware e dell'analisi. Per modificare il
// comportamento del programma, cambiare questi valori e ricompilare.

const int MAX_SAMPLES  = 1024;    // Numero di celle del chip DRS4 → 1024 campioni per canale.
                                  // NON MODIFICARE (è una proprietà del chip).

const int MAX_CHANNELS = 4;      // La DRS4 Eval Board ha esattamente 4 canali (SMA: CH1-CH4).
                                  // NON MODIFICARE.

const int NBL_SAMPLES  = 50;     // Quanti campioni all'inizio della finestra usare per
                                  // calcolare la baseline. A 5 GSPS, 50 campioni = ~10 ns.
                                  // Deve essere abbastanza piccolo da stare PRIMA dell'impulso.
                                  // Se gli impulsi arrivano molto presto nella finestra (< 10 ns),
                                  // ridurre questo valore. Se servono baselines più stabili, aumentare.

const double CFD_FRACTION = 0.5;  // Frazione del Constant Fraction Discriminator.
                                  // 0.5 = soglia al 50% dell'ampiezza dell'impulso.
                                  // Valori tipici: 0.2 (timing sul fronte iniziale, più sensibile
                                  // al rumore), 0.5 (compromesso standard), 0.8 (vicino al picco).

const double NOISE_THRESH = 5.0;  // Soglia minima di ampiezza [mV] per dichiarare che un
                                  // impulso è presente. Sotto questa soglia, il "segnale" è
                                  // compatibile con fluttuazioni di rumore. Abbassare per PMT
                                  // con segnali molto deboli; alzare se il rumore è forte.

const double CLIP_LEVEL = -499.0; // Livello di saturazione (clipping) della DRS4 [mV].
                                  // Il range di ingresso è ±500 mV; un segnale che raggiunge
                                  // -499 mV o meno è certamente saturato (es. segnali NIM a -800 mV).
const int MIN_CLIP_SAMPLES = 2;   // Numero minimo di campioni consecutivi a CLIP_LEVEL
                                  // per confermare il clipping (evita falsi positivi da
                                  // un singolo campione rumoroso).

// --- Parametri soglia fissa + correzione walk (metodo B: segnali clippati) ---
const double FIXED_THRESHOLD = -100.0; // Soglia fissa assoluta [mV].
                                      // Deve essere sopra il rumore e sotto il clipping.
                                      // Con baseline ~0 mV e rumore ~1 mV, -30 mV è
                                      // ben nel fronte di salita per impulsi > 30 mV.

const int DERIV_SMOOTH_HALFWIDTH = 2; // Semi-larghezza filtro smoothing per la derivata.
                                      // 2*N+1 = 5 campioni. A 2.5 GSPS (~400 ps/bin),
                                      // la finestra è ~2 ns, piccola vs rise time (~3 ns).

// ========================= GEOMETRIA TELESCOPI ===============================
// Posizioni dei tre telescopi nel sistema di riferimento del laboratorio.
// Origine: proiezione al suolo del centro del telescopio 08.
// Asse x: verso il corridoio (vedi planimetria)
// Asse y: perpendicolare a x nel piano orizzontale
// Asse z: verticale verso l'alto
// Unità: metri

const double POS_T08[3] = { 0.000,  0.000, 1.950};
const double POS_T06[3] = { 4.172,  3.278, 1.760};
const double POS_T04[3] = {-1.192,  6.407, 2.685};

// Vettori di baseline (differenze di posizione rispetto a T08) [m].
// Servono per l'equazione del fronte piano: d_ij · ŝ = c × Δt_ij^phys
const double BASELINE_06_08[3] = { 4.172,  3.278, -0.190};  // r_T06 - r_T08
const double BASELINE_04_08[3] = {-1.192,  6.407,  0.735};  // r_T04 - r_T08

// Velocità della luce [m/ns] — serve per convertire Δt in distanza
const double C_LIGHT = 0.299792458;

// Mappatura canali DRS4 → telescopi (0-based, corrispondente a evt.ch[])
//   CH1 (idx 0) = PMT2 del telescopio 08  → telescopio di RIFERIMENTO
//   CH2 (idx 1) = PMT2 del telescopio 06
//   CH3 (idx 2) = PMT2 del telescopio 04
//   CH4 (idx 3) = NIM trigger (non usato per la ricostruzione)
const int CH_IDX_T08 = 0;  // Indice in evt.ch[] del telescopio di riferimento
const int CH_IDX_T06 = 1;  // Indice del telescopio 06
const int CH_IDX_T04 = 2;  // Indice del telescopio 04

// Colori ROOT per i 4 canali. kBlue+1 è un blu leggermente più scuro di kBlue, ecc.
// Ordine: canale 0 (CHN1) = blu, 1 (CHN2) = verde, 2 (CHN3) = rosso, 3 (CHN4) = magenta.
const int CH_COLORS[] = {kBlue+1, kGreen+2, kRed+1, kMagenta+1};

const int BL_COLOR  = kGray+1;    // Colore della linea di baseline (grigio)
const int CFD_COLOR = kOrange+1;   // Colore dei marker CFD (arancione)
const int WALK_COLOR = kCyan+2;    // Colore dei marker soglia fissa + walk correction (ciano)


// ========================= STRUTTURE DATI ==================================
// Il programma usa due struct per organizzare i dati:
//   ChannelData: un singolo canale di un singolo evento
//   EventData:   un evento completo con tutti i suoi canali

/// ChannelData: dati grezzi + quantità estratte di UN canale di UN evento.
///
/// Le forme d'onda grezze (time[] e voltage[]) sono in float per risparmiare
/// memoria (50% rispetto a double). La precisione float è ~7 cifre significative,
/// più che sufficiente per ~200 ns di range con step da ~0.2 ns (3 cifre bastano).
///
/// Le quantità estratte (baseline, t_cfd, ...) sono in double perché il CFD
/// usa interpolazione sub-campione che richiede più precisione aritmetica.
struct ChannelData {
    int    nsamples;                  // Numero di campioni effettivamente letti (di solito 1024)
    float  time[MAX_SAMPLES];         // Tempo di ciascun campione [ns], già calibrato da DRSOsc
    float  voltage[MAX_SAMPLES];      // Tensione di ciascun campione [mV]

    // --- Quantità estratte da AnalyzeChannel() ---
    double baseline;      // Media della tensione nei primi NBL_SAMPLES campioni [mV]
    double baseline_rms;  // Deviazione standard della baseline [mV] (= rumore del canale)
    double amplitude;     // Ampiezza dell'impulso: baseline - V_min [mV] (sempre positiva)
    double v_min;         // Tensione minima (= picco negativo dell'impulso) [mV]
    double t_min;         // Tempo del campione con tensione minima [ns]
    double t_cfd;         // Tempo CFD con interpolazione sub-campione [ns]
                          //   -999.0 se non calcolato (es. segnale clippato)
    double integral;      // Integrale dell'impulso ∫(BL-V)dt [mV·ns], proporzionale alla carica

    // --- Timing soglia fissa + walk correction (per segnali clippati) ---
    double t_le;            // Tempo di attraversamento della soglia fissa [ns]
                            //   (Leading Edge, NON corretto per il walk). -999.0 se fallito.
    double slew_rate_max;   // Massima |dV/dt| sul fronte di discesa [mV/ns]
                            //   Proporzionale all'ampiezza reale del segnale.
    double walk_correction; // Correzione del walk [ns]: |V_th - BL| / slew_rate_max
    double t_walk_corrected;// t_le - walk_correction [ns]. -999.0 se fallito.

    // --- Timing unificato ---
    double t_timing;      // MIGLIOR TEMPO DISPONIBILE [ns]:
                          //   = t_cfd            se non clippato (metodo A)
                          //   = t_walk_corrected se clippato     (metodo B)
                          //   = -999.0           se fallito
    int    timing_method; // 0 = nessuno, 1 = CFD, 2 = walk corrected

    // --- Flag booleani ---
    bool   has_pulse;     // true se amplitude > NOISE_THRESH (c'è un impulso reale)
    bool   is_clipped;    // true se V_min <= CLIP_LEVEL (segnale saturato, es. NIM)
    bool   cfd_ok;        // true se il CFD ha trovato un crossing valido
};

/// EventData: header + dati completi di UN evento della DRS4.
///
/// Un evento corrisponde a un singolo trigger della DRS4. Al trigger,
/// l'onda domino si ferma e il contenuto di tutte le 1024 celle di tutti
/// i canali attivi viene letto e scritto nel file XML.
struct EventData {
    int         serial;                    // Numero progressivo nel file (parte da 1)
    std::string timestamp;                 // Data/ora dal PC: "YYYY/MM/DD HH:MM:SS.mmm"
    int         trigger_cell;              // Cella del DRS4 dove l'onda domino si è fermata (0-1023)
                                           // Nel formato XML si può ignorare (tempi già calibrati).
    int         board_serial;              // Numero seriale della scheda (es. 2856)
    int         scaler[MAX_CHANNELS];      // Rate scaler per canale [Hz] (dal discriminatore on-board)
    int         nchannels;                 // Quanti canali sono presenti in questo evento (1-4)
    int         channel_ids[MAX_CHANNELS]; // Numeri dei canali presenti (1-based: CHN1=1, CHN2=2, ...)
                                           // Esempio: se solo CHN1 e CHN3 sono attivi,
                                           //   nchannels = 2
                                           //   channel_ids[] = {1, 3, ?, ?}
                                           //   ch[0] = dati di CHN1
                                           //   ch[1] = dati di CHN3
    ChannelData ch[MAX_CHANNELS];          // Array dei dati dei canali (indicizzato 0-based)

    // --- Ricostruzione angolare dello sciame ---
    // Calcolati da ReconstructDirection() dopo l'analisi dei canali.
    // Usano le differenze temporali inter-telescopio corrette per i cable delays.
    double reco_theta;      // Angolo zenitale [rad]: 0 = verticale, π/2 = orizzonte
    double reco_phi;        // Angolo azimutale [rad]: nel sistema di rif. del laboratorio
    double reco_l;          // Coseno direttore x della PROVENIENZA: sinθ cosφ
    double reco_m;          // Coseno direttore y della PROVENIENZA: sinθ sinφ
    double reco_dt1_phys;   // Δt fisico T06-T08 dopo correzione cable [ns]
    double reco_dt2_phys;   // Δt fisico T04-T08 dopo correzione cable [ns]
    bool   reco_ok;         // true se la ricostruzione è riuscita
};


// ============================ ANALISI ======================================
// ComputeSmoothedDerivative(): calcola la derivata numerica dV/dt con filtro
// a media mobile. Invece di differenze tra campioni adiacenti (rumorose),
// calcola la differenza tra le medie di due finestre adiacenti di (hw) campioni.
//
// Per impulsi NEGATIVI, il fronte di discesa ha dV/dt < 0.
// La funzione cerca il MINIMO di dV/dt (= massimo slew rate in modulo)
// nella regione [search_start, search_end].
//
// Restituisce: indice del minimo della derivata, o -1 se non trovato.

int ComputeSmoothedDerivative(const float *t, const float *v, int ns, int hw,
                              double *dVdt, int search_start, int search_end) {

    for (int i = 0; i < ns; i++) dVdt[i] = 0.0;

    int i_first = hw;
    int i_last  = ns - hw - 1;
    if (i_first > i_last) return -1;

    for (int i = i_first; i <= i_last; i++) {
        // Media finestra "indietro": campioni [i-hw ... i-1]
        double sum_v_back = 0.0, sum_t_back = 0.0;
        for (int j = i - hw; j < i; j++) {
            sum_v_back += v[j];
            sum_t_back += t[j];
        }
        double avg_v_back = sum_v_back / hw;
        double avg_t_back = sum_t_back / hw;

        // Media finestra "avanti": campioni [i+1 ... i+hw]
        double sum_v_fwd = 0.0, sum_t_fwd = 0.0;
        for (int j = i + 1; j <= i + hw; j++) {
            sum_v_fwd += v[j];
            sum_t_fwd += t[j];
        }
        double avg_v_fwd = sum_v_fwd / hw;
        double avg_t_fwd = sum_t_fwd / hw;

        double dt_avg = avg_t_fwd - avg_t_back;
        if (fabs(dt_avg) > 1e-9)
            dVdt[i] = (avg_v_fwd - avg_v_back) / dt_avg;
    }

    // Cerca il MINIMO della derivata nella regione specificata
    int clamp_start = std::max(search_start, i_first);
    int clamp_end   = std::min(search_end, i_last);
    if (clamp_start > clamp_end) return -1;

    int    i_min_deriv = clamp_start;
    double min_deriv   = dVdt[clamp_start];
    for (int i = clamp_start + 1; i <= clamp_end; i++) {
        if (dVdt[i] < min_deriv) {
            min_deriv = dVdt[i];
            i_min_deriv = i;
        }
    }
    return i_min_deriv;
}


// AnalyzeChannel(): cuore dell'analisi fisica. Viene chiamata su OGNI canale
// di OGNI evento subito dopo il parsing.
//
// STRATEGIA DI TIMING DUALE:
//   (A) NON clippato → CFD al 35%: la soglia scala con l'ampiezza,
//       eliminando il time walk. Metodo preferito quando A è nota.
//   (B) CLIPPATO → Soglia fissa + correzione walk: quando il picco è
//       troncato, il CFD non può calcolare la soglia corretta. Si usa
//       una soglia fissa e si corregge il ritardo sistematico (walk)
//       usando lo slew rate massimo del fronte come proxy dell'ampiezza.
//
// Il campo t_timing contiene sempre il miglior tempo disponibile.

void AnalyzeChannel(ChannelData &cd) {

    // --- Inizializzazione ---
    cd.has_pulse = false;
    cd.is_clipped = false;
    cd.cfd_ok = false;
    cd.baseline = 0;
    cd.baseline_rms = 0;
    cd.amplitude = 0;
    cd.v_min = 0;
    cd.t_min = 0;
    cd.t_cfd = -999.0;
    cd.t_le = -999.0;
    cd.slew_rate_max = 0.0;
    cd.walk_correction = 0.0;
    cd.t_walk_corrected = -999.0;
    cd.t_timing = -999.0;
    cd.timing_method = 0;
    cd.integral = 0;

    int ns = cd.nsamples;
    if (ns < NBL_SAMPLES + 10) return;

    float *t = cd.time;
    float *v = cd.voltage;

    // =============================================
    //  PASSO 1: BASELINE
    // =============================================
    double sum = 0, sum2 = 0;
    for (int i = 0; i < NBL_SAMPLES; i++) {
        sum  += v[i];
        sum2 += (double)v[i] * v[i];
    }
    double bl     = sum / NBL_SAMPLES;
    double bl_rms = sqrt(fabs(sum2/NBL_SAMPLES - bl*bl));
    cd.baseline = bl;
    cd.baseline_rms = bl_rms;

    // =============================================
    //  PASSO 2: RICERCA DEL MINIMO
    // =============================================
    double vmin = v[0];
    int    imin = 0;
    for (int i = 1; i < ns; i++) {
        if (v[i] < vmin) { vmin = v[i]; imin = i; }
    }
    cd.v_min = vmin;
    cd.t_min = t[imin];

    // =============================================
    //  PASSO 3: AMPIEZZA E RILEVAMENTO CLIPPING
    // =============================================
    double amp = bl - vmin;
    cd.amplitude = amp;

    if (amp < NOISE_THRESH) return;
    cd.has_pulse = true;

    // --- Rilevamento clipping robusto ---
    // Richiede almeno MIN_CLIP_SAMPLES campioni consecutivi a CLIP_LEVEL
    // per confermare il clipping (evita falsi positivi da rumore).
    bool clipped = false;
    int  clip_count = 0;
    int  clip_start_idx = -1;

    for (int i = 0; i < ns; i++) {
        if (v[i] <= CLIP_LEVEL) {
            clip_count++;
            if (clip_count >= MIN_CLIP_SAMPLES) {
                clipped = true;
                if (clip_start_idx < 0)
                    clip_start_idx = i - MIN_CLIP_SAMPLES + 1;
            }
        } else {
            clip_count = 0;
        }
    }
    cd.is_clipped = clipped;


    // =============================================
    //  PASSO 4.1: LEADING EDGE TIME (soglia fissa) — calcolato SEMPRE
    // =============================================
    // La soglia fissa FIXED_THRESHOLD è un valore assoluto [mV].
    // Risalendo all'indietro dal minimo (o dall'inizio del clipping),
    // cerchiamo il primo crossing e interpoliamo linearmente.
    {
        int search_from = clipped ? std::max(clip_start_idx, 1) : imin;

        for (int i = search_from; i > 0; i--) {
            if (v[i] <= FIXED_THRESHOLD && v[i-1] > FIXED_THRESHOLD) {
                double dv = (double)v[i] - v[i-1];
                if (fabs(dv) > 1e-6) {
                    cd.t_le = t[i-1] + (FIXED_THRESHOLD - v[i-1]) / dv * (t[i] - t[i-1]);
                }
                break;
            }
        }
    }


    // =============================================
    //  PASSO 4.2: DERIVATA SMOOTHED E SLEW RATE MASSIMO — calcolato SEMPRE
    // =============================================
    // La derivata viene calcolata sul fronte di discesa (prima del clipping).
    // Lo slew rate massimo (= |dV/dt|_max) è proporzionale all'ampiezza reale.
    {
        double dVdt[MAX_SAMPLES];

        int deriv_search_start = NBL_SAMPLES;
        int deriv_search_end   = clipped ? std::max(clip_start_idx - 1, 0) : imin;

        int i_max_deriv = ComputeSmoothedDerivative(
            t, v, ns, DERIV_SMOOTH_HALFWIDTH,
            dVdt, deriv_search_start, deriv_search_end);

        if (i_max_deriv >= 0 && fabs(dVdt[i_max_deriv]) > 1e-6) {
            cd.slew_rate_max = fabs(dVdt[i_max_deriv]);
        }
    }


    // =============================================
    //  PASSO 4.3: CORREZIONE WALK — calcolato SEMPRE (quando possibile)
    // =============================================
    // walk_correction = |V_th - baseline| / slew_rate_max
    // Il tempo corretto: t_walk_corrected = t_le - walk_correction
    //
    // Perché sottraiamo: il walk è un RITARDO (il segnale piccolo
    // raggiunge la soglia più tardi). La correzione lo compensa.
    {
        if (cd.t_le > -900.0 && cd.slew_rate_max > 1e-3) {
            double delta_V = fabs(FIXED_THRESHOLD - bl);
            cd.walk_correction = delta_V / cd.slew_rate_max;
            cd.t_walk_corrected = cd.t_le - cd.walk_correction;
        }
    }


    // =============================================
    //  PASSO 4.4: CFD — calcolato SOLO se NON clippato
    // =============================================
    // Se il segnale è clippato, V_min è il valore di clip (non il vero picco),
    // e la soglia CFD cadrebbe nella zona saturata → risultato privo di senso.
    {
        if (!clipped) {
            double v_thr = bl + CFD_FRACTION * (vmin - bl);

            for (int i = imin; i > 0; i--) {
                if (v[i] <= v_thr && v[i-1] > v_thr) {
                    double dv = (double)v[i] - v[i-1];
                    if (fabs(dv) > 1e-6) {
                        cd.t_cfd = t[i-1] + (v_thr - v[i-1]) / dv * (t[i] - t[i-1]);
                        cd.cfd_ok = true;
                    }
                    break;
                }
            }
        }
    }


    // =============================================
    //  PASSO 4.5: TIMING UNIFICATO
    // =============================================
    // Priorità: CFD (se non clippato) > walk corrected > nessuno
    {
        if (!clipped && cd.cfd_ok) {
            cd.t_timing = cd.t_cfd;
            cd.timing_method = 1;   // CFD
        }
        else if (cd.t_walk_corrected > -900.0) {
            cd.t_timing = cd.t_walk_corrected;
            cd.timing_method = 2;   // Walk corrected
        }
    }


    // =============================================
    //  PASSO 5: INTEGRALE DELL'IMPULSO
    // =============================================
    double integ = 0;
    double thr = bl - 3.0 * bl_rms;
    for (int i = 0; i < ns-1; i++) {
        if (v[i] < thr)
            integ += (bl - v[i]) * (t[i+1] - t[i]);
    }
    cd.integral = integ;
}


// ==================== CABLE DELAY OFFSETS (modificabili) ======================
// I ritardi dei cavi [ns] determinano l'offset sistematico tra i tempi misurati
// dal DRS4 e i tempi fisici di arrivo delle particelle dello sciame.
//
// Δt_phys(T_j - T_ref) = Δt_misurato(CH_j - CH_ref) - cable_offset_j
//
// dove cable_offset_j = delay_j - delay_ref.
//
// Valori iniziali: misurati con il metodo dei rimbalzi all'oscilloscopio.
//   T08: 17 ns,  T06: 205 ns,  T04: 67.5 ns
//
// Questi possono essere aggiornati con SetCableOffsets() o CalibrateCableOffsets().
// NOTA: i cable delays includono solo la propagazione nei cavi; i ritardi
// aggiuntivi (PMT transit time, amplificatore, ecc.) NON sono inclusi.
// Per una calibrazione completa, usare CalibrateCableOffsets() che stima
// gli offset totali dalla media dei Δt misurati.

static double gCableOffset_06_08 = 205.0 - 17.0 ;   // delay_T06 - delay_T08 = 205 - 17 [ns]
static double gCableOffset_04_08 =  67.5 - 17.0 ;   // delay_T04 - delay_T08 = 67.5 - 17 [ns]


// ==================== RICOSTRUZIONE ANGOLARE ==================================
// ReconstructDirection(): dato un evento con 3 timing inter-telescopio validi,
// ricostruisce la direzione di provenienza dello sciame usando l'approssimazione
// di fronte d'onda piano.
//
// MODELLO FISICO:
// ---------------
// Lo sciame esteso produce un fronte di particelle approssimabile come un piano
// che si propaga alla velocità c nella direzione n̂ (verso il basso).
// L'istante di arrivo al telescopio i è:
//
//   t_i = t_0 - r_i · n̂ / c
//
// La direzione di PROVENIENZA dello sciame è ŝ = -n̂ (verso l'alto).
// In termini di ŝ, le differenze temporali fisiche sono:
//
//   Δt_ij^phys = t_j - t_i = (r_j - r_i) · ŝ / c = d_ij · ŝ / c
//
// dove d_ij = r_j - r_i è il vettore di baseline.
//
// METODO DI SOLUZIONE:
// --------------------
// Con 3 telescopi, abbiamo 2 equazioni indipendenti in 3 incognite (ŝ):
//
//   d_1 · ŝ = c × Δt_1    (baseline T06-T08)
//   d_2 · ŝ = c × Δt_2    (baseline T04-T08)
//
// più il vincolo |ŝ| = 1.
//
// Il sistema lineare M·ŝ = b (M è 2×3) ha uno spazio di soluzioni 1D
// (una retta nel 3D). L'intersezione con la sfera unitaria dà 0 o 2 punti.
//
// Soluzione:
//   1. ŝ₀ = Mᵀ(MMᵀ)⁻¹ b    (pseudoinversa: soluzione a norma minima)
//   2. k̂ = (d_1 × d_2)/|d_1 × d_2|   (direzione dello spazio nullo di M)
//   3. ŝ = ŝ₀ + λ k̂   con  |ŝ₀|² + λ² = 1  (perché ŝ₀ ⊥ k̂)
//   4. λ = ±√(1 - |ŝ₀|²),  si sceglie ŝ con s_z > 0 (viene dall'alto)
//
// Se |ŝ₀|² > 1, i Δt misurati sono incompatibili con un fronte piano a
// velocità c (segnale "superluminale") → ricostruzione fallita.
//
// PARAMETRI E RESTITUZIONI:
//   dt1_phys, dt2_phys: differenze temporali fisiche [ns] (corrette per i cavi)
//   theta: angolo zenitale [rad]  (output)
//   phi:   angolo azimutale [rad] (output)
//   l, m:  coseni direttori orizzontali della provenienza (output)
//   Restituisce: true se la ricostruzione è riuscita
//
// ===========================================================================

bool ReconstructDirection(double dt1_phys, double dt2_phys,
                          double &theta, double &phi,
                          double &l_out, double &m_out) {

    // Alias per le baseline [m]
    const double *d1 = BASELINE_06_08;
    const double *d2 = BASELINE_04_08;

    // --- Prodotti scalari delle baseline ---
    double d1d1 = d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2];  // |d1|²
    double d2d2 = d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2];  // |d2|²
    double d1d2 = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2];  // d1 · d2

    // --- Determinante di M·Mᵀ (matrice 2×2) ---
    //   M·Mᵀ = |d1·d1  d1·d2|
    //           |d1·d2  d2·d2|
    double det = d1d1 * d2d2 - d1d2 * d1d2;

    if (fabs(det) < 1e-12) {
        // Le baseline sono parallele → ricostruzione impossibile
        // (non dovrebbe succedere con la nostra geometria)
        return false;
    }

    // --- Vettore dei termini noti: b = c × [Δt1, Δt2] ---
    double b1 = C_LIGHT * dt1_phys;   // [m]
    double b2 = C_LIGHT * dt2_phys;   // [m]

    // --- Inversa di M·Mᵀ ---
    double inv00 =  d2d2 / det;
    double inv01 = -d1d2 / det;
    double inv10 = -d1d2 / det;
    double inv11 =  d1d1 / det;

    // --- Moltiplicazione: (MMᵀ)⁻¹ × b → vettore 2D alpha ---
    double alpha0 = inv00 * b1 + inv01 * b2;
    double alpha1 = inv10 * b1 + inv11 * b2;

    // --- Pseudoinversa: ŝ₀ = Mᵀ × alpha → vettore 3D ---
    //   Mᵀ = |d1[0] d2[0]|
    //         |d1[1] d2[1]|     ŝ₀ = Mᵀ × alpha = d1 × alpha0 + d2 × alpha1
    //         |d1[2] d2[2]|
    double s0[3];
    for (int k = 0; k < 3; k++)
        s0[k] = d1[k] * alpha0 + d2[k] * alpha1;

    // --- Norma quadra di ŝ₀ ---
    double s0_sq = s0[0]*s0[0] + s0[1]*s0[1] + s0[2]*s0[2];

    if (s0_sq > 1.0) {
        // I Δt misurati sono incompatibili con un fronte piano a velocità c.
        // Questo succede quando i ritardi sono troppo grandi rispetto alle
        // baseline (segnale "superluminale"): |Δt| > |d|/c.
        // Possibili cause: errore nella calibrazione dei cable delays,
        // rumore sul timing, o evento non fisico (coincidenza casuale).
        return false;
    }

    // --- Vettore del kernel: k = d1 × d2 (prodotto vettoriale), normalizzato ---
    double kk[3];
    kk[0] = d1[1]*d2[2] - d1[2]*d2[1];   //  3.278×0.735 - (-0.19)×6.407
    kk[1] = d1[2]*d2[0] - d1[0]*d2[2];   // (-0.19)×(-1.192) - 4.172×0.735
    kk[2] = d1[0]*d2[1] - d1[1]*d2[0];   //  4.172×6.407 - 3.278×(-1.192)

    double kk_norm = sqrt(kk[0]*kk[0] + kk[1]*kk[1] + kk[2]*kk[2]);
    if (kk_norm < 1e-12) return false;

    for (int k = 0; k < 3; k++) kk[k] /= kk_norm;  // Normalizza

    // --- Lambda dalla condizione |ŝ|² = 1 ---
    // |ŝ₀ + λk̂|² = |ŝ₀|² + λ² = 1   (perché ŝ₀ ⊥ k̂ per costruzione)
    double lambda_sq = 1.0 - s0_sq;
    double lambda = sqrt(lambda_sq);   // Sempre ≥ 0 (già controllato s0_sq ≤ 1)

    // --- Due soluzioni: ŝ = ŝ₀ ± λ k̂ ---
    // Scegliamo quella con s_z > 0 (la provenienza è dall'alto).
    double s_plus[3], s_minus[3];
    for (int k = 0; k < 3; k++) {
        s_plus[k]  = s0[k] + lambda * kk[k];
        s_minus[k] = s0[k] - lambda * kk[k];
    }

    // Scegli la soluzione con s_z > 0
    double *s_best = nullptr;
    if (s_plus[2] > 0 && s_minus[2] > 0) {
        // Entrambe puntano verso l'alto: prendi quella più verticale
        s_best = (s_plus[2] > s_minus[2]) ? s_plus : s_minus;
    } else if (s_plus[2] > 0) {
        s_best = s_plus;
    } else if (s_minus[2] > 0) {
        s_best = s_minus;
    } else {
        // Nessuna soluzione con s_z > 0 (lo sciame verrebbe dal basso)
        return false;
    }

    // --- Estrazione degli angoli ---
    // θ = angolo zenitale: cos θ = s_z  →  θ = arccos(s_z)
    // φ = angolo azimutale: φ = atan2(s_y, s_x)
    // l = s_x = sinθ cosφ (coseno direttore x della provenienza)
    // m = s_y = sinθ sinφ (coseno direttore y della provenienza)

    double sz = s_best[2];
    if (sz > 1.0) sz = 1.0;     // Protezione numerica
    if (sz < 0.0) sz = 0.0;

    theta = acos(sz);              // [rad], 0 = verticale
    phi   = atan2(s_best[1], s_best[0]);  // [rad], [-π, π]
    if (phi < 0) phi += 2.0 * M_PI;      // Converte in [0, 2π]

    l_out = s_best[0];   // sinθ cosφ
    m_out = s_best[1];   // sinθ sinφ

    return true;
}


// ==================== RICOSTRUZIONE PER UN EVENTO ============================
// ReconstructEvent(): applica la ricostruzione angolare a un singolo evento.
// Controlla che i 3 canali PMT abbiano timing valido, sottrae i cable offsets,
// e chiama ReconstructDirection().

void ReconstructEvent(EventData &evt) {
    // Inizializzazione
    evt.reco_ok = false;
    evt.reco_theta = 0;
    evt.reco_phi = 0;
    evt.reco_l = 0;
    evt.reco_m = 0;
    evt.reco_dt1_phys = 0;
    evt.reco_dt2_phys = 0;

    // Servono almeno 3 canali (i primi 3 sono i PMT)
    if (evt.nchannels < 3) return;

    // Tutti e 3 i canali PMT devono avere timing valido
    if (evt.ch[CH_IDX_T08].t_timing < -900.0 ||
        evt.ch[CH_IDX_T06].t_timing < -900.0 ||
        evt.ch[CH_IDX_T04].t_timing < -900.0) return;

    // --- Differenze temporali misurate ---
    double dt1_meas = evt.ch[CH_IDX_T06].t_timing - evt.ch[CH_IDX_T08].t_timing;
    double dt2_meas = evt.ch[CH_IDX_T04].t_timing - evt.ch[CH_IDX_T08].t_timing;

    // --- Correzione cable delays ---
    // Δt_phys = Δt_meas - cable_offset
    evt.reco_dt1_phys = dt1_meas - gCableOffset_06_08;
    evt.reco_dt2_phys = dt2_meas - gCableOffset_04_08;

    // --- Ricostruzione ---
    evt.reco_ok = ReconstructDirection(
        evt.reco_dt1_phys, evt.reco_dt2_phys,
        evt.reco_theta, evt.reco_phi,
        evt.reco_l, evt.reco_m
    );
}



// ParseXML(): legge un file XML della DRS4 e popola il vettore di eventi.
//
// Il parser usa una MACCHINA A STATI FINITI con 4 stati:
//   IDLE       → attesa di un nuovo <Event>
//   IN_EVENT   → lettura header (Serial, Time, HUnit, VUnit)
//   IN_BOARD   → lettura dati scheda (Trigger_Cell, Scaler, apertura canali)
//   IN_CHANNEL → lettura campioni <Data>tempo,tensione</Data>
//
// Ogni riga del file viene letta con getline(), copiata in un buffer C,
// e analizzata con sscanf() che è molto più veloce di un parser XML generico.
//
// Restituisce: il numero di eventi letti con successo, o -1 in caso di errore.

int ParseXML(const char* filename, std::vector<EventData> &events) {

    // Apri il file in lettura
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "[ERRORE] Impossibile aprire: " << filename << std::endl;
        return -1;
    }
    std::cout << "[INFO] Parsing " << filename << " ..." << std::flush;

    // Definizione degli stati della macchina
    enum State { IDLE, IN_EVENT, IN_BOARD, IN_CHANNEL };
    State state = IDLE;                 // Stato corrente: partiamo in IDLE

    EventData current_evt;              // Evento in costruzione
    int ch_idx = -1;                    // Indice del canale corrente in ch[] (0-based)
    int sample_count = 0;               // Contatore campioni nel canale corrente
    int scaler_idx = 0;                 // Contatore scaler (per mappatura progressiva)
    std::string line;                   // Riga corrente del file
    char buf[512];                      // Buffer C per sscanf (più veloce di stringstream)

    // Lettura riga per riga — memoria costante O(1), non carica tutto in RAM
    while (std::getline(infile, line)) {

        // Rimuovi '\r' se il file ha line ending Windows (\r\n)
        if (!line.empty() && line.back() == '\r') line.pop_back();

        // Copia la stringa C++ in un buffer C per sscanf
        strncpy(buf, line.c_str(), sizeof(buf)-1);
        buf[sizeof(buf)-1] = '\0';   // Assicura terminazione

        // ========================
        //  STATO: IDLE
        // ========================
        // In attesa di un tag <Event>. Tutto il resto viene ignorato
        // (commenti XML, header <?xml ...?>, tag <DRSOSC>, ecc.)
        if (state == IDLE) {
            if (line.find("<Event>") != std::string::npos) {
                // Inizio di un nuovo evento: resetta la struttura
                current_evt = EventData();      // Costruttore default (azzera tutto)
                current_evt.nchannels = 0;
                memset(current_evt.scaler, 0, sizeof(current_evt.scaler));
                scaler_idx = 0;
                state = IN_EVENT;               // Transizione: IDLE → IN_EVENT
            }
            continue;   // Passa alla riga successiva
        }

        // ========================
        //  STATO: IN_EVENT
        // ========================
        // Legge i tag dell'header: Serial, Time, HUnit (ignorato), VUnit (ignorato).
        // Quando trova <Board_NNNN>, transita a IN_BOARD.
        // Quando trova </Event>, analizza e salva l'evento.
        if (state == IN_EVENT) {
            int val;
            char ts[64];    // Buffer per il timestamp

            // <Serial>123</Serial>
            if (sscanf(buf, " <Serial>%d</Serial>", &val) == 1) {
                current_evt.serial = val;
            }
            // <Time>2026/03/17 16:39:01.576</Time>
            // %63[^<] = leggi fino a 63 caratteri qualsiasi tranne '<'
            else if (sscanf(buf, " <Time>%63[^<]</Time>", ts) == 1) {
                current_evt.timestamp = ts;
            }
            // <Board_2856> → ingresso nel blocco della scheda
            else if (sscanf(buf, " <Board_%d>", &val) == 1) {
                current_evt.board_serial = val;
                state = IN_BOARD;               // Transizione: IN_EVENT → IN_BOARD
            }
            // </Event> → evento completo!
            else if (line.find("</Event>") != std::string::npos) {
                // Analizza TUTTI i canali di questo evento
                for (int i = 0; i < current_evt.nchannels-1; i++) {
                    AnalyzeChannel(current_evt.ch[i]);
                }
                // Ricostruzione angolare (se 3 canali PMT hanno timing valido)
                ReconstructEvent(current_evt);
                // Aggiungi l'evento al vettore (copia l'intera struttura)
                events.push_back(current_evt);
                state = IDLE;                   // Transizione: IN_EVENT → IDLE
            }
            // HUnit, VUnit: ignorati (sono sempre "ns" e "mV")
            continue;
        }

        // ========================
        //  STATO: IN_BOARD
        // ========================
        // Legge Trigger_Cell, Scaler, e i tag di apertura dei canali <CHNn>.
        if (state == IN_BOARD) {
            int val, val2;

            // <Trigger_Cell>735</Trigger_Cell>
            if (sscanf(buf, " <Trigger_Cell>%d</Trigger_Cell>", &val) == 1) {
                current_evt.trigger_cell = val;
            }
            // <Scaler0>0</Scaler0> oppure <Scaler1>1234</Scaler1>
            // Il formato "%d>%d" legge sia l'indice dello scaler che il valore.
            // Usiamo scaler_idx (contatore progressivo) anziché il numero nel tag
            // per gestire sia la numerazione 0-based che 1-based.
            else if (sscanf(buf, " <Scaler%d>%d</Scaler", &val, &val2) == 2) {
                if (scaler_idx < MAX_CHANNELS) {
                    current_evt.scaler[scaler_idx] = val2;
                }
                scaler_idx++;
            }
            // <CHN1> → inizio dati del canale 1
            else if (sscanf(buf, " <CHN%d>", &val) == 1) {
                int idx = current_evt.nchannels;  // Prossimo slot libero
                if (idx < MAX_CHANNELS) {
                    current_evt.channel_ids[idx] = val;    // Salva il numero CHN (1-based)
                    current_evt.ch[idx].nsamples = 0;
                    ch_idx = idx;                          // Indice corrente per i campioni
                    current_evt.nchannels++;               // Un canale in più
                    sample_count = 0;
                    state = IN_CHANNEL;                    // Transizione: IN_BOARD → IN_CHANNEL
                }
            }
            // </Board_2856> → torna a IN_EVENT (potrebbe esserci un secondo board,
            // ma nel nostro caso c'è sempre una sola scheda)
            else if (line.find("</Board_") != std::string::npos) {
                state = IN_EVENT;                          // Transizione: IN_BOARD → IN_EVENT
            }
            continue;
        }

        // ========================
        //  STATO: IN_CHANNEL
        // ========================
        // Legge le coppie <Data>tempo,tensione</Data> una per riga.
        // Quando trova </CHNn>, torna a IN_BOARD per il canale successivo.
        if (state == IN_CHANNEL) {
            // </CHN1> → fine del canale corrente
            if (line.find("</CHN") != std::string::npos) {
                if (ch_idx >= 0) {
                    current_evt.ch[ch_idx].nsamples = sample_count;  // Registra quanti campioni letti
                }
                state = IN_BOARD;   // Transizione: IN_CHANNEL → IN_BOARD
                continue;
            }

            // <Data>12.500,-35.4</Data>
            // sscanf con "%f,%f" legge due float separati da virgola
            float tv, vv;
            if (ch_idx >= 0 && sample_count < MAX_SAMPLES &&
                sscanf(buf, " <Data>%f,%f</Data>", &tv, &vv) == 2) {
                current_evt.ch[ch_idx].time[sample_count]    = tv;  // Tempo [ns]
                current_evt.ch[ch_idx].voltage[sample_count] = vv;  // Tensione [mV]
                sample_count++;
            }
            continue;
        }
    }

    infile.close();
    std::cout << " " << events.size() << " eventi." << std::endl;
    return events.size();
}


// ===================== VARIABILI GLOBALI ===================================
// Queste variabili sono globali perché devono essere accessibili da tutte
// le funzioni di navigazione (Next, Prev, GoTo, Overlay, ...) che l'utente
// chiama dal prompt ROOT. Il prefisso 'g' le distingue dalle variabili locali.

static std::vector<EventData> gEvents;           // Vettore con TUTTI gli eventi del file
static int gCurrentEvent = 0;                    // Indice dell'evento attualmente visualizzato
static std::string gFilename;                    // Nome del file XML (per i titoli)
static std::map<int, std::string> gChannelLabels;// Mappa CHN_number → etichetta personalizzata
static TCanvas *gCanvas = nullptr;               // Canvas principale del browser
static TPad *gPadWF[MAX_CHANNELS] = {};          // Pad per le forme d'onda (1 per canale)
static TPad *gPadInfo = nullptr;                 // Pad del pannello informativo in basso


// ===================== DRAW EVENT ==========================================
// DrawEvent(): disegna un evento nel canvas multi-pad del browser.
// Viene chiamata da Next(), Prev(), GoTo() e all'avvio.
//
// Per ogni canale:
//  1. Crea un TGraph con i campioni (tempo, tensione)
//  2. Disegna la linea di baseline tratteggiata
//  3. Se c'è un impulso: disegna il marker CFD (linea verticale + orizzontale)
//  4. Aggiunge annotazioni (ampiezza, t_CFD, carica, flag CLIPPED)
//  5. Aggiorna il pannello info in basso con header e Δt tra canali

void DrawEvent(int ievt) {

    // Controllo limiti
    if (ievt < 0 || ievt >= (int)gEvents.size()) {
        std::cout << "[WARNING] Indice " << ievt << " fuori range [0, "
                  << gEvents.size()-1 << "]" << std::endl;
        return;
    }

    gCurrentEvent = ievt;                   // Aggiorna l'indice globale
    EventData &evt = gEvents[ievt];         // Riferimento all'evento corrente

    // === Loop sui canali dell'evento ===
    for (int i = 0; i < evt.nchannels; i++) {
        if (!gPadWF[i]) continue;           // Safety: skip se il pad non esiste

        gPadWF[i]->cd();                   // Rendi questo pad il "canvas corrente" per ROOT
        gPadWF[i]->Clear();                // Cancella tutto il contenuto precedente
                                            // (ROOT distrugge automaticamente gli oggetti
                                            //  creati con 'new' che erano nel pad)

        ChannelData &cd = evt.ch[i];        // Riferimento ai dati del canale
        int ch_id = evt.channel_ids[i];     // Numero del canale (1-based: 1, 2, 3, o 4)

        // --- Crea il TGraph della forma d'onda ---
        // TGraph(n, x[], y[]) con 'new' → ROOT ne gestisce la ownership.
        // Verrà distrutto dal Clear() alla prossima chiamata di DrawEvent.
        TGraph *gr = new TGraph();
        for (int s = 0; s < cd.nsamples; s++)
            gr->SetPoint(s, cd.time[s], cd.voltage[s]);   // Aggiungi punto per punto

        // Etichetta canale: usa quella personalizzata se esiste, altrimenti "CHNn"
        std::string label = gChannelLabels.count(ch_id)
                            ? gChannelLabels[ch_id]         // Es. "PMT top"
                            : Form("CHN%d", ch_id);         // Es. "CHN1"

        // Titolo del pad: "PMT top | Ev. 1 (#1/20);asse X;asse Y"
        // Il formato ROOT usa ';' per separare titolo, asse X, asse Y
        gr->SetTitle(Form("%s  |  Ev. %d (#%d/%d);t [ns];V [mV]",
                          label.c_str(), evt.serial,
                          ievt+1, (int)gEvents.size()));
        gr->SetLineColor(CH_COLORS[i % 4]);   // Colore del canale
        gr->SetLineWidth(1);                   // Linea sottile per la forma d'onda

        // --- Range Y "intelligente" ---
        // Cerca il minimo e massimo della forma d'onda per impostare il range
        // dell'asse Y in modo che l'impulso sia ben visibile
        double vmin_d = cd.voltage[0], vmax_d = cd.voltage[0];
        for (int s = 1; s < cd.nsamples; s++) {
            if (cd.voltage[s] < vmin_d) vmin_d = cd.voltage[s];
            if (cd.voltage[s] > vmax_d) vmax_d = cd.voltage[s];
        }
        double vrange = vmax_d - vmin_d;
        if (vrange < 5) vrange = 5;    // Protezione: se il segnale è piatto, evita range=0
        // Margini: 10% sotto il minimo, 18% sopra il massimo (spazio per annotazioni)
        double ylo = vmin_d - 0.10 * vrange;
        double yhi = vmax_d + 0.18 * vrange;

        gr->Draw("AL");    // "A" = disegna gli assi, "L" = collega i punti con linee
        gr->GetYaxis()->SetRangeUser(ylo, yhi);  // Imposta il range Y calcolato

        // --- Linea tratteggiata della baseline ---
        if (cd.nsamples > NBL_SAMPLES) {
            double tmax = cd.time[cd.nsamples - 1];   // Tempo dell'ultimo campione
            TLine *lbl = new TLine(0, cd.baseline, tmax, cd.baseline);
            lbl->SetLineColor(BL_COLOR);      // Grigio
            lbl->SetLineStyle(7);             // Tratteggiata lunga
            lbl->SetLineWidth(1);
            lbl->Draw("same");                // "same" = disegna nel pad corrente
        }

        // --- Marker CFD (arancione, solo se CFD ok e non clippato) ---
        if (cd.cfd_ok && !cd.is_clipped) {
            double v_thr = cd.baseline + CFD_FRACTION * (cd.v_min - cd.baseline);

            // Linea verticale al tempo CFD
            TLine *lcfd = new TLine(cd.t_cfd, ylo*0.95, cd.t_cfd, cd.baseline);
            lcfd->SetLineColor(CFD_COLOR);
            lcfd->SetLineStyle(2);
            lcfd->SetLineWidth(2);
            lcfd->Draw("same");

            // Linea orizzontale alla soglia CFD
            double dt_vis = cd.time[cd.nsamples-1] * 0.04;
            TLine *lhr = new TLine(cd.t_cfd - dt_vis, v_thr,
                                   cd.t_cfd + dt_vis, v_thr);
            lhr->SetLineColor(CFD_COLOR);
            lhr->SetLineWidth(2);
            lhr->Draw("same");
        }

        // --- Marker soglia fissa + walk correction (ciano, solo se clippato e walk ok) ---
        if (cd.is_clipped && cd.t_walk_corrected > -900.0) {
            // Linea verticale tratteggiata al tempo walk-corretto
            TLine *lwc = new TLine(cd.t_walk_corrected, ylo*0.95,
                                   cd.t_walk_corrected, cd.baseline);
            lwc->SetLineColor(WALK_COLOR);
            lwc->SetLineStyle(2);
            lwc->SetLineWidth(2);
            lwc->Draw("same");

            // Linea orizzontale alla soglia fissa
            double dt_vis = cd.time[cd.nsamples-1] * 0.04;
            TLine *lhr_le = new TLine(cd.t_walk_corrected - dt_vis, FIXED_THRESHOLD,
                                      cd.t_walk_corrected + dt_vis, FIXED_THRESHOLD);
            lhr_le->SetLineColor(WALK_COLOR);
            lhr_le->SetLineWidth(2);
            lhr_le->Draw("same");

            // Linea sottile al tempo LE NON corretto (per vedere la correzione)
            if (cd.t_le > -900.0) {
                TLine *lle = new TLine(cd.t_le, FIXED_THRESHOLD - 15,
                                       cd.t_le, FIXED_THRESHOLD + 15);
                lle->SetLineColor(WALK_COLOR);
                lle->SetLineStyle(3);
                lle->SetLineWidth(1);
                lle->Draw("same");
            }
        }

        // --- Annotazioni di testo nel pad ---
        // NDC = Normalized Device Coordinates: (0,0) = angolo in basso a sinistra,
        // (1,1) = angolo in alto a destra del pad
        TLatex *ltx = new TLatex();
        ltx->SetNDC();
        ltx->SetTextFont(42);     // Font 42 = Helvetica medium

        if (cd.has_pulse) {
            // ltx->SetTextSize(0.050);
            // ltx->SetTextColor(CH_COLORS[i % 4]);
            // ltx->DrawLatex(0.13, 0.87, Form("A = %.1f mV", cd.amplitude));

            // Mostra il tempo e il metodo usato
            if (cd.timing_method == 1) {
                // CFD (arancione)
                ltx->SetTextColor(CFD_COLOR);
                ltx->DrawLatex(0.13, 0.79, Form("t_{CFD} = %.2f ns", cd.t_cfd));
            } else if (cd.timing_method == 2) {
                // Walk corrected (ciano)
                ltx->SetTextColor(WALK_COLOR);
                ltx->DrawLatex(0.13, 0.79, Form("t_{WC} = %.2f ns", cd.t_walk_corrected));
                ltx->SetTextSize(0.038);
                ltx->DrawLatex(0.13, 0.72, Form("SR=%.0f mV/ns  #Delta_{w}=%.2f ns",
                               cd.slew_rate_max, cd.walk_correction));
            }

            ltx->SetTextSize(0.050);
            ltx->SetTextColor(CH_COLORS[i % 4]);
            // // #upoint = unicode ·  (punto centrato, per l'unità mV·ns)
            // ltx->DrawLatex(0.13, cd.timing_method == 2 ? 0.64 : 0.71,
            //                Form("Q = %.0f mV#upointns", cd.integral));

            // Indicatore rosso CLIPPED (se il segnale ha saturato la DRS4)
            if (cd.is_clipped) {
                ltx->SetTextSize(0.055);
                ltx->SetTextColor(kRed);
                ltx->DrawLatex(0.62, 0.87, "CLIPPED");
                // Mostra il metodo usato
                // ltx->SetTextSize(0.040);
                // ltx->SetTextColor(WALK_COLOR);
                // ltx->DrawLatex(0.62, 0.79, "Walk corr.");
            }
        } else {
            // ltx->SetTextSize(0.050);
            // ltx->SetTextColor(kGray+1);
            // ltx->DrawLatex(0.40, 0.50, "No pulse detected");
        }

        // // Baseline info (piccolo, in basso a destra)
        // ltx->SetTextSize(0.038);
        // ltx->SetTextColor(BL_COLOR);
        // ltx->DrawLatex(0.55, 0.20, Form("BL = %.1f #pm %.1f mV",
        //                                  cd.baseline, cd.baseline_rms));

        // Forza il ridisegno del pad
        gPadWF[i]->Modified();
        gPadWF[i]->Update();
    }

    // === Pannello informativo in basso ===
    gPadInfo->cd();
    gPadInfo->Clear();

    // TPaveText: box di testo con sfondo e bordo
    TPaveText *info = new TPaveText(0.01, 0.05, 0.99, 0.95, "NDC");
    info->SetFillColor(kWhite);
    info->SetBorderSize(0);
    info->SetTextAlign(12);    // Allineamento: sinistra, centrato verticalmente
    info->SetTextFont(42);
    info->SetTextSize(0.18);   // Dimensione relativa al pad (che è piccolo)

    // // Riga 1: informazioni dell'header
    // info->AddText(Form("Evento %d/%d   Serial: %d   %s   "
    //                    "Board: %d   TrigCell: %d",
    //                    ievt+1, (int)gEvents.size(),
    //                    evt.serial, evt.timestamp.c_str(),
    //                    evt.board_serial, evt.trigger_cell));
    
    // Riga 1: informazioni dell'header ridotta a quello che ci interessa
    info->AddText(Form("Evento %d/%d",
                       ievt+1, (int)gEvents.size()));
                       
                                       

    // Riga 2: differenze temporali usando t_timing (funziona anche per clippati)
    std::string dt_str = "#Deltat: ";
    bool has_dt = false;
    for (int i = 0; i < evt.nchannels; i++) {
        for (int j = i+1; j < evt.nchannels; j++) {
            // Salta CH4 (trigger NIM, nessun significato fisico nel dt)
            if (evt.channel_ids[i] == 4 || evt.channel_ids[j] == 4) continue;
            if (evt.ch[i].t_timing > -900.0 && evt.ch[j].t_timing > -900.0) {
                double dt = evt.ch[j].t_timing - evt.ch[i].t_timing;
                // Indica i metodi usati: C=CFD, W=Walk
                const char* m_i = evt.ch[i].timing_method == 1 ? "C" : "W";
                const char* m_j = evt.ch[j].timing_method == 1 ? "C" : "W";
                dt_str += Form("(%d-%d)=%+.2f ns [%s,%s]   ",
                               evt.channel_ids[j], evt.channel_ids[i], dt,
                               m_j, m_i);
                has_dt = true;
            }
        }
    }
    if (!has_dt) dt_str += "(N/A)";
    info->AddText(dt_str.c_str());
    info->AddText("\n");

    // Riga 3: risultato della ricostruzione angolare
    if (evt.reco_ok) {
        double theta_deg = evt.reco_theta * 180.0 / M_PI;
        double phi_deg   = evt.reco_phi   * 180.0 / M_PI;
        info->AddText(Form("Direzione: #theta = %.1f#circ   #phi = %.1f#circ   "
                           "(l=%.3f, m=%.3f)   "
                           "#Deltat_{phys}: (06-08)=%+.1f ns  (04-08)=%+.1f ns",
                           theta_deg, phi_deg,
                           evt.reco_l, evt.reco_m,
                           evt.reco_dt1_phys, evt.reco_dt2_phys));
    } else {
        info->AddText("\nDirezione: ricostruzione non riuscita \n");
    }

    // // Riga 4: riepilogo comandi
    // info->AddText("Comandi: Next()  Prev()  GoTo(n)  Overlay()  "
    //               "OverlaySelect(\"1,2,3\")  Summary()  Skymap()");
    
    // Riga 4: riepilogo comandi semplificato
    TText *legenda_terminale = info->AddText("Comandi nel terminale");
    legenda_terminale->SetTextFont(42);    // solo questa riga
    legenda_terminale->SetTextSize(0.15);  // solo questa riga
    legenda_terminale->SetTextAlign(32); // solo questa riga
    //info->AddText("Comandi nel terminale");              


    info->Draw();
    gPadInfo->Modified();
    gPadInfo->Update();

    // Aggiorna l'intero canvas
    gCanvas->Modified();
    gCanvas->Update();
}


// ================ FUNZIONI DI NAVIGAZIONE ==================================
// Funzioni globali chiamabili dal prompt ROOT. Devono essere globali (non
// dentro una classe) perché ROOT le riconosce automaticamente come comandi.

void Next()  { DrawEvent(gCurrentEvent + 1); }   // Evento successivo
void Prev()  { DrawEvent(gCurrentEvent - 1); }   // Evento precedente
void GoTo(int n) { DrawEvent(n); }                // Vai all'evento n (0-based)

/// GoToSerial: cerca l'evento con un dato numero seriale (diverso dall'indice!)
void GoToSerial(int serial) {
    for (int i = 0; i < (int)gEvents.size(); i++)
        if (gEvents[i].serial == serial) { DrawEvent(i); return; }
    std::cout << "[WARNING] Serial " << serial << " non trovato." << std::endl;
}

/// Save: salva il canvas corrente come immagine
void Save(const char* f = "") {
    if (strlen(f) == 0)
        gCanvas->SaveAs(Form("DRS4_ev%d.png", gEvents[gCurrentEvent].serial));
    else
        gCanvas->SaveAs(f);   // ROOT deduce il formato dall'estensione (.png, .pdf, .eps, ...)
}

/// SaveAll: salva TUTTI gli eventi in un PDF multipagina
void SaveAll(const char* pdf = "DRS4_all_events.pdf") {
    std::cout << "[INFO] Salvataggio " << gEvents.size() << " eventi ..." << std::endl;
    for (int i = 0; i < (int)gEvents.size(); i++) {
        DrawEvent(i);
        // Convenzione ROOT per PDF multipagina:
        //   "file.pdf(" = prima pagina (apre il file)
        //   "file.pdf"  = pagine intermedie
        //   "file.pdf)" = ultima pagina (chiude il file)
        if (i == 0) gCanvas->Print(Form("%s(", pdf));
        else if (i == (int)gEvents.size()-1) gCanvas->Print(Form("%s)", pdf));
        else gCanvas->Print(pdf);
    }
    std::cout << "[INFO] Salvato: " << pdf << std::endl;
}

/// Summary: stampa media e sigma delle differenze temporali
/// tra i primi 3 canali, usando t_timing (CFD o walk-corrected).
void Summary() {
    int ntot = gEvents.size(), ngood = 0;
    int n_allCFD = 0, n_mixed = 0, n_allWalk = 0;
    std::vector<double> d12, d13, d23;

    for (auto &e : gEvents) {
        if (e.nchannels < 3) continue;
        bool ok = true;
        for (int i = 0; i < 3; i++)
            if (e.ch[i].t_timing < -900.0) ok = false;
        if (!ok) continue;
        ngood++;

        // Conta i tipi di combinazioni
        bool all_cfd = true, all_walk = true;
        for (int i = 0; i < 3; i++) {
            if (e.ch[i].timing_method != 1) all_cfd = false;
            if (e.ch[i].timing_method != 2) all_walk = false;
        }
        if (all_cfd) n_allCFD++;
        else if (all_walk) n_allWalk++;
        else n_mixed++;

        d12.push_back(e.ch[1].t_timing - e.ch[0].t_timing);
        d13.push_back(e.ch[2].t_timing - e.ch[0].t_timing);
        d23.push_back(e.ch[2].t_timing - e.ch[1].t_timing);
    }

    std::cout << "\n=== SOMMARIO (timing unificato) ===" << std::endl;
    std::cout << "  File: " << gFilename << std::endl;
    std::cout << "  Eventi: " << ntot << "  (con 3 timing buoni: " << ngood << ")" << std::endl;
    std::cout << "    tutti CFD: " << n_allCFD
              << "  misti: " << n_mixed
              << "  tutti Walk: " << n_allWalk << std::endl;

    // Statistiche clipping per canale
    for (int ch = 0; ch < std::min(3, (int)gEvents[0].nchannels); ch++) {
        int n_clip = 0, n_cfd = 0, n_walk = 0;
        for (auto &e : gEvents) {
            if (ch < e.nchannels) {
                if (e.ch[ch].is_clipped) n_clip++;
                if (e.ch[ch].timing_method == 1) n_cfd++;
                if (e.ch[ch].timing_method == 2) n_walk++;
            }
        }
        std::cout << "  CH" << ch+1 << ": clip=" << n_clip << "/" << ntot
                  << " (" << (int)(100.0*n_clip/ntot) << "%)  CFD=" << n_cfd
                  << " Walk=" << n_walk << std::endl;
    }

    if (ngood >= 2) {
        auto stats = [](const std::vector<double> &v) {
            double s=0, s2=0;
            for (double x : v) { s += x; s2 += x*x; }
            double m = s / v.size();
            double sd = sqrt(fabs(s2/v.size() - m*m));
            return std::make_pair(m, sd);
        };
        std::pair<double,double> r12 = stats(d12);
        std::pair<double,double> r13 = stats(d13);
        std::pair<double,double> r23 = stats(d23);
        std::cout << "\n  dt(2-1) = " << Form("%+.2f +/- %.2f ns", r12.first, r12.second) << std::endl;
        std::cout << "  dt(3-1) = " << Form("%+.2f +/- %.2f ns", r13.first, r13.second) << std::endl;
        std::cout << "  dt(3-2) = " << Form("%+.2f +/- %.2f ns", r23.first, r23.second) << std::endl;

        std::cout << "\n  NOTA: dt misti (un canale CFD, l'altro Walk) contengono un" << std::endl;
        std::cout << "  offset sistematico ~0.5-1 ns. Per correggerlo, usare gli eventi" << std::endl;
        std::cout << "  non clippati per calibrare <t_cfd - t_walk> per canale." << std::endl;
    }
    std::cout << "====================================\n" << std::endl;
}


// ================== OVERLAY: tutti i canali =================================
// Overlay(): sovrappone TUTTI i canali dell'evento corrente su un unico grafico.
//
// Parametro zoom_pmt:
//   false (default) → il range Y include tutto (anche il NIM a -500 mV)
//   true → il range Y esclude i canali clippati → zoom sui soli PMT
//
// Apre un canvas SEPARATO che non sovrascrive il browser principale.

void Overlay(bool zoom_pmt = false) {
    if (gEvents.empty()) return;
    EventData &evt = gEvents[gCurrentEvent];
    int nch = evt.nchannels;

    // Crea un canvas separato
    TCanvas *covl = new TCanvas("covl",
        Form("Overlay  |  Ev. %d (serial %d)  |  %s",
             gCurrentEvent+1, evt.serial, gFilename.c_str()),
        1000, 600);
    covl->SetGrid(1, 1);     // Griglia X e Y
    covl->SetTickx(1);       // Tick anche sull'asse superiore
    covl->SetTicky(1);       // Tick anche sull'asse destro
    covl->SetMargin(0.08, 0.03, 0.10, 0.08);  // margini: sx, dx, basso, alto

    // --- Calcolo range globali (scandendo TUTTI i campioni di TUTTI i canali) ---
    double tmin_g = 1e9, tmax_g = -1e9;       // Range tempo globale
    double vmin_g = 1e9, vmax_g = -1e9;       // Range tensione globale (tutti)
    double vmin_pmt = 1e9, vmax_pmt = -1e9;   // Range tensione solo PMT (non clippati)

    for (int i = 0; i < nch; i++) {
        ChannelData &cd = evt.ch[i];
        for (int s = 0; s < cd.nsamples; s++) {
            if (cd.time[s]    < tmin_g) tmin_g = cd.time[s];
            if (cd.time[s]    > tmax_g) tmax_g = cd.time[s];
            if (cd.voltage[s] < vmin_g) vmin_g = cd.voltage[s];
            if (cd.voltage[s] > vmax_g) vmax_g = cd.voltage[s];
            // Range solo PMT: esclude i canali clippati (es. NIM a -500 mV)
            if (!evt.ch[i].is_clipped) {
                if (cd.voltage[s] < vmin_pmt) vmin_pmt = cd.voltage[s];
                if (cd.voltage[s] > vmax_pmt) vmax_pmt = cd.voltage[s];
            }
        }
    }

    // Scegli il range Y in base a zoom_pmt
    double ylo, yhi;
    if (zoom_pmt && vmin_pmt < 1e8) {     // vmin_pmt < 1e8 = c'è almeno un PMT
        double vr = vmax_pmt - vmin_pmt;
        if (vr < 5) vr = 5;
        ylo = vmin_pmt - 0.15 * vr;
        yhi = vmax_pmt + 0.20 * vr;
    } else {
        double vr = vmax_g - vmin_g;
        if (vr < 5) vr = 5;
        ylo = vmin_g - 0.10 * vr;
        yhi = vmax_g + 0.20 * vr;
    }

    // --- Disegno delle forme d'onda sovrapposte ---
    bool first = true;     // Il primo grafico viene disegnato con "AL" (crea gli assi),
                           // i successivi con "L same" (aggiungono al grafico esistente)

    TLegend *leg = new TLegend(0.70, 0.75, 0.95, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);      // Sfondo trasparente
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    for (int i = 0; i < nch; i++) {
        ChannelData &cd = evt.ch[i];
        int ch_id = evt.channel_ids[i];

        TGraph *gr = new TGraph();
        for (int s = 0; s < cd.nsamples; s++)
            gr->SetPoint(s, cd.time[s], cd.voltage[s]);

        std::string label = gChannelLabels.count(ch_id)
                            ? gChannelLabels[ch_id] : Form("CHN%d", ch_id);

        gr->SetLineColor(CH_COLORS[i % 4]);
        gr->SetLineWidth(2);

        if (first) {
            gr->SetTitle(Form("Overlay  |  Ev. %d  |  %s;t [ns];V [mV]",
                              evt.serial, evt.timestamp.c_str()));
            gr->Draw("AL");
            gr->GetYaxis()->SetRangeUser(ylo, yhi);
            first = false;
        } else {
            gr->Draw("L same");
        }

        leg->AddEntry(gr, label.c_str(), "l");

        // Marker di timing (colore del canale, stile dipende dal metodo)
        if (cd.timing_method == 1 && cd.cfd_ok) {
            // CFD: linea tratteggiata corta
            TLine *lcfd = new TLine(cd.t_cfd, ylo, cd.t_cfd, cd.baseline);
            lcfd->SetLineColor(CH_COLORS[i % 4]);
            lcfd->SetLineStyle(2);
            lcfd->SetLineWidth(2);
            lcfd->Draw("same");
        } else if (cd.timing_method == 2) {
            // Walk corrected: linea tratteggiata lunga
            TLine *lwc = new TLine(cd.t_walk_corrected, ylo,
                                   cd.t_walk_corrected, cd.baseline);
            lwc->SetLineColor(CH_COLORS[i % 4]);
            lwc->SetLineStyle(7);
            lwc->SetLineWidth(2);
            lwc->Draw("same");
        }
    }

    leg->Draw("same");

    // Annotazioni dt CFD in basso a sinistra (coordinate NDC)
    TLatex *ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.032);
    ltx->SetTextColor(14);   // Grigio scuro

    double ytxt = 0.28;
    for (int i = 0; i < nch; i++) {
        for (int j = i+1; j < nch; j++) {
            // Salta CH4 (trigger NIM)
            if (evt.channel_ids[i] == 4 || evt.channel_ids[j] == 4) continue;
            if (evt.ch[i].t_timing > -900.0 && evt.ch[j].t_timing > -900.0) {
                double dt = evt.ch[j].t_timing - evt.ch[i].t_timing;
                const char* m_i = evt.ch[i].timing_method == 1 ? "C" : "W";
                const char* m_j = evt.ch[j].timing_method == 1 ? "C" : "W";
                ltx->DrawLatex(0.12, ytxt,
                    Form("#Deltat(CHN%d-CHN%d) = %+.2f ns [%s,%s]",
                         evt.channel_ids[j], evt.channel_ids[i], dt, m_j, m_i));
                ytxt -= 0.04;
            }
        }
    }

    covl->Modified();
    covl->Update();
}


// ============ OVERLAY NORMALIZZATO: tutti i canali alla stessa altezza ======
// OverlayNorm(): come Overlay(), ma ogni forma d'onda viene riscalata
// in modo che baseline → 0 e picco → -1. Questo permette di confrontare
// SOLO i tempi d'arrivo, indipendentemente dalle ampiezze.
//
// Formula di normalizzazione:
//   V_norm(t) = (V(t) - baseline) / amplitude
//
// dove amplitude = baseline - V_min (positiva).

void OverlayNorm() {
    if (gEvents.empty()) return;
    EventData &evt = gEvents[gCurrentEvent];
    int nch = evt.nchannels;

    TCanvas *cnrm = new TCanvas("cnrm",
        Form("Overlay normalizzato  |  Ev. %d (serial %d)",
             gCurrentEvent+1, evt.serial),
        1000, 600);
    cnrm->SetGrid(1, 1);
    cnrm->SetTickx(1);
    cnrm->SetTicky(1);
    cnrm->SetMargin(0.08, 0.03, 0.10, 0.08);

    bool first = true;
    TLegend *leg = new TLegend(0.70, 0.72, 0.95, 0.92);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    for (int i = 0; i < nch; i++) {
        ChannelData &cd = evt.ch[i];
        int ch_id = evt.channel_ids[i];

        TGraph *gr = new TGraph();

        if (cd.has_pulse && cd.amplitude > 0.1) {
            // Normalizzazione: baseline diventa 0, picco diventa -1
            for (int s = 0; s < cd.nsamples; s++) {
                double v_norm = (cd.voltage[s] - cd.baseline) / cd.amplitude;
                gr->SetPoint(s, cd.time[s], v_norm);
            }
        } else {
            // Se non c'è impulso, mostra il segnale grezzo diviso per 100
            // (approssimazione per non avere un grafico vuoto)
            for (int s = 0; s < cd.nsamples; s++)
                gr->SetPoint(s, cd.time[s], cd.voltage[s] / 100.0);
        }

        std::string label = gChannelLabels.count(ch_id)
                            ? gChannelLabels[ch_id] : Form("CHN%d", ch_id);

        gr->SetLineColor(CH_COLORS[i % 4]);
        gr->SetLineWidth(2);

        if (first) {
            gr->SetTitle(Form("Overlay normalizzato  |  Ev. %d  |  %s;"
                              "t [ns];(V - baseline) / amplitude",
                              evt.serial, evt.timestamp.c_str()));
            gr->Draw("AL");
            gr->GetYaxis()->SetRangeUser(-1.3, 0.35);  // Range fisso per la vista normalizzata
            first = false;
        } else {
            gr->Draw("L same");
        }

        leg->AddEntry(gr, label.c_str(), "l");

        // Marker timing al livello normalizzato
        if (cd.timing_method == 1 && cd.cfd_ok) {
            TLine *lcfd = new TLine(cd.t_cfd, -1.2, cd.t_cfd, 0.0);
            lcfd->SetLineColor(CH_COLORS[i % 4]);
            lcfd->SetLineStyle(2);
            lcfd->SetLineWidth(2);
            lcfd->Draw("same");
        } else if (cd.timing_method == 2 && cd.t_walk_corrected > -900.0) {
            TLine *lwc = new TLine(cd.t_walk_corrected, -1.2, cd.t_walk_corrected, 0.0);
            lwc->SetLineColor(CH_COLORS[i % 4]);
            lwc->SetLineStyle(7);
            lwc->SetLineWidth(2);
            lwc->Draw("same");
        }
    }

    // Linea orizzontale al livello CFD 50% (= -0.5 nella scala normalizzata)
    TLine *lcfd_level = new TLine(0, -CFD_FRACTION, 200, -CFD_FRACTION);
    lcfd_level->SetLineColor(CFD_COLOR);
    lcfd_level->SetLineStyle(7);       // Tratteggiata lunga
    lcfd_level->SetLineWidth(1);
    lcfd_level->Draw("same");

    // Etichetta "CFD 50%" vicino alla linea
    TLatex *lbl_cfd = new TLatex();
    lbl_cfd->SetNDC();
    lbl_cfd->SetTextFont(42);
    lbl_cfd->SetTextSize(0.030);
    lbl_cfd->SetTextColor(CFD_COLOR);
    lbl_cfd->DrawLatex(0.85, 0.43, "CFD 50%");

    leg->Draw("same");

    // Annotazioni dt (stessa logica di Overlay)
    TLatex *ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.032);
    ltx->SetTextColor(14);

    double ytxt = 0.28;
    for (int i = 0; i < nch; i++) {
        for (int j = i+1; j < nch; j++) {
            // Salta CH4 (trigger NIM)
            if (evt.channel_ids[i] == 4 || evt.channel_ids[j] == 4) continue;
            if (evt.ch[i].t_timing > -900.0 && evt.ch[j].t_timing > -900.0) {
                double dt = evt.ch[j].t_timing - evt.ch[i].t_timing;
                const char* m_i = evt.ch[i].timing_method == 1 ? "C" : "W";
                const char* m_j = evt.ch[j].timing_method == 1 ? "C" : "W";
                ltx->DrawLatex(0.12, ytxt,
                    Form("#Deltat(CHN%d-CHN%d) = %+.2f ns [%s,%s]",
                         evt.channel_ids[j], evt.channel_ids[i], dt, m_j, m_i));
                ytxt -= 0.04;
            }
        }
    }

    cnrm->Modified();
    cnrm->Update();
}


// ============ OVERLAY SELETTIVO con auto-scaling ottimale ==================
// OverlaySelect(): la funzione più flessibile. Permette di scegliere
// esattamente quali canali sovrapporre, se normalizzare, e se zoomare.
//
// Parametri:
//   channels  → stringa "1,2,3" con i numeri dei canali (1-based)
//   normalize → se true, normalizza (BL→0, picco→-1)
//   zoom      → se true, zoom asse X sulla regione dell'impulso

void OverlaySelect(const char* channels = "1,2,3",
                   bool normalize = false,
                   bool zoom = false) {

    if (gEvents.empty()) {
        std::cout << "[ERRORE] Nessun evento caricato." << std::endl;
        return;
    }
    EventData &evt = gEvents[gCurrentEvent];

    // --- 1. PARSING DELLA SELEZIONE CANALI ---
    // Converte "1,2,3" → vettore di indici in evt.ch[]

    std::vector<int> sel_idx;   // Indici 0-based nel vettore evt.ch[]
    std::vector<int> sel_id;    // Numeri CHN originali (1-based)

    // Parsing manuale della stringa separata da virgole
    // (evita istringstream e il problema del "most vexing parse")
    std::string ch_str(channels);
    std::string token;
    std::string::size_type pos = 0, found;

    while (pos < ch_str.size()) {
        found = ch_str.find(',', pos);                       // Cerca la prossima virgola
        if (found == std::string::npos) found = ch_str.size(); // Ultimo token
        token = ch_str.substr(pos, found - pos);              // Estrai il token
        pos = found + 1;                                      // Avanza dopo la virgola

        int ch_num = atoi(token.c_str());                    // Converti "3" → 3
        if (ch_num < 1) continue;                             // Ignora valori non validi

        // Cerca questo numero di canale nell'evento corrente
        bool found_ch = false;
        for (int i = 0; i < evt.nchannels; i++) {
            if (evt.channel_ids[i] == ch_num) {
                sel_idx.push_back(i);         // Salva l'indice 0-based
                sel_id.push_back(ch_num);     // Salva il numero 1-based
                found_ch = true;
                break;
            }
        }
        if (!found_ch) {
            std::cout << "[WARNING] CHN" << ch_num
                      << " non presente nell'evento " << evt.serial << std::endl;
        }
    }

    if (sel_idx.empty()) {
        std::cout << "[ERRORE] Nessun canale valido selezionato." << std::endl;
        return;
    }

    // --- 2. CALCOLO DEI RANGE OTTIMALI ---

    double vmin_sel = 1e9, vmax_sel = -1e9;           // Range Y dei canali selezionati
    double tmin_pulse = 1e9, tmax_pulse = -1e9;       // Range X della regione impulso
    double t_start = 1e9, t_end = -1e9;               // Range X globale dei dati

    for (int k = 0; k < (int)sel_idx.size(); k++) {
        int i = sel_idx[k];
        ChannelData &cd = evt.ch[i];

        for (int s = 0; s < cd.nsamples; s++) {
            double v_val = cd.voltage[s];
            double t_val = cd.time[s];
            if (v_val < vmin_sel) vmin_sel = v_val;
            if (v_val > vmax_sel) vmax_sel = v_val;
            if (t_val < t_start) t_start = t_val;
            if (t_val > t_end)   t_end = t_val;
        }

        // Per il zoom X: trova dove c'è segnale significativo
        if (cd.has_pulse) {
            double thr = cd.baseline - 5.0 * cd.baseline_rms;
            if (thr > cd.baseline - NOISE_THRESH) thr = cd.baseline - NOISE_THRESH;
            for (int s = 0; s < cd.nsamples; s++) {
                if (cd.voltage[s] < thr) {
                    if (cd.time[s] < tmin_pulse) tmin_pulse = cd.time[s];
                    if (cd.time[s] > tmax_pulse) tmax_pulse = cd.time[s];
                }
            }
        }
    }

    // Range Y
    double ylo, yhi;
    if (normalize) {
        ylo = -1.30;    // Fisso per la vista normalizzata
        yhi = 0.35;
    } else {
        double vr = vmax_sel - vmin_sel;
        if (vr < 5.0) vr = 5.0;
        ylo = vmin_sel - 0.08 * vr;    // 8% sotto
        yhi = vmax_sel + 0.22 * vr;    // 22% sopra (spazio per legenda)
    }

    // Range X
    double xlo, xhi;
    if (zoom && tmin_pulse < 1e8) {   // tmin_pulse < 1e8 = trovato almeno un impulso
        double pulse_width = tmax_pulse - tmin_pulse;
        if (pulse_width < 10.0) pulse_width = 10.0;
        double margin = pulse_width * 1.5;      // 150% di margine per lato
        if (margin < 20.0) margin = 20.0;        // Minimo 20 ns per lato
        xlo = tmin_pulse - margin;
        xhi = tmax_pulse + margin;
        if (xlo < t_start) xlo = t_start;        // Non andare sotto il primo campione
        if (xhi > t_end)   xhi = t_end;          // Non andare oltre l'ultimo
    } else {
        xlo = t_start;
        xhi = t_end;
    }

    // --- 3. CREAZIONE CANVAS E DISEGNO ---
    static int canvas_counter = 0;    // Contatore per nomi canvas unici
    canvas_counter++;

    TCanvas *csel = new TCanvas(
        Form("csel_%d", canvas_counter),
        Form("Overlay [CHN %s] %s%s |  Ev. %d (serial %d)",
             channels,
             normalize ? "NORM " : "",
             zoom ? "ZOOM " : "",
             gCurrentEvent+1, evt.serial),
        1100, 650);
    csel->SetGrid(1, 1);
    csel->SetTickx(1);
    csel->SetTicky(1);
    csel->SetMargin(0.08, 0.03, 0.10, 0.06);

    TLegend *leg = new TLegend(0.72, 0.78, 0.96, 0.93);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    bool first = true;

    for (int k = 0; k < (int)sel_idx.size(); k++) {
        int i = sel_idx[k];
        int ch_id = sel_id[k];
        ChannelData &cd = evt.ch[i];

        TGraph *gr = new TGraph();

        if (normalize) {
            if (cd.has_pulse && cd.amplitude > 0.1) {
                for (int s = 0; s < cd.nsamples; s++) {
                    double v_norm = (cd.voltage[s] - cd.baseline) / cd.amplitude;
                    gr->SetPoint(s, cd.time[s], v_norm);
                }
            } else {
                for (int s = 0; s < cd.nsamples; s++)
                    gr->SetPoint(s, cd.time[s], (cd.voltage[s] - cd.baseline) / 100.0);
            }
        } else {
            for (int s = 0; s < cd.nsamples; s++)
                gr->SetPoint(s, cd.time[s], cd.voltage[s]);
        }

        std::string label = gChannelLabels.count(ch_id)
                            ? gChannelLabels[ch_id] : Form("CHN%d", ch_id);

        gr->SetLineColor(CH_COLORS[i % 4]);
        gr->SetLineWidth(2);

        if (first) {
            std::string title_str = "Overlay";
            if (normalize) title_str += " normalizzato";
            if (zoom) title_str += " (zoom)";
            title_str += Form("  |  Ev. %d  |  %s", evt.serial, evt.timestamp.c_str());

            const char* yaxis_title = normalize ? "(V - baseline) / amplitude" : "V [mV]";

            gr->SetTitle(Form("%s;t [ns];%s", title_str.c_str(), yaxis_title));
            gr->Draw("AL");
            gr->GetYaxis()->SetRangeUser(ylo, yhi);
            gr->GetXaxis()->SetRangeUser(xlo, xhi);
            first = false;
        } else {
            gr->Draw("L same");
        }

        // In legenda: etichetta + ampiezza + metodo timing
        const char* meth_tag = cd.timing_method == 1 ? "CFD" :
                               cd.timing_method == 2 ? "WC"  : "--";
        leg->AddEntry(gr, Form("%s (A=%.0f mV, %s)", label.c_str(),
                               cd.amplitude, meth_tag), "l");

        // Marker timing (stile dipende dal metodo)
        if (cd.timing_method == 1 && cd.cfd_ok) {
            // CFD: tratteggio corto
            double bot = normalize ? -1.2 : ylo;
            double top = normalize ? 0.0  : cd.baseline;
            TLine *lcfd = new TLine(cd.t_cfd, bot, cd.t_cfd, top);
            lcfd->SetLineColor(CH_COLORS[i % 4]);
            lcfd->SetLineStyle(2);
            lcfd->SetLineWidth(2);
            lcfd->Draw("same");
        } else if (cd.timing_method == 2 && cd.t_walk_corrected > -900.0) {
            // Walk corrected: tratteggio lungo
            double bot = normalize ? -1.2 : ylo;
            double top = normalize ? 0.0  : cd.baseline;
            TLine *lwc = new TLine(cd.t_walk_corrected, bot,
                                   cd.t_walk_corrected, top);
            lwc->SetLineColor(CH_COLORS[i % 4]);
            lwc->SetLineStyle(7);
            lwc->SetLineWidth(2);
            lwc->Draw("same");
        }
    }

    // Linea soglia CFD (solo normalizzato)
    if (normalize) {
        TLine *lcfd_lev = new TLine(xlo, -CFD_FRACTION, xhi, -CFD_FRACTION);
        lcfd_lev->SetLineColor(CFD_COLOR);
        lcfd_lev->SetLineStyle(7);
        lcfd_lev->SetLineWidth(1);
        lcfd_lev->Draw("same");

        TLatex *lcfd_txt = new TLatex();
        lcfd_txt->SetNDC();
        lcfd_txt->SetTextFont(42);
        lcfd_txt->SetTextSize(0.028);
        lcfd_txt->SetTextColor(CFD_COLOR);
        lcfd_txt->DrawLatex(0.88, 0.42, "CFD 50%");
    }

    leg->Draw("same");

    // --- 4. ANNOTAZIONI dt (timing unificato) ---
    TLatex *ltx = new TLatex();
    ltx->SetNDC();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.033);

    double ytxt = 0.30;
    for (int k1 = 0; k1 < (int)sel_idx.size(); k1++) {
        for (int k2 = k1+1; k2 < (int)sel_idx.size(); k2++) {
            // Salta CH4 (trigger NIM)
            if (sel_id[k1] == 4 || sel_id[k2] == 4) continue;
            int i1 = sel_idx[k1], i2 = sel_idx[k2];
            if (evt.ch[i1].t_timing > -900.0 && evt.ch[i2].t_timing > -900.0) {
                double dt = evt.ch[i2].t_timing - evt.ch[i1].t_timing;
                const char* m1 = evt.ch[i1].timing_method == 1 ? "C" : "W";
                const char* m2 = evt.ch[i2].timing_method == 1 ? "C" : "W";
                ltx->SetTextColor(CH_COLORS[i2 % 4]);
                ltx->DrawLatex(0.10, ytxt,
                    Form("#Deltat(CHN%d #minus CHN%d) = %+.2f ns [%s,%s]",
                         sel_id[k2], sel_id[k1], dt, m2, m1));
                ytxt -= 0.042;
            }
        }
    }

    csel->Modified();
    csel->Update();
}


// ==================== FUNZIONI CABLE OFFSET ==================================

/// SetCableOffsets(): imposta manualmente gli offset dei cavi [ns].
/// I due argomenti sono Δ(delay_T06 - delay_T08) e Δ(delay_T04 - delay_T08).
/// Dopo la modifica, ricalcola la ricostruzione angolare per TUTTI gli eventi.
void SetCableOffsets(double offset_06_08, double offset_04_08) {
    gCableOffset_06_08 = offset_06_08;
    gCableOffset_04_08 = offset_04_08;
    std::cout << "[INFO] Cable offsets aggiornati:" << std::endl;
    std::cout << "  Δ(T06-T08) = " << offset_06_08 << " ns" << std::endl;
    std::cout << "  Δ(T04-T08) = " << offset_04_08 << " ns" << std::endl;

    // Ricalcola la ricostruzione per tutti gli eventi
    int n_ok = 0;
    for (auto &e : gEvents) {
        ReconstructEvent(e);
        if (e.reco_ok) n_ok++;
    }
    std::cout << "  Ricostruzione ricalcolata: " << n_ok << "/" << gEvents.size()
              << " eventi OK" << std::endl;

    // Aggiorna il display
    DrawEvent(gCurrentEvent);
}

/// CalibrateCableOffsets(): stima gli offset totali dalla media dei Δt misurati.
/// L'idea è che la distribuzione angolare degli sciami è centrata sullo zenith,
/// quindi ⟨Δt_phys⟩ ≈ Δt_verticale ≈ Δz/c (piccolo, ~1 ns), e:
///    cable_offset ≈ ⟨Δt_misurato⟩ - Δz/c
/// Questa stima include TUTTI i ritardi strumentali (cavi + PMT + elettronica).
void CalibrateCableOffsets() {
    if (gEvents.empty()) {
        std::cout << "[ERRORE] Nessun evento caricato." << std::endl;
        return;
    }

    double sum_dt1 = 0, sum_dt2 = 0;
    int n = 0;
    for (auto &e : gEvents) {
        if (e.nchannels < 3) continue;
        if (e.ch[CH_IDX_T08].t_timing < -900.0 ||
            e.ch[CH_IDX_T06].t_timing < -900.0 ||
            e.ch[CH_IDX_T04].t_timing < -900.0) continue;
        sum_dt1 += e.ch[CH_IDX_T06].t_timing - e.ch[CH_IDX_T08].t_timing;
        sum_dt2 += e.ch[CH_IDX_T04].t_timing - e.ch[CH_IDX_T08].t_timing;
        n++;
    }

    if (n < 2) {
        std::cout << "[ERRORE] Troppi pochi eventi con 3 timing validi (" << n << ")." << std::endl;
        return;
    }

    double mean_dt1 = sum_dt1 / n;
    double mean_dt2 = sum_dt2 / n;

    // Correzione per lo sciame verticale medio:
    // Δt_vert = -Δz / c  (segno: dal modello ŝ·d = c×Δt con ŝ=(0,0,1) per verticale)
    double dt_vert_1 = BASELINE_06_08[2] / C_LIGHT;  // -0.19/0.3 ≈ -0.63 ns
    double dt_vert_2 = BASELINE_04_08[2] / C_LIGHT;  //  0.735/0.3 ≈ +2.45 ns

    double new_offset_1 = mean_dt1 - dt_vert_1;
    double new_offset_2 = mean_dt2 - dt_vert_2;

    std::cout << "\n[INFO] Calibrazione cable offsets da " << n << " eventi:" << std::endl;
    std::cout << "  ⟨Δt_meas(06-08)⟩ = " << Form("%.2f ns", mean_dt1) << std::endl;
    std::cout << "  ⟨Δt_meas(04-08)⟩ = " << Form("%.2f ns", mean_dt2) << std::endl;
    std::cout << "  Δt_vert(06-08) = " << Form("%.2f ns", dt_vert_1)
              << "  (da Δz = " << BASELINE_06_08[2] << " m)" << std::endl;
    std::cout << "  Δt_vert(04-08) = " << Form("%.2f ns", dt_vert_2)
              << "  (da Δz = " << BASELINE_04_08[2] << " m)" << std::endl;
    std::cout << "  Vecchi offsets: " << gCableOffset_06_08 << ", " << gCableOffset_04_08 << std::endl;

    // Applica i nuovi offsets
    SetCableOffsets(new_offset_1, new_offset_2);
}


// ==================== SKYMAP E DISTRIBUZIONI ANGOLARI =========================
// Skymap(): produce un canvas con 4 plot della distribuzione angolare
// degli sciami ricostruiti:
//   (1) Mappa del cielo in coordinate (l, m) — proiezione zenitale
//   (2) dN/dθ con fit ∝ cos²θ sinθ
//   (3) dN/d(cosθ) con fit ∝ cosⁿθ
//   (4) dN/dφ — distribuzione azimutale
//
// La proiezione zenitale (l, m) è lo standard nella fisica dei raggi cosmici:
//   l = sinθ cosφ,  m = sinθ sinφ
// Il centro (0,0) è lo zenith, il cerchio di raggio 1 è l'orizzonte.
// Cerchi concentrici corrispondono ad angoli zenitali costanti.

void Skymap() {
    if (gEvents.empty()) {
        std::cout << "[ERRORE] Nessun evento caricato." << std::endl;
        return;
    }

    // --- Raccolta eventi ricostruiti ---
    std::vector<double> v_theta, v_phi, v_l, v_m, v_costheta;
    for (auto &e : gEvents) {
        if (!e.reco_ok) continue;
        v_theta.push_back(e.reco_theta);
        v_phi.push_back(e.reco_phi);
        v_l.push_back(e.reco_l);
        v_m.push_back(e.reco_m);
        v_costheta.push_back(cos(e.reco_theta));
    }

    int nreco = v_theta.size();
    std::cout << "\n[INFO] Skymap: " << nreco << "/" << gEvents.size()
              << " eventi ricostruiti" << std::endl;

    if (nreco < 1) {
        std::cout << "[ERRORE] Nessun evento ricostruito. Controllare i cable offsets" << std::endl;
        std::cout << "  (usare CalibrateCableOffsets() o SetCableOffsets(a,b))." << std::endl;
        return;
    }

    // --- Canvas 2×2 ---
    TCanvas *csky = new TCanvas("csky",
        Form("Distribuzioni angolari  |  %d eventi  |  %s",
             nreco, gFilename.c_str()),
        1400, 1000);

    // ===================================================================
    //  PAD 1 (alto-sinistra): SKYMAP in coordinate (l, m)
    // ===================================================================
    TPad *p1 = new TPad("p_skymap", "", 0.0, 0.5, 0.5, 1.0);
    p1->SetMargin(0.12, 0.04, 0.12, 0.10);
    p1->SetGrid(1, 1);
    p1->SetTickx(1);
    p1->SetTicky(1);
    p1->Draw();
    p1->cd();

    // Scatter plot dei punti ricostruiti
    TGraph *gr_sky = new TGraph(nreco);
    for (int i = 0; i < nreco; i++)
        gr_sky->SetPoint(i, v_l[i], v_m[i]);

    gr_sky->SetTitle("Mappa del cielo (proiezione zenitale);l = sin#theta cos#phi;m = sin#theta sin#phi");
    gr_sky->SetMarkerStyle(20);    // Cerchio pieno
    gr_sky->SetMarkerSize(1.2);
    gr_sky->SetMarkerColor(kBlue+1);
    gr_sky->Draw("AP");

    // Range simmetrico centrato sull'origine
    double lm_max = 0.3;
    for (int i = 0; i < nreco; i++) {
        if (fabs(v_l[i]) > lm_max) lm_max = fabs(v_l[i]) * 1.3;
        if (fabs(v_m[i]) > lm_max) lm_max = fabs(v_m[i]) * 1.3;
    }
    if (lm_max > 1.0) lm_max = 1.0;
    gr_sky->GetXaxis()->SetRangeUser(-lm_max, lm_max);
    gr_sky->GetYaxis()->SetRangeUser(-lm_max, lm_max);

    // Cerchi di θ costante (guide visive)
    for (double theta_circ : {10.0, 20.0, 30.0, 45.0, 60.0}) {
        double r = sin(theta_circ * M_PI / 180.0);
        if (r > lm_max) continue;
        TEllipse *ell = new TEllipse(0, 0, r, r);
        ell->SetFillStyle(0);        // Trasparente
        ell->SetLineColor(kGray);
        ell->SetLineStyle(2);        // Tratteggiata
        ell->Draw("same");
        // Etichetta
        TLatex *lt = new TLatex(r + 0.01, 0.01, Form("%.0f#circ", theta_circ));
        lt->SetTextSize(0.028);
        lt->SetTextColor(kGray+1);
        lt->Draw("same");
    }

    // Punto dello zenith
    TGraph *zenith = new TGraph(1);
    zenith->SetPoint(0, 0, 0);
    zenith->SetMarkerStyle(5);     // Croce
    zenith->SetMarkerSize(1.5);
    zenith->SetMarkerColor(kGray+2);
    zenith->Draw("P same");

    p1->Modified();
    p1->Update();

    // ===================================================================
    //  PAD 2 (alto-destra): dN/dθ
    // ===================================================================
    csky->cd();
    TPad *p2 = new TPad("p_theta", "", 0.5, 0.5, 1.0, 1.0);
    p2->SetMargin(0.12, 0.04, 0.12, 0.10);
    p2->SetGrid(1, 1);
    p2->Draw();
    p2->cd();

    // Istogramma: 18 bin da 0° a 90° (5° per bin)
    TH1D *h_theta = new TH1D("h_theta", "dN/d#theta;#theta [#circ];Eventi / 5#circ",
                              18, 0, 90);
    h_theta->SetLineColor(kBlue+1);
    h_theta->SetLineWidth(2);
    h_theta->SetFillColor(kBlue-9);
    h_theta->SetFillStyle(1001);
    for (double th : v_theta)
        h_theta->Fill(th * 180.0 / M_PI);
    h_theta->Draw("HIST");

    // Fit con f(θ) = A × cos²θ × sinθ  (forma attesa per raggi cosmici)
    // In gradi: f(x) = A × cos²(x·π/180) × sin(x·π/180)
    if (nreco >= 5) {
        TF1 *f_theta = new TF1("f_theta",
            "[0] * pow(cos(x*TMath::Pi()/180.0), [1]) * sin(x*TMath::Pi()/180.0)",
            0, 90);
        f_theta->SetParameters(nreco * 5.0 / 18.0, 2.0);  // [0]=normalizzazione, [1]=esponente n
        f_theta->SetParNames("Norm", "n");
        f_theta->SetLineColor(kRed);
        f_theta->SetLineWidth(2);
        h_theta->Fit(f_theta, "QR");   // Q=quiet, R=range
    }

    p2->Modified();
    p2->Update();

    // ===================================================================
    //  PAD 3 (basso-sinistra): dN/d(cosθ)
    // ===================================================================
    //  In questa rappresentazione il fattore geometrico sinθ dell'angolo
    //  solido è assorbito nel cambio di variabile (dΩ = d(cosθ) dφ),
    //  quindi il flusso di raggi cosmici appare come una semplice legge
    //  di potenza: dN/d(cosθ) ∝ cosⁿθ, con n ≈ 2 atteso.
    csky->cd();
    TPad *p3 = new TPad("p_costheta", "", 0.0, 0.0, 0.5, 0.5);
    p3->SetMargin(0.12, 0.04, 0.12, 0.10);
    p3->SetGrid(1, 1);
    p3->Draw();
    p3->cd();

    // 20 bin uniformi in cosθ da 0 a 1
    TH1D *h_costheta = new TH1D("h_costheta",
        "dN/d(cos#theta);cos#theta;Eventi / 0.05",
        20, 0, 1);
    h_costheta->SetLineColor(kGreen+2);
    h_costheta->SetLineWidth(2);
    h_costheta->SetFillColor(kGreen-9);
    h_costheta->SetFillStyle(1001);
    for (double ct : v_costheta)
        h_costheta->Fill(ct);
    h_costheta->Draw("HIST");

    // Fit con f(cosθ) = A × cosⁿθ
    if (nreco >= 5) {
        TF1 *f_costh = new TF1("f_costh", "[0] * pow(x, [1])", 0, 1);
        f_costh->SetParameters(nreco * 0.05, 2.0);
        f_costh->SetParNames("Norm", "n");
        f_costh->SetLineColor(kRed);
        f_costh->SetLineWidth(2);
        h_costheta->Fit(f_costh, "QR");
    }

    p3->Modified();
    p3->Update();

    // ===================================================================
    //  PAD 4 (basso-destra): dN/dφ
    // ===================================================================
    //  Per raggi cosmici, la distribuzione azimutale è quasi uniforme
    //  (isotropia in φ), salvo piccoli effetti dell'asimmetria Est-Ovest
    //  dovuta al campo geomagnetico.
    csky->cd();
    TPad *p4 = new TPad("p_phi", "", 0.5, 0.0, 1.0, 0.5);
    p4->SetMargin(0.12, 0.04, 0.12, 0.10);
    p4->SetGrid(1, 1);
    p4->Draw();
    p4->cd();

    // 12 bin da 0° a 360° (30° per bin, corrispondente a ~ore di un orologio)
    TH1D *h_phi = new TH1D("h_phi", "dN/d#phi;#phi [#circ];Eventi / 30#circ",
                            12, 0, 360);
    h_phi->SetLineColor(kMagenta+1);
    h_phi->SetLineWidth(2);
    h_phi->SetFillColor(kMagenta-9);
    h_phi->SetFillStyle(1001);
    for (double ph : v_phi)
        h_phi->Fill(ph * 180.0 / M_PI);
    h_phi->Draw("HIST");

    // Linea orizzontale al valore medio atteso (isotropo: N/12 per bin)
    if (nreco > 0) {
        double expected = (double)nreco / 12.0;
        TLine *l_flat = new TLine(0, expected, 360, expected);
        l_flat->SetLineColor(kRed);
        l_flat->SetLineStyle(2);
        l_flat->SetLineWidth(2);
        l_flat->Draw("same");
    }

    p4->Modified();
    p4->Update();

    // --- Sommario a terminale ---
    std::cout << "\n=== DISTRIBUZIONI ANGOLARI ===" << std::endl;
    std::cout << "  Eventi ricostruiti: " << nreco << "/" << gEvents.size() << std::endl;
    std::cout << "  Cable offsets usati: Δ(06-08) = " << gCableOffset_06_08
              << " ns, Δ(04-08) = " << gCableOffset_04_08 << " ns" << std::endl;
    if (nreco > 0) {
        double sum_th = 0, sum_th2 = 0;
        for (double th : v_theta) {
            double thd = th * 180.0 / M_PI;
            sum_th += thd;
            sum_th2 += thd * thd;
        }
        double mean_th = sum_th / nreco;
        double rms_th  = sqrt(fabs(sum_th2 / nreco - mean_th * mean_th));
        std::cout << "  ⟨θ⟩ = " << Form("%.1f ± %.1f°", mean_th, rms_th) << std::endl;
    }
    std::cout << "==============================\n" << std::endl;

    csky->Modified();
    csky->Update();
}



// DRS4Browser_v3(): funzione principale. ROOT la chiama automaticamente
// quando si lancia:   root -l 'DRS4Browser_v3.cpp("file.xml")'
//
// Il nome DEVE corrispondere al nome del file senza estensione!

void DRS4Browser(const char* xmlfile, const char* labels = "") {

    std::cout << "=============================================" << std::endl;
    std::cout << "  DRS4 Waveform Browser v3                   " << std::endl;
    std::cout << "  Timing: CFD + Walk-Corrected Leading Edge  " << std::endl;
    std::cout << "  Laboratorio Sciami Estesi, Marzo 2026      " << std::endl;
    std::cout << "=============================================" << std::endl;

    gFilename = xmlfile;
    gEvents.clear();     // Svuota eventuali dati da un lancio precedente

    // --- Parsing delle etichette personalizzate ---
    // Formato atteso: "CHN1=PMT top,CHN2=PMT mid,CHN3=PMT bot,CHN4=NIM trigger"
    gChannelLabels[1] = "CHN1";    // Default
    gChannelLabels[2] = "CHN2";
    gChannelLabels[3] = "CHN3";
    gChannelLabels[4] = "CHN4";

    if (strlen(labels) > 0) {
        // Crea una copia std::string per evitare il "most vexing parse"
        // (std::istringstream ss(std::string(labels)) verrebbe interpretato
        //  come dichiarazione di funzione dal compilatore C++)
        std::string labels_str(labels);
        std::istringstream ss(labels_str);
        std::string token;
        while (std::getline(ss, token, ',')) {
            int ch;
            char lbl[128];
            // Formato atteso: "CHN1=testo libero"
            if (sscanf(token.c_str(), "CHN%d=%127[^\n]", &ch, lbl) == 2)
                gChannelLabels[ch] = lbl;
        }
    }

    // --- Parsing del file XML ---
    if (ParseXML(xmlfile, gEvents) <= 0) return;   // Esce se il file è vuoto o non leggibile

    // --- Setup stile grafico ROOT ---
    gStyle->SetOptStat(0);              // Disabilita il box delle statistiche
    gStyle->SetTitleSize(0.06, "t");    // Titolo principale più grande
    gStyle->SetLabelSize(0.050, "xy");  // Etichette degli assi
    gStyle->SetTitleSize(0.050, "xy");  // Titoli degli assi
    gStyle->SetTitleOffset(0.9, "y");   // Distanza del titolo Y dall'asse

    // --- Determina il layout della griglia di pad ---
    // Il layout si adatta al numero di canali presenti nel primo evento
    int nch = gEvents[0].nchannels;
    int ncols, nrows;
    if (nch <= 2)      { ncols = nch; nrows = 1; }   // 1 o 2 canali: una riga
    else if (nch == 3) { ncols = 3;   nrows = 1; }   // 3 canali: 3 affiancati
    else               { ncols = 2;   nrows = 2; }   // 4 canali: griglia 2×2

    // --- Creazione del canvas principale ---
    gCanvas = new TCanvas("DRS4Browser",
        Form("DRS4 Waveform Browser  |  %s  |  %d eventi",
             xmlfile, (int)gEvents.size()),
        1400, 900);    // Larghezza 1400px, altezza 900px

    // Pannello informativo in basso (occupa il 16% dell'altezza)
    gPadInfo = new TPad("info", "", 0.0, 0.0, 1.0, 0.16);
    gPadInfo->SetFillColor(kWhite);
    gPadInfo->SetBorderSize(1);
    gPadInfo->Draw();

    // Pad per le forme d'onda (occupano la parte superiore)
    // Le coordinate sono in frazioni del canvas: (0,0) = angolo basso-sinistra, (1,1) = alto-destra
    float top = 0.98, bot = 0.16;       // Limiti verticali della zona waveform
    float left = 0.005, right = 0.995;  // Limiti orizzontali (con piccolo margine)
    float dw = (right - left) / ncols;  // Larghezza di ogni pad
    float dh = (top - bot) / nrows;     // Altezza di ogni pad

    for (int i = 0; i < nch; i++) {
        int col = i % ncols;             // Colonna (0, 1, ...)
        int row = i / ncols;             // Riga (0, 1)
        float x1 = left + col*dw;       // Bordo sinistro del pad
        float x2 = x1 + dw;             // Bordo destro
        float y2 = top - row*dh;        // Bordo superiore
        float y1 = y2 - dh;             // Bordo inferiore

        gPadWF[i] = new TPad(Form("wf%d",i), "", x1, y1, x2, y2);
        gPadWF[i]->SetMargin(0.10, 0.03, 0.12, 0.09);  // Margini interni (sx, dx, basso, alto)
        gPadWF[i]->SetGrid(1, 1);       // Griglia
        gPadWF[i]->SetTickx(1);         // Tick anche in alto
        gPadWF[i]->SetTicky(1);         // Tick anche a destra
        gPadWF[i]->Draw();
    }

    // Mostra il primo evento
    DrawEvent(0);

    // Stampa le istruzioni al terminale
    std::cout << "\n  Timing: CFD (arancione) per segnali non clippati" << std::endl;
    std::cout << "          Walk-corrected (ciano) per segnali clippati" << std::endl;
    std::cout << "          [C]=CFD, [W]=Walk nelle annotazioni dt\n" << std::endl;
    std::cout << "  Navigazione:  Next()  Prev()  GoTo(n)  GoToSerial(s)" << std::endl;
    std::cout << "  Overlay:      Overlay()  Overlay(true)  OverlayNorm()" << std::endl;
    std::cout << "                OverlaySelect(\"1,2,3\")  OverlaySelect(\"1,2,3\",true,true)" << std::endl;
    std::cout << "  Angoli:       Skymap()                 distribuzioni angolari" << std::endl;
    std::cout << "  Calibrazione: CalibrateCableOffsets()   stima offset dalla media dei dt" << std::endl;
    std::cout << "                SetCableOffsets(a,b)      imposta offset manuali [ns]" << std::endl;
    std::cout << "  Altro:        Save()  SaveAll()  Summary()\n" << std::endl;
    std::cout << "  Cable offsets attuali: d(06-08)=" << gCableOffset_06_08
              << " ns, d(04-08)=" << gCableOffset_04_08 << " ns" << std::endl;
}


// ========================== MAIN (compilato) ===============================
// Questa sezione viene inclusa SOLO quando si compila il file come
// programma standalone (g++ -o ...). Quando ROOT lo carica come macro
// con .L o con root -l 'file.cpp("...")', le guardie #ifndef escludono main().

#ifndef __CINT__       // Guardia per il vecchio interprete CINT (ROOT 5)
#ifndef __CLING__      // Guardia per il nuovo interprete Cling (ROOT 6)
int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Utilizzo: " << argv[0] << " <file.xml> [\"CHN1=...,CHN2=...\"]" << std::endl;
        return 1;
    }
    // TApplication è necessaria per la versione compilata:
    // gestisce il loop degli eventi grafici (senza di essa la finestra
    // si aprirebbe e chiuderebbe immediatamente)
    TApplication app("app", &argc, argv);
    DRS4Browser_v3(argv[1], argc > 2 ? argv[2] : "");
    app.Run();     // Entra nel loop degli eventi (blocca finché non si chiude la finestra)
    return 0;
}
#endif
#endif
