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
#include "TH1F.h"           // TH1F: istogramma 1D a precisione float (per risoluzioni)


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

};


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

static double gCableOffset_06_08 = 0. ;   // delay_T06 - delay_T08 = 205 - 17 [ns]
static double gCableOffset_04_08 =  0.;   // delay_T04 - delay_T08 = 67.5 - 17 [ns]

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
                //ReconstructEvent(current_evt);
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
/// Istogramma risoluzione temporale di ogni pinao del telescopio
/// funzione Hist_resolution(): crea un istogramma con le differenze temporali tra i canali 1-2, 1-3, 2-3
/// gEvents: vettore globale con tutti gli eventi, già popolato da ParseXML()
/// attributo t_timing di ogni canale: contiene il tempo unificato (CFD o walk-corrected) da usare per le differenze
/// ngood = numero di eventi con almeno 3 canali e timing valido (t_timing > -900 ns)
void Hist_resolution() { 
    int ntot = gEvents.size(), ngood = 0;
    std::vector<double> d12, d13, d23; // vettori per le differenze temporali

    TH1F *hist12 = new TH1F("h1", "Istogramma #Delta t 1-2", 100, -20.0, 20.0); // istogramma con 100 bin da -20 a +20 ns, delle differenze temporali tra canali 1-2, 1-3, 2-3
    TH1F *hist13 = new TH1F("h2", "Istogramma #Delta t 1-3", 100, -20.0, 20.0);
    TH1F *hist23 = new TH1F("h3", "Istogramma #Delta t 2-3", 100, -20.0, 20.0);


    //ciclo for su &e di gEvents
    //auto serve a iterare su tutti gli eventi del vettore gEvents, e &e è un riferimento all'evento corrente (per evitare copie inutili)
    for(auto &e : gEvents) {
         if(e.nchannels < 3) continue; // se l'evento ha meno di 3 canali, salta
        bool ok = true; // flag per verificare se tutti e 3 i canali hanno timing valido
        for(int i=0; i<3; i++) {
             if(e.ch[i].t_timing < -900.0) ok = false; // se il timing del canale i è < -900 ns, non è valido
    }
        if(!ok) continue; // se almeno un canale non ha timing valido, salta l'evento
        ngood++; // conta gli eventi buoni
        d12.push_back(e.ch[1].t_timing - e.ch[0].t_timing); // differenza temporale canale 2 - canale 1
        d13.push_back(e.ch[2].t_timing - e.ch[0].t_timing); // differenza temporale canale 3 - canale 1
        d23.push_back(e.ch[2].t_timing - e.ch[1].t_timing); // differenza temporale canale 3 - canale 2
        hist12->Fill(d12.back()); // riempi l'istogramma con l'ultima differenza calcolata
        hist13->Fill(d13.back());
        hist23->Fill(d23.back());

    }
    // ora faccio l'istogramma con d12, d13, d23 usando TH1D di ROOT
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600); //primo canvas per il primo istogramma
    hist12->GetXaxis()->SetTitle("#Delta t_1_2 [ns]"); // primo istogramma: differenza temporale tra canali 1-2 
    hist12->GetYaxis()->SetTitle("Counts");
    hist12->SetFillColor(kBlue-7);
    hist12->Draw();
    TCanvas *c2 = new TCanvas("c2", "Canvas", 800, 600); // secondo canvas per il secondo istogramma
    hist13->GetXaxis()->SetTitle("#Delta t_1_3 [ns]"); // secondo istogramma: differenza temporale tra canali 1-3 (1-3)
    hist13->GetYaxis()->SetTitle("Counts");
    hist13->SetFillColor(kRed-7);
    hist13->Draw();
    TCanvas *c3 = new TCanvas("c3", "Canvas", 800, 600); // terzo canvas per il terzo istogramma
    hist23->GetXaxis()->SetTitle("#Delta t_2_3 [ns]"); // terzo istogramma: differenza temporale tra canali 2-3 (2-3)
    hist23->GetYaxis()->SetTitle("Counts");
    hist23->SetFillColor(kGreen-7); 
    hist23->Draw();     
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

    // Aggiorna il display
    DrawEvent(gCurrentEvent);
}
// DRS4Browser_v3(): funzione principale. ROOT la chiama automaticamente
// quando si lancia:   root -l 'DRS4Browser_v3.cpp("file.xml")'
//
// Il nome DEVE corrispondere al nome del file senza estensione!
//funzione derivata da DRS4Browser_v3() 
//serve per mostrare la risoluzione temporale dei canali, usando CFD e walk-corrected leading edge
//realizza anche il grafico delle differenze di tempo tra i piani dello stesso telescopio in tre cannvas separati
void Resolution(const char* xmlfile, const char* labels = "") {

    std::cout << "=============================================" << std::endl;
    std::cout << "  DRS4 Waveform Browser v3                   " << std::endl;
    std::cout << "  Timing: CFD + Walk-Corrected Leading Edge  " << std::endl;
    std::cout << "  Laboratorio Sciami Estesi, Marzo 2026      " << std::endl;
    std::cout << "=============================================" << std::endl;

    gFilename = xmlfile;
    gEvents.clear();     // Svuota eventuali dati da un lancio precedente

    // --- Parsing delle etichette personalizzate ---
    // Formato atteso: "CHN1=PMT top,CHN2=PMT mid,CHN3=PMT bot,CHN4=NIM trigger"
    gChannelLabels[1] = "PMT1";    // Default
    gChannelLabels[2] = "PMT2";
    gChannelLabels[3] = "PMT3";
    gChannelLabels[4] = "TRG_tripl_nim";     // Default per CH4 (trigger)

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
