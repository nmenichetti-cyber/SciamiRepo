info su come lanciare i codici per l'analisi dati della DRS.

#############################
#                           #
#     DRS4Browser_v4.cpp    #
#                           #
#############################

=====================================================
Questo codice fa il parsing dei dati della DRS4 in modo da inserirli all'interno di una struttura 

-ChannelData: un singolo canale di un singolo evento
-EventData:   un evento completo con tutti i suoi canali

Una volta effettuato il parsing,a seconda che sia stato rilevato clipping o meno del segnale campionat vengono utliizzati due meccanismi differenti di ricerca
dell'istante temporale in cui arriva l'impulso:

-nel caso no clipping si utilizza il metodo di costant fraction discrimination e si ricerca il momento in cui il segnale raggiunge metà dell'ampiezza 
-nel caso ci sia clipping si utlizza una soglia costante e si corregge per il walk time usando la derivata del segnale nella discesa e calcolando lo slewrate massimo
 
 Una volta ricavati gli istanti temporali si calcolano tutti gli intervalli tra i vari telescopi 08-06 08-04 06-04



Caricare su root con .L DRS4Browser_v4.cpp
alternativamente con  .L DRS4Browser_v4.cpp+

La funzione principale è DRS4Browser("/percorso_dati/dati.xml", "custom labels")

l'opzione custom labels è praticamente inutile si può lanciare semplicemente con RS4Browser_v4("/percorso_dati/dati.xml")

Questo aprirà il browser per visualizzare le 4 forme d'onda campionate evento per evento,
nel canvas sono riportati i delta_t calcolati e il metodo con cui sono stati ricavati C (constant fraction)   W (threshold-walk-time corrected) 

Una volta aperto il browser è possibile utilizzare una serie di comandi da terminale per interagire:

  Navigazione:  Next()  Prev() scorri gli eventi
                GoTo(n) vai all'evento desiderato
                GoToSerial(s)   



  Overlay:      Overlay()  Overlay(true)  OverlayNorm() vari tipi di visualizzazione, scala normale uguale per tutti o rinormalizzazione per avere una visualizzazione migliore dei segnali 
                OverlaySelect("1,2,3")  OverlaySelect("1,2,3",true,true)    scelta dei canali da visualizzare in un univo plot 



  Angoli:       Skymap()                 distribuzioni angolari in dtheta dcostheta piano (l,m) 

  Calibrazione: CalibrateCableOffsets()   stima offset dalla media dei dt in maniera da correggere errori sulla misura grossolana del ritardo
                SetCableOffsets(a,b)      imposta offset manuali [ns]


  Altro:        Save()  SaveAll()  salva il canvas corrente in pdf o salva tutti i canvas in un unico pdf
                Summary()   riepilogo di tutte le info 

Prima è necessario lanciare la funzione CalibrateCableOffsets questa corregge i ritardi introdotto dai cavi facendo la media dei dt,
i dt =  dt_sciame (dovuto al diverso tempo di arrivo ordine di pochi nanosecondi 30cm/ns ) dt_cavi (ordine decine di nanosecondi 206 per lo 06 e 67 per il 04)
in media si troverà il valore del ritardo dovuto alla linea di trasmissione, in questo modo è possibile correggere e trovare la vera distribuzione angolare.


Una volta fatta la calibrazione dei ritardi dei cavi si può lanciare l'opzione Skymap() 
questa genera 4 grafici per le distribuzioni in angolo degli sciami

-la prima la proiezione zenitale (l,m) m=sin(theta)cos(phi) l=sin(theta)sin(phi)
-la distribuzione dN/dtheta proporzionale  a sin(theta)cos(theta)^nanosecondi
-la distribuzione dN/dcos(theta) proporzionale a cos(theta)^n 
-la distribuzione dN/dphi che dovrebbe essere circa piatta 


#############################
#                           #
#     Resolution.cpp        #
#                           #
#############################

=====================================================

Il codice ha come scopo principale ricavare la distribuzione dei ritardi tra i piani di un telescopio in maniera da indagare sulla risoluzione intrinseca di ognuno.
il programma dopo aver fatto il parsing dei dati li carica nelle stesse strutture dati del precedente programma fa ricerca del clipping 
e calcola in maniera analoga i tempi di arrivo e gli intervalli temporali tra i diversi piani del telescopio.
In seguito realizza un istogramma con le occorrenze degli intervalli.

Per lanciare il programma:

Caricare su root con .L Resolution.cpp
alternativamente con  .L Resolution.cpp+

La funzione principale è Resolution("/percorso_dati/dati.xml", "custom labels")
Questo aprirà il browser per visualizzare le 4 forme d'onda campionate evento per evento,
nel canvas sono riportati i delta_t calcolati e il metodo con cui sono stati ricavati C (constant fraction)   W (threshold-walk-time corrected)

Per visualizzare l'istogramma delle differenze di tempo è necessario lanciare la funzione 

Hist_resolution();

questa funzione genera 3 canvas con i grafici desiderati.