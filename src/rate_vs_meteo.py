#!/usr/bin/env python3
"""
Stem plot del rate orario degli eventi XML sovrapposto alla media oraria
di variabili meteo da file stazione Davis (downld*.txt).

Uso:
    python3 rate_vs_meteo.py events.xml weather.txt --start "24/03/26 16:00"
    python3 rate_vs_meteo.py events.xml weather.txt --start "24/03/26 16:00" -o plots/
"""

import re
import sys
import argparse
from collections import Counter, defaultdict
from pathlib import Path
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np


# -----------------------------------------------------------------------
# Colonne nel file meteo (0-based, dopo split())
# -----------------------------------------------------------------------
COL_DATE         = 0
COL_TIME         = 1
COL_TEMP_OUT     = 2
COL_BAR          = 16
COL_SOLAR_ENERGY = 20
COL_AIR_DENSITY  = 32
COL_HUM_OUT      = 5

DATE_FMT_METEO  = "%d/%m/%y %H:%M"
DATE_FMT_EVENTS = "%Y/%m/%d %H:%M:%S.%f"


# -----------------------------------------------------------------------
# Parsing
# -----------------------------------------------------------------------

def parse_start(start_str: str) -> datetime:
    """Converte la stringa di start in datetime. Accetta 'DD/MM/YY HH:MM'."""
    for fmt in ("%d/%m/%y %H:%M", "%Y-%m-%d %H:%M", "%d/%m/%Y %H:%M"):
        try:
            return datetime.strptime(start_str, fmt)
        except ValueError:
            pass
    raise ValueError(f"Formato data non riconosciuto: '{start_str}'. Usa DD/MM/YY HH:MM")


def parse_events(xml_path: str, start: datetime) -> dict:
    """Restituisce Counter {ora: n_eventi} per eventi >= start."""
    pattern = re.compile(
        r"<Time>(\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}\.\d+)</Time>"
    )
    counts = Counter()
    with open(xml_path, "r", encoding="ISO-8859-1") as f:
        for line in f:
            m = pattern.search(line)
            if m:
                dt = datetime.strptime(m.group(1), DATE_FMT_EVENTS)
                if dt >= start:
                    counts[dt.replace(minute=0, second=0, microsecond=0)] += 1
    return counts


def parse_meteo(txt_path: str, start: datetime) -> dict:
    """
    Restituisce dizionario:
      { ora: {'temp': [...], 'bar': [...], 'solar': [...], 'density': [...]} }
    per righe con timestamp >= start.
    """
    data = defaultdict(lambda: {"temp": [], "bar": [], "solar": [], "density": [], "hum": []})

    with open(txt_path, "r", encoding="latin-1") as f:
        for line in f:
            line = line.rstrip("\r\n")
            cols = line.split()
            # Salta header e separatore
            if len(cols) < 33:
                continue
            try:
                dt = datetime.strptime(cols[COL_DATE] + " " + cols[COL_TIME],
                                       DATE_FMT_METEO)
            except ValueError:
                continue
            if dt < start:
                continue

            hour = dt.replace(minute=0, second=0, microsecond=0)
            try:
                data[hour]["temp"].append(float(cols[COL_TEMP_OUT]))
                data[hour]["bar"].append(float(cols[COL_BAR]))
                data[hour]["solar"].append(float(cols[COL_SOLAR_ENERGY]))
                data[hour]["density"].append(float(cols[COL_AIR_DENSITY]))
                data[hour]["hum"].append(float(cols[COL_HUM_OUT]))
            except (ValueError, IndexError):
                continue

    # Calcola medie orarie
    averages = {}
    for hour, vals in data.items():
        averages[hour] = {k: float(np.mean(v)) if v else float("nan")
                          for k, v in vals.items()}
    return averages


# -----------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------

METEO_VARS = [
    ("temp",    "Temperatura esterna (°C)",  "#e74c3c"),
    ("bar",     "Pressione (hPa)",            "#2980b9"),
    ("solar",   "Energia solare (MJ/m²)",     "#f39c12"),
    ("density", "Densità aria (kg/m³)",       "#27ae60"),
    ("hum",     "Umidità esterna (%)",   "#8e44ad"),
]


def stem_ax(ax, hours, counts, color="#555555", label="Rate [ev/h]"):
    """Disegna uno stem plot sull'asse ax."""
    if not hours:
        return
    for x, y in zip(hours, counts):
        ax.plot([x, x], [0, y], color=color, linewidth=2)
    ax.plot(hours, counts, marker="o", markersize=7, color=color,
            linestyle="none", label=label, zorder=5)
    ax.set_ylim(bottom=0, top=max(counts) * 1.3 if counts else 1)


def make_plot(event_hours, event_counts, meteo_hours, meteo_vals,
              var_label, var_color, title, outfile=None):
    """Crea un singolo plot con asse sinistro (rate) e destro (variabile meteo)."""

    import datetime as _dt

    # --- Taglia i dati meteo prima del primo evento ---
    first_event = min(event_hours) if event_hours else None
    if first_event and meteo_hours:
        paired      = [(h, v) for h, v in zip(meteo_hours, meteo_vals) if h >= first_event]
        meteo_hours = [h for h, v in paired]
        meteo_vals  = [v for h, v in paired]

    # ---------------------------------------------------------------
    # STILE PLOT — modifica qui per cambiare l'aspetto dei grafici:
    #   fig size    → riga "figsize"          (larghezza, altezza in pollici)
    #   colore stem → color="#444444"         in stem_ax(...)
    #   colore meteo→ var_color               passato da METEO_VARS
    #   marker stem → markersize=7            in stem_ax(...)
    #   marker meteo→ markersize=5, marker="s" nella ax2.plot(...)
    #   titolo      → fig.suptitle(...)
    #   dpi output  → plt.savefig(..., dpi=150)
    # ---------------------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(13, 5))

    # --- Stem plot rate (asse sinistro) ---
    stem_ax(ax1, event_hours, event_counts,
            color="#444444", label="Rate orario [ev/h]")
    ax1.set_xlabel("Data e ora", fontsize=12)
    ax1.set_ylabel("Rate orario [ev/h]", fontsize=12, color="#444444")
    ax1.tick_params(axis="y", labelcolor="#444444")

    # --- Asse X: etichette orizzontali ogni 6 ore a partire dal primo evento ---
    if first_event:
        x_max = max(event_hours + (meteo_hours if meteo_hours else []))
        ax1.set_xlim(first_event - _dt.timedelta(minutes=30),
                     x_max      + _dt.timedelta(minutes=30))
        # Genera i tick ogni 6 ore a partire esattamente da first_event
        ticks = []
        t = first_event
        while t <= x_max + _dt.timedelta(hours=6):
            ticks.append(t)
            t += _dt.timedelta(hours=6)
        ax1.set_xticks(ticks)

    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M\n%d/%m"))
    plt.xticks(rotation=0, ha="center")   # etichette orizzontali
    ax1.grid(True, linestyle="--", alpha=0.3)

    # --- Media meteo (asse destro) ---
    ax2 = ax1.twinx()
    if meteo_hours:
        ax2.plot(meteo_hours, meteo_vals, color=var_color, linewidth=0,
                 marker="o", markersize=5, label=var_label, zorder=4)
        ax2.set_ylabel(var_label, fontsize=12, color=var_color)
        ax2.tick_params(axis="y", labelcolor=var_color)
        margin = (max(meteo_vals) - min(meteo_vals)) * 0.3 or 1
        ax2.set_ylim(min(meteo_vals) - margin, max(meteo_vals) + margin)

    # Legenda combinata
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               loc="upper left", fontsize=9)

    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.01)
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        print(f"  Salvato: {outfile}")
    else:
        plt.show()

    plt.close(fig)

def make_scatter(event_hours, event_counts, meteo_hours, meteo_vals,
                 var_label, var_color, title, outfile=None):
    """Scatter plot: rate orario [ev/h] vs variabile meteo, con colorbar temporale."""

    # Trova i timestamp in comune
    ev_dict  = dict(zip(event_hours, event_counts))
    met_dict = dict(zip(meteo_hours, meteo_vals))
    common   = sorted(set(ev_dict) & set(met_dict))

    if not common:
        print(f"  [scatter] Nessuna ora in comune per '{var_label}', skip.")
        return

    x = [met_dict[h] for h in common]   # variabile meteo
    y = [ev_dict[h]  for h in common]   # rate orario

    # Colore per punto = posizione temporale (0 → primo, 1 → ultimo)
    t_num = mdates.date2num(common)
    norm  = plt.Normalize(t_num.min(), t_num.max())
    cmap  = plt.cm.plasma

    fig, ax = plt.subplots(figsize=(7, 5))
    sc = ax.scatter(x, y, c=t_num, cmap=cmap, norm=norm,
                    s=70, edgecolors="white", linewidths=0.5, zorder=5)

    # Etichette ora su ogni punto
    # for xi, yi, hi in zip(x, y, common):
    #     ax.annotate(hi.strftime("%H:%M\n%d/%m"), (xi, yi),
    #                 textcoords="offset points", xytext=(5, 4),
    #                 fontsize=6.5, color="#333333")

    # Retta di regressione lineare
    # if len(x) >= 2:
    #     coeffs = np.polyfit(x, y, 1)
    #     x_line = np.linspace(min(x), max(x), 200)
    #     ax.plot(x_line, np.polyval(coeffs, x_line),
    #             color=var_color, linewidth=1.5, linestyle="--",
    #             label=f"fit lineare  (m={coeffs[0]:.3g})")
    #     ax.legend(fontsize=9)


    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Tempo", fontsize=10)
    cbar.ax.yaxis.set_major_formatter(
        mdates.DateFormatter("%H:%M\n%d/%m")
    )

    ax.set_xlabel(var_label, fontsize=12, color=var_color)
    ax.set_ylabel("Rate orario [ev/h]", fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.3)
    fig.suptitle(title, fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        print(f"  Salvato: {outfile}")
    else:
        plt.show()

    plt.close(fig)

def make_scatter_3d(event_hours, event_counts, meteo_hours, meteo_vals,
                    var_label, var_color, title, outfile=None):
    """Scatter 3D: X = variabile meteo, Y = rate orario, Z = tempo."""

    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

    # Ore in comune
    ev_dict  = dict(zip(event_hours, event_counts))
    met_dict = dict(zip(meteo_hours, meteo_vals))
    common   = sorted(set(ev_dict) & set(met_dict))

    if not common:
        print(f"  [scatter3d] Nessuna ora in comune per '{var_label}', skip.")
        return

    x = np.array([met_dict[h] for h in common])          # variabile meteo
    y = np.array([ev_dict[h]  for h in common])           # rate orario
    z = np.array(mdates.date2num(common), dtype=float)    # tempo (numerico)

    # Colore = tempo, stessa colormap del 2D
    norm  = plt.Normalize(z.min(), z.max())
    cmap  = plt.cm.plasma
    colors = cmap(norm(z))

    fig = plt.figure(figsize=(10, 7))
    ax  = fig.add_subplot(111, projection="3d")

    sc = ax.scatter(x, y, z, c=z, cmap=cmap, norm=norm,
                    s=80, edgecolors="white", linewidths=0.4, depthshade=True)

    # Linee verticali dal piano XY fino al punto (aiutano a leggere la Z)
    for xi, yi, zi in zip(x, y, z):
        ax.plot([xi, xi], [yi, yi], [z.min(), zi],
                color="grey", linewidth=0.6, alpha=0.4)

    # # Etichette ora su ogni punto
    # for xi, yi, zi, hi in zip(x, y, z, common):
    #     ax.text(xi, yi, zi, hi.strftime("%H:%M\n%d/%m"),
    #             fontsize=6, color="#222222")

    # Colorbar con etichette temporali
    cbar = fig.colorbar(sc, ax=ax, pad=0.12, shrink=0.6)
    cbar.set_label("Tempo", fontsize=10)
    cbar.ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M\n%d/%m"))

    # Asse Z: sostituisce i numeri matplotlib con etichette leggibili
    z_ticks = np.linspace(z.min(), z.max(), min(6, len(common)))
    ax.set_zticks(z_ticks)
    ax.set_zticklabels(
        [mdates.num2date(t).strftime("%H:%M\n%d/%m") for t in z_ticks],
        fontsize=7
    )

    ax.set_xlabel(var_label, fontsize=10, color=var_color, labelpad=10)
    ax.set_ylabel("Rate orario [ev/h]", fontsize=10, labelpad=10)
    ax.set_zlabel("Tempo", fontsize=10, labelpad=10)

    fig.suptitle(title, fontsize=13, fontweight="bold")
    plt.tight_layout()

    if outfile:
        plt.savefig(outfile, dpi=150, bbox_inches="tight")
        print(f"  Salvato: {outfile}")
    else:
        plt.show()

    plt.close(fig)

def make_daily_plots(event_hours, event_counts, meteo_hours, meteo_vals,
                     var_label, var_color, title_prefix, outdir=None):

    ev_dict  = dict(zip(event_hours, event_counts))
    met_dict = dict(zip(meteo_hours, meteo_vals))
    common   = sorted(set(ev_dict) & set(met_dict))

    if not common:
        print(f"  [daily] Nessuna ora in comune per '{title_prefix}', skip.")
        return

    by_day = defaultdict(list)
    for h in common:
        by_day[h.date()].append(h)

    days = sorted(by_day.keys())

    for idx, day in enumerate(days):
        hours = by_day[day]
        x = [ev_dict[h]  for h in hours]
        y = [met_dict[h] for h in hours]

        t_num = mdates.date2num(hours)
        norm  = plt.Normalize(t_num.min(), t_num.max())
        cmap  = plt.cm.plasma

        def _fill_ax(ax):                                   # ← ANNIDATA QUI
            sc = ax.scatter(x, y, c=t_num, cmap=cmap, norm=norm,
                            s=70, edgecolors="white", linewidths=0.5, zorder=5)
            # for xi, yi, hi in zip(x, y, hours):
            #     ax.annotate(hi.strftime("%H:%M"), (xi, yi),
            #                 textcoords="offset points", xytext=(5, 4),
            #                 fontsize=7, color="#333333")
            if len(x) >= 2:
                r = float(np.corrcoef(x, y)[0, 1])
                ax.text(1.02, 0.10, f"Pearson r = {r:.2f}",
                        transform=ax.transAxes,
                        fontsize=10, fontweight="bold",
                        verticalalignment="top",
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7))
            ax.set_xlabel("Rate orario [ev/h]", fontsize=10)
            ax.set_ylabel(var_label, fontsize=10, color=var_color)
            ax.grid(True, linestyle="--", alpha=0.3)
            ax.set_title(day.strftime("%d/%m/%Y"), fontsize=11, fontweight="bold")
            return sc

        # --- Plot singolo ---
        fig_single, ax_single = plt.subplots(figsize=(7, 5))
        sc = _fill_ax(ax_single)
        cbar = fig_single.colorbar(sc, ax=ax_single, pad=0.12, shrink=0.6)
        cbar.set_label("Tempo", fontsize=10)
        cbar.ax.yaxis.set_major_formatter(mdates.DateFormatter("%H:%M\n%d/%m"))
        fig_single.suptitle(f"{title_prefix}  vs  Rate  —  {day.strftime('%d/%m/%Y')}",
                            fontsize=13, fontweight="bold", y=1.01)
        plt.tight_layout()
        if outdir:
            fname = outdir / f"daily_{day.strftime('%Y%m%d')}_{title_prefix}.png"
            plt.savefig(fname, dpi=150, bbox_inches="tight")
            print(f"  Salvato: {fname}")
        else:
            plt.show()
        plt.close(fig_single)
        
# -----------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Stem plot rate orario eventi + variabili meteo."
    )
    parser.add_argument("xml",     help="File XML degli eventi (DRSOSC)")
    parser.add_argument("meteo",   help="File TXT della stazione meteo")
    parser.add_argument("--start", required=True,
                        help="Data/ora di inizio: 'DD/MM/YY HH:MM'")
    parser.add_argument("-o", "--outdir",
                        help="Cartella di output per i PNG (default: mostra a schermo)",
                        default=None)
    args = parser.parse_args()

    for p in [args.xml, args.meteo]:
        if not Path(p).exists():
            print(f"Errore: file non trovato: {p}", file=sys.stderr)
            sys.exit(1)

    try:
        start = parse_start(args.start)
    except ValueError as e:
        print(f"Errore: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Data/ora di inizio: {start.strftime('%d/%m/%Y %H:%M')}")

    # --- Parsing ---
    print("Lettura eventi XML...")
    ev_counts = parse_events(args.xml, start)
    if not ev_counts:
        print("Nessun evento trovato dopo la data di inizio.", file=sys.stderr)
        sys.exit(1)

    print(f"  {sum(ev_counts.values())} eventi in {len(ev_counts)} ore.")

    print("Lettura dati meteo...")
    meteo = parse_meteo(args.meteo, start)
    print(f"  {len(meteo)} ore di dati meteo.")

    # Ore eventi ordinate
    ev_hours = sorted(ev_counts.keys())
    ev_vals  = [ev_counts[h] for h in ev_hours]

    # Ore meteo: usa solo quelle presenti nei dati meteo
    all_hours = sorted(meteo.keys())

    # Cartella output
    outdir = Path(args.outdir) if args.outdir else None
    if outdir:
        outdir.mkdir(parents=True, exist_ok=True)

    # --- 4 plot ---
    print("Generazione plot...")
    for var_key, var_label, var_color in METEO_VARS:
        met_hours = [h for h in all_hours if not np.isnan(meteo[h][var_key])]
        met_vals  = [meteo[h][var_key] for h in met_hours]

        # title = f"Rate orario  vs  {var_label}"
        # outfile = str(outdir / f"rate_vs_{var_key}.png") if outdir else None
        # make_plot(ev_hours, ev_vals, met_hours, met_vals, var_label, var_color, title, outfile=outfile)
        
        # scatter_title   = f"Rate orario  vs  {var_label}  [scatter]"
        # scatter_outfile = str(outdir / f"scatter_vs_{var_key}.png") if outdir else None
        # make_scatter(ev_hours, ev_vals, met_hours, met_vals, var_label, var_color, scatter_title, outfile=scatter_outfile)

        make_daily_plots(ev_hours, ev_vals, met_hours, met_vals, var_label, var_color, var_key, outdir=outdir)
        
        # scatter3d_title   = f"Rate orario  vs  {var_label}  [3D]"
        # scatter3d_outfile = str(outdir / f"scatter3d_vs_{var_key}.png") if outdir else None
        # make_scatter_3d(ev_hours, ev_vals, met_hours, met_vals, var_label, var_color, scatter3d_title, outfile=scatter3d_outfile)

    print("Fatto.")


if __name__ == "__main__":
    main()