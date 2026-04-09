#!/usr/bin/env python3
"""
Event Counter per file XML DRSOSC
Legge i tag <Time> da ogni evento e produce uno stem plot
con il numero di eventi per ora.

Uso:
    python3 event_counter.py trelescopio_06.xml
    python3 event_counter.py trelescopio_06.xml -o plot.png
"""

import re
import sys
import argparse
from collections import Counter
from pathlib import Path
from datetime import datetime

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np


def parse_timestamps(filepath: str) -> list: 
    timestamps = []
    pattern = re.compile(r"<Time>(\d{4}/\d{2}/\d{2} \d{2}:\d{2}:\d{2}\.\d+)</Time>")
    with open(filepath, "r", encoding="ISO-8859-1") as f:
        for line in f:
            m = pattern.search(line)
            if m:
                dt = datetime.strptime(m.group(1), "%Y/%m/%d %H:%M:%S.%f")
                timestamps.append(dt)
    return timestamps


def count_per_hour(timestamps: list) -> dict:
    hours = [dt.replace(minute=0, second=0, microsecond=0) for dt in timestamps]
    return Counter(hours)


def plot_stem(counts: dict, output_path=None):
    sorted_hours = sorted(counts.keys())
    sorted_counts = [counts[h] for h in sorted_hours]

    fig, ax = plt.subplots(figsize=(12, 5))

    # Stem plot nativo di matplotlib
    markerline, stemlines, baseline = ax.stem(sorted_hours, sorted_counts)

    plt.setp(stemlines,  color="#000000", linewidth=1)
    plt.setp(markerline, color="#0b004d", marker="o", markersize=5, zorder=5)
    plt.setp(baseline,   color="black",   linewidth=1)

    # Asse X con solo ora (HH:00)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
    plt.xticks(rotation=75, ha="right")

    ax.set_xlabel("Ora", fontsize=10)
    ax.set_ylabel("Numero di eventi", fontsize=8)
    ax.set_title("Distribuzione oraria degli eventi", fontsize=8, fontweight="bold")
    ax.set_ylim(bottom=0, top=max(sorted_counts) * 1.25)
    ax.grid(True, linestyle="--", alpha=0.4)

    # Totale eventi
    total = sum(sorted_counts)
    ax.text(0.98, 0.97, f"Totale eventi: {total}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightyellow", alpha=0.8))

    #plt.tight_layout()

    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Plot salvato in: {output_path}")
    else:
        plt.show()

def plot_distribuzione(counts: dict, output_path=None):
    valori = list(counts.values())
    min_v, max_v = min(valori), max(valori)
    bins = range(min_v, max_v + 2)  # +2 per includere l'ultimo bin

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(valori, bins=bins, align="left", color="#2980b9",
            edgecolor="#1a5276", alpha=0.85, rwidth=0.7)

    ax.set_xlabel("Numero di eventi per ora", fontsize=12)
    ax.set_ylabel("Numero di ore", fontsize=12)
    ax.set_title("Distribuzione dei conteggi orari", fontsize=14, fontweight="bold")
    ax.xaxis.set_major_locator(plt.MultipleLocator(1))
    ax.grid(True, linestyle="--", alpha=0.4, axis="y")

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150)
        print(f"Plot salvato in: {output_path}")
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description="Stem plot eventi XML per ora.")
    parser.add_argument("file", help="Percorso del file XML di input")
    parser.add_argument("-o", "--output", help="Salva il plot (es. plot.png)", default=None)
    args = parser.parse_args()

    if not Path(args.file).exists():
        print(f"Errore: file '{args.file}' non trovato.", file=sys.stderr)
        sys.exit(1)

    print(f"Lettura del file: {args.file}")
    timestamps = parse_timestamps(args.file)

    if not timestamps:
        print("Nessun timestamp trovato.", file=sys.stderr)
        sys.exit(1)

    print(f"Trovati {len(timestamps)} eventi totali.")
    counts = count_per_hour(timestamps)

    print("\nEventi per ora:")
    for hour in sorted(counts):
        print(f"  {hour.strftime('%Y-%m-%d %H:00')} -> {counts[hour]} eventi")

    plot_stem(counts, output_path=args.output)

    plot_distribuzione(counts, output_path=args.output.replace(".png", "_dist.png") if args.output else None)




if __name__ == "__main__":
    main()