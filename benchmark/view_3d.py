"""Create a 3D view of matched C-alpha atoms after the global fit."""

import argparse
from pathlib import Path

import plotly.graph_objects as go

from .snapshot import ROOT, verify_snapshot
from .structure import compare


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "results" / "pgk1_superposition_3d.html")
    args = parser.parse_args()
    verify_snapshot()
    comparison = compare()
    figure = go.Figure()

    for name, coordinates, color in (
        ("AlphaFold", comparison.global_fit, "#116d83"),
        ("3PGK", comparison.experimental, "#d66c32"),
    ):
        x, y, z = [], [], []
        previous = None
        for position, coordinate in zip(comparison.positions, coordinates, strict=True):
            if previous is not None and position != previous + 1:
                x.append(None)
                y.append(None)
                z.append(None)
            x.append(coordinate[0])
            y.append(coordinate[1])
            z.append(coordinate[2])
            previous = position
        figure.add_trace(go.Scatter3d(
            x=x, y=y, z=z, mode="lines", name=name,
            line={"color": color, "width": 4}, hoverinfo="name",
        ))
        mask = (comparison.positions >= 248) & (comparison.positions <= 258)
        region = coordinates[mask]
        positions = comparison.positions[mask]
        labels = [f"{name} · {comparison.canonical[int(i) - 1]}{i}" for i in positions]
        sizes = [9 if i == 254 else 4 for i in positions]
        figure.add_trace(go.Scatter3d(
            x=region[:, 0], y=region[:, 1], z=region[:, 2],
            mode="markers", name=f"{name} · residues 248–258",
            marker={"color": color, "size": sizes},
            text=labels, hovertemplate="%{text}<extra></extra>",
        ))

    figure.update_layout(
        title="PGK1: AlphaFold and 3PGK after global superposition",
        scene={"aspectmode": "data", "xaxis_title": "X (Å)", "yaxis_title": "Y (Å)", "zaxis_title": "Z (Å)"},
        margin={"l": 0, "r": 0, "t": 55, "b": 0},
        legend={"orientation": "h", "y": 0.02},
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(args.output, include_plotlyjs=True, full_html=True)
    print(f"Saved 3D view: {args.output}")
    print("Blue: AlphaFold. Orange: 3PGK. Large markers: S254.")


if __name__ == "__main__":
    main()
