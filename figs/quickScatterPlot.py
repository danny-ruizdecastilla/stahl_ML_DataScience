# scatter plot generator
# Daniel Ruiz de Castilla | 05.15.2026

import plotly.graph_objects as go
import pandas as pd
from pathlib import Path
import sys

parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))

def scatterTemplate():
    template = go.layout.Template()

    template.layout.font = dict(
        family="Arial",
        size=18,
        color="black"
    )

    template.layout.plot_bgcolor = "white"
    template.layout.width = 1200
    template.layout.height = 600

    # x-axis formatting
    template.layout.xaxis.showline = True
    template.layout.xaxis.linewidth = 5
    template.layout.xaxis.linecolor = "black"
    template.layout.xaxis.showgrid = False
    template.layout.xaxis.mirror = True
    template.layout.xaxis.ticks = "outside"
    template.layout.xaxis.tickwidth = 4

    # y-axis formatting
    template.layout.yaxis.showline = True
    template.layout.yaxis.linewidth = 5
    template.layout.yaxis.linecolor = "black"
    template.layout.yaxis.showgrid = False
    template.layout.yaxis.mirror = True
    template.layout.yaxis.ticks = "outside"
    template.layout.yaxis.tickwidth = 4
    template.layout.xaxis.zeroline = False
    template.layout.yaxis.zeroline = False

    return template


def askingForAddition(askStr):

    response = input(f"{askStr}").strip().lower()

    if response in {"yes", "y"}:
        return True

    elif response in {"no", "n"}:
        return False

    else:
        print("Please enter a valid 'yes' or 'no' response.")
        return askingForAddition(askStr)


def checkCols(df, xCol=None, yCol=None):

    print("\nAvailable columns:")
    print(", ".join(df.columns))

    xInput = input("\nEnter x-column name: ").strip()
    yInput = input("Enter y-column name: ").strip()

    # validate existence
    if xInput not in df.columns:
        print(f"\n'{xInput}' is not a valid column.")
        return checkCols(df, xCol, yCol)

    if yInput not in df.columns:
        print(f"\n'{yInput}' is not a valid column.")
        return checkCols(df, xCol, yCol)

    # validate consistency with previous dataframes
    if xCol is not None and xInput != xCol:
        print(f"\nAll datasets must use the same x-column.")
        print(f"Expected: {xCol}")
        return checkCols(df, xCol, yCol)

    if yCol is not None and yInput != yCol:
        print(f"\nAll datasets must use the same y-column.")
        print(f"Expected: {yCol}")
        return checkCols(df, xCol, yCol)

    return xInput, yInput


def compareSources(fig, saveDir, saveStr, xCol=None, yCol=None):

    addSources = askingForAddition(
        "Do you want to add a dataframe to the scatter plot? (yes/no): "
    )

    if addSources:

        dfDir = input("\nEnter dataframe csv path: ").strip()

        try:
            df = pd.read_csv(dfDir)

        except Exception as e:
            print(f"\nError reading dataframe:\n{e}")
            return compareSources(fig, saveDir, saveStr, xCol, yCol)

        # get validated columns
        xColNew, yColNew = checkCols(df, xCol, yCol)

        # set master columns if first dataframe
        if xCol is None:
            xCol = xColNew

        if yCol is None:
            yCol = yColNew

        datasetName = input(
            "\nEnter label for this dataframe: "
        ).strip()

        # scatter plot
        fig.add_trace(
            go.Scatter(
                x=df[xCol],
                y=df[yCol],
                mode="markers",
                name=datasetName,
                marker=dict(
                    size=12,
                    opacity= 0.8,
                    line=dict(
                        width=1,
                        color="black"
                    )
                )
            )
        )

        fig.update_layout(
            template=scatterTemplate()
        )

        return compareSources(
            fig,
            saveDir,
            saveStr,
            xCol,
            yCol
        )

    else:

        titleStr = input("\nEnter title for scatter plot: ").strip()

        fig.update_layout(
            title=titleStr,
            xaxis_title=xCol,
            yaxis_title=yCol
        )

        fig.write_html(
            saveDir + f"/{saveStr}_Scatter.html"
        )

        print(
            f"\nScatter plot saved to:\n"
            f"{saveDir}/{saveStr}_Scatter.html"
        )


if __name__ == "__main__":

    saveDir = str(sys.argv[1])
    saveStr = str(sys.argv[2])

    savePath = Path(saveDir)
    savePath.mkdir(parents=True, exist_ok=True)

    fig = go.Figure()

    compareSources(fig, saveDir, saveStr)