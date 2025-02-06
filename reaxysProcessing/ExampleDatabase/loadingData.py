import csv
import pathlib
import os

import bigsmiles
import numpy as np
import plotly.graph_objects as go

import chemdraw


def get_parent_folder() -> pathlib.Path:
    return pathlib.Path(__file__).parent.absolute().parent


def save_plotly_image(fig: go.Figure, filename: str, output: bool = True):
    path_ = get_parent_folder() / 'fig' / f'{filename}.png'
    fig.write_image(path_)
    if output:
        print("image saved:", path_)


def filter_features(
        feature_labels: list[str],
        features: np.ndarray,
        remove_features: list[str] = None,
        limit_features: list[str] = None
):
    if ((remove_features is None and limit_features is None) or
            (remove_features is not None and limit_features is not None)):
        raise ValueError("You must specify either remove_features or limit_features, but not both")

    if remove_features is not None:
        remove_index = [i for i, feature in enumerate(feature_labels) if feature in remove_features]
    else:
        remove_index = [i for i, feature in enumerate(feature_labels) if feature not in limit_features]

    feature_labels = [feature_labels[i] for i in range(len(feature_labels)) if i not in remove_index]
    features = np.delete(features, remove_index, axis=1)
    return feature_labels, features


def filter_methods(
        methods_labels: list[str],
        yields: np.ndarray,
        remove_methods: list[str] = None,
        limit_methods: list[str] = None
):
    if ((remove_methods is None and limit_methods is None) or
            (remove_methods is not None and limit_methods is not None)):
        raise ValueError("You must specify either remove_features or limit_features, but not both")

    if remove_methods is not None:
        remove_index = [i for i, feature in enumerate(methods_labels) if feature in remove_methods]
    else:
        remove_index = [i for i, feature in enumerate(methods_labels) if feature not in limit_methods]

    methods_labels = [methods_labels[i] for i in range(len(methods_labels)) if i not in remove_index]
    yields = np.delete(yields, remove_index, axis=1)
    return methods_labels, yields


def categorize_values(array: np.ndarray):
    categorized = np.zeros_like(array, dtype=int)
    categorized[(array >= 15) & (array <= 40)] = 1
    categorized[array > 40] = 2
    return categorized


def _load_features(file_features: str | pathlib.Path):
    with open(file_features, 'r') as f:
        reader = csv.reader(f)
        feature_labels = next(reader)
        features = []
        for row in reader:
            features.append(list(float(i) for i in row[1:]))

    feature_labels = feature_labels[1:]
    features = np.array(features)

    return feature_labels, features


def _load_yield(file_yield: str | pathlib.Path):
    with open(file_yield, 'r') as f:
        reader = csv.reader(f)
        methods = next(reader)
        smiles = []
        yields = []
        for row in reader:
            smiles.append(row[0])
            yields.append(list(float(i) for i in row[1:]))

    methods = methods[1:]
    yields = np.array(yields)

    return smiles, methods, yields


def load_features(key: str = None):
    directory = get_parent_folder() / "data"
    from difflib import SequenceMatcher

    # find file with close match to key
    _file_key = 'features_'
    directory = get_parent_folder() / "data"
    files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f)) and f.startswith(_file_key) and f.endswith('.csv')]
    if not files:
        raise ValueError("No files found in the directory.")

    def similarity(a: str, b: str) -> float:
        return SequenceMatcher(None, a, b).ratio()
    closest_file = max(files, key=lambda file: similarity(file[len(_file_key):], key))

    if similarity(closest_file, _file_key+key) < 0.7:
        raise ValueError(f"File not found.\n\tkey: {key}\n\tclosest_file: {closest_file}")

    return _load_features(directory / closest_file)


def load_yield():
    file_yield = r"C:\Users\nicep\Desktop\pyth_proj\Wisc_ML\wisc_ML\Chlornation_Leah\data\yields.csv"
    return _load_yield(file_yield)


def load_yield_cat():
    file_yield = r"C:\Users\nicep\Desktop\pyth_proj\Wisc_ML\wisc_ML\Chlornation_Leah\data\yields.csv"
    smiles, methods, yields = _load_yield(file_yield)
    return smiles, methods, categorize_values(yields)


def generate_molecule_images(type_: str = "png"):
    smiles, methods, yields = load_yield_cat()

    for i, smile in enumerate(smiles):
        drawer = chemdraw.Drawer(smile)
        drawer.draw_img(get_parent_folder() / "fig"/ "mol" / f"mol_{i}.{type_}")
        print("saved:", i, "/", len(smiles))


def main():
    feature_labels, features = load_features("dft")
    smiles, methods, yields = load_yield()
    methods, yields = filter_methods(methods, yields, remove_methods=['Kajigaeshi', 'Chen', 'Wu', 'Lopez_Stahl', 'Ariafard_Chan', 'Kanai'])
    print(smiles)
    print(methods)
    print(yields)
    print(features)

    print()
    print()
    print("number of molecules:", len(smiles))
    print("number of methods:", len(methods))
    print("number of features:", len(features))
    print("number of points:", yields.size)

    drawer = chemdraw.GridDrawer(smiles)
    drawer.draw_png(get_parent_folder() / "fig" / f"mol_grid.png")

    # check for valid smiles
    smiles_ = [bigsmiles.BigSMILES(smile) for smile in smiles]

    yields = categorize_values(yields)
    unique, counts = np.unique(yields, return_counts=True)
    print()
    print("categorize_values")
    print("category | count")
    for i in range(len(unique)):
        print(unique[i], counts[i])


if __name__ == "__main__":
    main()
    # generate_molecule_svg()
