#from data.loading_data import *
#import utils_plotly  # don't delete; add plotly presets
import pandas as pd 
import numpy as np
#import plotly.graph_objects as go


def correlation_analysis(X, feature_labels: list[str], threshold: float = 0.98, output: bool = False): #X is a Pandas array
    
    numFeatures = X.shape[1]
    scaleWindow = 500 * numFeatures / 6

    correlation_matrix = np.corrcoef(X, rowvar=False)

    # Identify features with high correlation and deleting them
    # TODO: need to fix multiple correlation
    group_features = dict()
    drop_set = set()
    for i in range(numFeatures):
        feature = feature_labels[i]
        if i in drop_set:
            continue

        group_features[feature] = []
        for j in range(i + 1, numFeatures):
            if abs(correlation_matrix[i, j]) > threshold:
                drop_set.add(j)
                group_features[feature].append(feature_labels[j])

    mask = np.zeros(numFeatures, dtype=bool)
    mask[list(drop_set)] = True
    X = np.delete(X, mask, axis=1)

    print("Features removed due to correlation:", len(drop_set))
    drop_features = [label for i, label in enumerate(feature_labels) if i in drop_set]
    print("\t" + ", ".join(drop_features))

    return X, list(group_features.keys()), group_features


def spearmanr_correlation(X: np.ndarray, feature_labels: list[str], threshold: float = 0.98, output: bool = False):
    """
    The Spearman rank-order correlation coefficient is a nonparametric measure of the monotonicity of the relationship
    between two datasets. Like other correlation coefficients, this one varies between -1 and +1 with 0 implying
    no correlation. Correlations of -1 or +1 imply an exact monotonic relationship. Positive correlations imply
    that as x increases, so does y. Negative correlations imply that as x increases, y decreases.
    """
    numFeatures = X.shape[1]
    scaleWindow = 500 * numFeatures / 6

    from scipy import stats
    result = stats.spearmanr(X, axis=0)
    correlation_matrix = result.correlation


    # Identify features with high correlation and deleting them
    # TODO: need to fix multiple correlation
    group_features = dict()
    drop_set = set()
    for i in range(numFeatures):
        feature = feature_labels[i]
        if i in drop_set:
            continue

        group_features[feature] = []
        for j in range(i + 1, numFeatures):
            if abs(correlation_matrix[i, j]) > threshold:
                drop_set.add(j)
                group_features[feature].append(feature_labels[j])

    mask = np.zeros(numFeatures, dtype=bool)
    mask[list(drop_set)] = True
    X = np.delete(X, mask, axis=1)

    print("Features removed due to spearmanr_correlation:", len(drop_set))
    drop_features = [label for i, label in enumerate(feature_labels) if i in drop_set]
    print("\t" + ", ".join(drop_features))

    xMAST = pd.DataFrame(X, columns=group_features)

    return xMAST , drop_features


def remove_by_variance(X: np.ndarray, feature_labels: list[str], threshold: float = 0) -> tuple[np.ndarray, list[str], list[str]]:
    from sklearn.feature_selection import VarianceThreshold
    sel = VarianceThreshold(threshold=threshold)
    X = sel.fit_transform(X)
    keep = sel.get_feature_names_out(feature_labels)
    drop = [label for label in feature_labels if label not in keep]

    print(f"Features removed due to low variance <{threshold}: {len(drop)}")
    print("\t" + ", ".join(drop))
    return X, keep, drop


def main():
    feature_label = "mordred"  # "dft" "mordred" "rdkit"   <<<<<<
    feature_labels, features = load_features(feature_label)
    smiles, methods, yields = load_yield()

    X = features
    Y = yields

    print("Initial features:", len(feature_labels))
    with open(get_parent_folder() / 'data' / f'features_{feature_label}_removed.txt', 'w') as f:
        f.write(f"Total starting feature count: {len(feature_labels)}")
        f.write("".join([f'\n\t{label}' for label in feature_labels]))
        
        X, feature_labels, dropped_features = remove_by_variance(X, feature_labels)
        text = "\n\n\nFeatures drop due to low variance: " + "".join([f'\n\t{label}' for label in dropped_features])
        f.write(text)

        X, feature_labels, drop_group = correlation_analysis(X, feature_labels, threshold=0.95)
        text = "\n\n\nFeatures drop due correlation: (STILL HAS ISSUES) "
        import json
        text += json.dumps(drop_group, indent=4).replace('\n', '\n\t')
        f.write(text)

        X, feature_labels, drop_group = spearmanr_correlation(X, feature_labels, threshold=0.95)
        text = "\n\n\nFeatures drop due spearmanr_correlation: (STILL HAS ISSUES) "
        import json
        text += json.dumps(drop_group, indent=4).replace('\n', '\n\t')
        f.write(text)

    print("remaining features:", len(feature_labels))
    with open(get_parent_folder() / 'data' / f'features_{feature_label}_filter.csv', 'w', newline='') as f:
        csv_writer = csv.writer(f, delimiter=',')
        csv_writer.writerow(["smiles"] + feature_labels)
        csv_writer.writerows([[smiles[i]] + d for i, d in enumerate(X.tolist())])


if __name__ == "__main__":
    main()

