#from data.loading_data import *
#import utils_plotly  # don't delete; add plotly presets
import pandas as pd 
import numpy as np
#import plotly.graph_objects as go



def correlation_analysis(X: pd.DataFrame, threshold: float = 0.98, output: bool = False):
    feature_labels = list(X.columns)
    numFeatures = X.shape[1]

    # Pearson correlation matrix
    correlationMatrix = np.corrcoef(X.values, rowvar=False)
    pearsonDF = pd.DataFrame(correlationMatrix, index=feature_labels, columns=feature_labels)

    group_features = dict()
    drop_set = set()

    # Identify groups of highly correlated features
    for i in range(numFeatures):
        feature = feature_labels[i]
        if i in drop_set:
            continue

        group_features[feature] = []
        for j in range(i + 1, numFeatures):
            if abs(correlationMatrix[i, j]) > threshold:
                drop_set.add(j)
                group_features[feature].append(feature_labels[j])

    drop_features = [feature_labels[i] for i in drop_set]

    X_reduced = X.drop(columns=drop_features)

    if output:
        print("Features removed due to correlation:", len(drop_set))
        if drop_features:
            print("\t" + ", ".join(drop_features))

    return X_reduced, list(group_features.keys()), group_features, pearsonDF


def spearmanr_correlation(X: pd.DataFrame, threshold: float = 0.98, output: bool = False):
    from scipy import stats
    feature_labels = list(X.columns)
    """
    The Spearman rank-order correlation coefficient is a nonparametric measure of the monotonicity of the relationship
    between two datasets. Like other correlation coefficients, this one varies between -1 and +1 with 0 implying
    no correlation. Correlations of -1 or +1 imply an exact monotonic relationship. Positive correlations imply
    that as x increases, so does y. Negative correlations imply that as x increases, y decreases.
    """
    numFeatures = X.shape[1]

    # Compute correlation matrix
    result = stats.spearmanr(X.values, axis=0)
    correlationMatrix = result.correlation
    spearmanDF = pd.DataFrame(correlationMatrix, index=feature_labels, columns=feature_labels)

    # Identify features with high correlation
    group_features = dict()
    drop_set = set()
    for i in range(numFeatures):
        feature = feature_labels[i]
        if i in drop_set:
            continue

        group_features[feature] = []
        for j in range(i + 1, numFeatures):
            if abs(correlationMatrix[i, j]) > threshold:
                drop_set.add(j)
                group_features[feature].append(feature_labels[j])

    # Features to drop
    drop_features = [feature_labels[i] for i in drop_set]

    # Keep as DataFrame
    X_reduced = X.drop(columns=drop_features)

    if output:
        print("Features removed due to spearmanr_correlation:", len(drop_set))
        print("\t" + ", ".join(drop_features))

    return X_reduced, drop_features, spearmanDF


def remove_by_variance(X,feature_labels: list[str],threshold: float = 0) -> tuple[pd.DataFrame, list[str], list[str]]:

    from sklearn.feature_selection import VarianceThreshold
    import numpy as np
    import pandas as pd

    sel = VarianceThreshold(threshold=threshold)

    # Preserve index if X is a DataFrame
    if isinstance(X, pd.DataFrame):
        index = X.index
        X_values = X.values
    else:
        index = None
        X_values = X

    X_reduced = sel.fit_transform(X_values)
    keep = list(sel.get_feature_names_out(feature_labels))
    drop = [label for label in feature_labels if label not in keep]

    print(f"Features removed due to low variance <{threshold}: {len(drop)}")
    if drop:
        print("\t" + ", ".join(drop))

    XDF = pd.DataFrame(X_reduced, columns=keep, index=index)

    return XDF, keep, drop


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

