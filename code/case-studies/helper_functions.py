import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

def draw_cov_ellipse(ax, X, color, n_std=2.0, linewidth=2, **kwargs):
    """
    Draw covariance ellipse for points X (n_samples x 2)
    """
    cov = np.cov(X, rowvar=False)
    mean = X.mean(axis=0)

    # eigen decomposition
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]

    # angle of ellipse
    theta = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * n_std * np.sqrt(vals)

    ellipse = Ellipse(xy=mean, width=width, height=height,
                      edgecolor=color, facecolor='none', 
                      linewidth=linewidth, zorder=3,
                      angle=theta, fill=False, **kwargs)

    ax.add_patch(ellipse)
    return mean

def plot_with_ellipses(res_df, pcs=['PC1','PC2'], meta_col='diagnosis_binned'):
    fig, ax = plt.subplots(figsize=(6,5))

    classes = np.unique(res_df[meta_col])

    colors = ['tab:blue', 'tab:orange']

    for c, color in zip(classes, colors):
        subset = res_df[res_df[meta_col] == c]
        X = subset[pcs].values

        # scatter
        ax.scatter(X[:,0], X[:,1], label=str(c), alpha=0.7)

        # ellipse + centroid
        centroid = draw_cov_ellipse(ax, X, n_std=2, edgecolor=color, linewidth=2)

        # plot centroid
        ax.scatter(centroid[0], centroid[1],
                   color=color, marker='X', s=150, edgecolor='black')

    ax.set_xlabel(pcs[0])
    ax.set_ylabel(pcs[1])
    ax.legend()
    ax.set_title("Latent space with centroids and covariance ellipses")

    return fig, ax