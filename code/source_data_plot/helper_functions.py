import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib
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
    

omic_title = {'18S': '18S', '16S': '16S', 'mag': 'MAG', 'gene': 'Gene',
                  'gene_module': 'Gene Module', 'metabolomics': 'Metabolomics',
                  'joint': 'Joint-RPCA'}

def ordination_scatterplots(rpca_results_plt, x="PC1", y="PC2",
                            hue=None, hue_order=None, palette="tab10", 
                            markers=None, style=None, style_order=None, 
                            colorbar=False, cbar_label=None, point_size=100, 
                            title_ = omic_title, xlim=None, 
                            subplots=(2,3), figsize=(25, 10),
                            save_fig=False, save_path=None):

    if subplots is not None:
        fig, axn = plt.subplots(subplots[0], subplots[1], figsize=figsize)
        axn = axn.flatten()
    else:
        fig = plt.figure(figsize=figsize)
        axn = [plt.gca()]

    

    for ax, (tblid, ord_plt) in zip(axn, rpca_results_plt.items()):
        #plotting
        sns.scatterplot(x=x, y=y, hue=hue, hue_order=hue_order, 
                        palette=palette, style=style, 
                        style_order=style_order, markers=markers, 
                        data=ord_plt, s=point_size, ax=ax)
        ax.set_xlabel(x, color='black', weight='bold', fontsize=18, fontname='Arial')
        ax.set_ylabel(y, color='black', weight='bold', fontsize=18, fontname='Arial')
        if title_ is None:
            ax.set_title(tblid, color='black', weight='bold', 
                         fontsize=20, fontname='Arial')
        else:
            ax.set_title(title_[tblid], color='black', weight='bold', 
                        fontsize=20, fontname='Arial')

        if xlim is not None:
            ax.set_xlim(xlim)

        #fix backround
        ax.set_facecolor('white')
        ax.set_axisbelow(True)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.spines['bottom'].set_visible(True)
        ax.spines['top'].set_visible(False)
        for child in ax.get_children():
            if isinstance(child, matplotlib.spines.Spine):
                child.set_color('grey')

        for tick in ax.get_yticklabels():
            tick.set_fontproperties('arial')
            tick.set_color("black")
            tick.set_fontsize(15)
    
        for tick in ax.get_xticklabels():
            tick.set_fontproperties('arial')
            tick.set_color("black")
            tick.set_fontsize(15)
        
        hue_dtype = ord_plt[hue].dtype
        if (hue_dtype!='float64') and (hue_dtype!='int64'):
            ax.legend_.remove()
        
        if (subplots is not None) and (colorbar):
            norm = plt.Normalize(ord_plt[hue].min(), ord_plt[hue].max())
            sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
            sm.set_array([])
            ax.legend_.remove()
            cbar = ax.figure.colorbar(sm, location='right', shrink=0.8,
                                      anchor=(0.0, 0.5), ax=ax)
            cbar.ax.tick_params(labelsize=14)
            cbar.set_label(cbar_label, labelpad=-70, y=0.5, 
                           fontname='Arial', fontsize=15, weight='bold')
            
    plt.tight_layout()

    if (subplots is not None) and (not colorbar):
        ncats = ord_plt[hue].nunique()
        legend = ax.legend(loc=2, bbox_to_anchor=(-1.1, 2.6),
                           prop={'size':15, 'family':'Arial'}, 
                           title="", fancybox=True, framealpha=.0,
                           ncol=ncats, markerscale=1.5)
        legend.get_title().set_fontsize('15')
        #increase the line width in the legend 
        for line in legend.get_lines()[:]:
            line.set_linewidth(2.0)
        for line in legend.get_lines()[:]:
            line.set_linewidth(2.0)
    
    elif (subplots is None) and (not colorbar):
        legend = ax.legend(loc=2, bbox_to_anchor=(1, 1),
                           prop={'size':15, 'family':'Arial'}, 
                           title="", fancybox=True, framealpha=.0,
                           ncol=1, markerscale=1.5)
        legend.get_title().set_fontsize('15')
        #increase the line width in the legend 
        for line in legend.get_lines()[:]:
            line.set_linewidth(2.0)
        for line in legend.get_lines()[:]:
            line.set_linewidth(2.0)

    if (colorbar) and (subplots is None):
        norm = plt.Normalize(ord_plt[hue].min(), ord_plt[hue].max())
        sm = plt.cm.ScalarMappable(cmap=palette, norm=norm)
        sm.set_array([])

        # Remove the legend and add a colorbar
        ax.get_legend().remove()
        cbar = ax.figure.colorbar(sm, location='right', shrink=0.8,
                                  anchor=(0.0, 0.5), ax=axn[-1])
        cbar.ax.tick_params(labelsize=14)
        cbar.set_label(cbar_label, labelpad=-65, y=0.5, 
                       fontname='Arial', fontsize=15, weight='bold')

    if save_fig:
        plt.savefig(save_path, 
                    dpi=600, 
                    bbox_inches='tight',
                    facecolor=fig.get_facecolor(), 
                    edgecolor='none')
    plt.show()