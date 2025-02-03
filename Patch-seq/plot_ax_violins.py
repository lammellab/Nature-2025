import os
import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib import pyplot as plt
from sklearn.neighbors.kde import KernelDensity

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 0.5

def get_sequencing_data(celltypes):
    # for given set of celltypes, get sequencing data for doing violin plot
    
    # read in all datalabels
    file_loc = '/Users/davidlukacsovich/Documents/Jupyter/Reference Data/sequencing data'
    dataframes = []
    for datalabel in {cell.split('-')[0] for cell in celltypes}:
        print(datalabel)
        path = os.path.join(file_loc, '%s.txt' % datalabel)
        df = pd.read_csv(path, sep='\t')
        print(df.columns)
        df['CellType'] = datalabel + '-' + df['CellType']
        dataframes.append(df)

    # merge data
    df = pd.concat(dataframes)
    df = df.loc[df['CellType'].isin(celltypes)]
    
    return df

def get_electrophys_data(categories):
    # for a given set of cell categories, get electrophys data for doing violint plot
    print('function not yet implemented')
    assert 0
    
    # read in all datalabels
    file_loc = '/home/foldy_lab/Documents/Newest/Soma Paper/Produce Figures/created references/electrophys'
    dataframes = [pd.read_csv(os.path.join(file_loc, '%s.txt' % category), sep='\t') for category in categories]
    
    return dataframes

def define_axes(fig, top, num, left=.1, right=.9, row_count=6, height=.06, dw=.045, dh=.02):
    # define n axes to fill up rows of row_count each
    dw = dw
    width = (right-left) / row_count - dw
    dh = dh
    bot = top - height
    axes = [fig.add_axes([left+dw+(dw+width)*(i%row_count), bot-(dh+height)*(i//row_count), width, height]) for i in range(num)]
    return axes, bot-(dh+height)*((num-1)//row_count)

def plot_ax_violin(ax, data, xval, color, width=.25, show_violin=True, show_error=False, show_scatter=True):
    # do a violin plot of points around xval position
    
    # define points to run KDE on and KDE step size
    data = data[np.isfinite(data)]
    if data.size==0:
        return
    low, high = np.min(data), np.max(data)
    if high == 0:
        high = 1
    if low == high:
        low = high*.9
    diff = high-low
    step = diff/5
    yvals = np.linspace(low,high,100)
    
    # generate KDE
    kde = KernelDensity(kernel='gaussian', bandwidth=step).fit(data[:,np.newaxis])
    log_dens = np.exp(kde.score_samples(yvals[:,np.newaxis]))
    log_dens = log_dens / np.max(log_dens)*width
    
    # generate randomly scattered points
    # should be randomly distributed around xval (so they don't overlap), where variation scales with kde
    xvals = xval + (np.random.rand(data.size)-.5)*log_dens[np.searchsorted(yvals,data)]
    
    if show_violin:
        # plot KDE
        ax.fill_betweenx(yvals, log_dens+xval, xval-log_dens, alpha=.5, color=color, linewidth=0, edgecolor='face')
    if show_error:
        mean, error = np.mean(data), np.std(data) / np.sqrt(data.size)
        ax.errorbar(xval, mean, yerr=error, color=color, capsize=4)
        ax.errorbar(xval, mean, yerr=0, color=color, capsize=6)
    
    # plot points
    if show_scatter:
        ax.scatter(xvals, data, s=1, color=color)
    
    return

def plot_dataframe(fig, top, df_violin, df_colors, left=.1, right=.9, row_count=6, rotation=0, ticklabels=False, limits=[], show_violin=True, show_error=False, dh=.02):
    # initialize variables
    if len(limits) < df_violin.shape[1]:
        limits += [(None,None)] * (df_violin.shape[1] - len(limits))
        
    if not ticklabels:
        ticklabels = df_colors.index.tolist()
    
    positions = np.arange(df_colors.shape[0]) + .5
    
    # set up plots
    axes, bot = define_axes(fig, top, df_violin.shape[1], left=left, right=right, row_count=row_count, dh=dh)
    
    # do plots
    for ax, label, limit in zip(axes, df_violin.columns, limits):
        for position, celltype, color in zip(positions, df_colors.index, df_colors.Color):
            data = df_violin.loc[celltype, label]
            plot_ax_violin(ax, data, position, color, width=.4, show_violin=show_violin, show_error=show_error)
        ax.set_ylabel(label, fontsize=8, labelpad=1, family='Arial')
        ax.set_xlim(0, df_colors.shape[0])
        ax.set_ylim(limit[0], limit[1])
        ax.tick_params(size=1, labelsize=7, pad=1)
        ax.set_xticks(positions)
        ax.set_xticklabels(ticklabels, rotation=rotation, family='Arial')
    
    return bot
    

def plot_sequencing(fig, top, celltypes, celltype_labels, colors=['blue', 'red'], tissue=False, left=.1, right=.9, use_cols='All', row_count=6, rotation=0, ticklabels=False, limits=None, show_violin=True, show_error=False):
    # do plotting of sequencing on figure
    
    if len(colors) < len(celltypes):
        colors += [np.random.rand(3) for i in range(len(celltypes)-len(colors))]
    
    # read in data
    df = get_sequencing_data(celltypes)
    
    # create sequencing data
    label_order = ['Reads', 'Aligned', 'Mapped', 'GeneCount', 'AlignRate', 'MapRate']
    label_names = ['Reads (M)', 'Aligned (M)', 'Mapped (M)', 'Gene Count (K)', 'Align Rate (%)', 'Mapping Rate (%)']
    normalization = [1/1000000, 1/1000000, 1/1000000, 1/1000, 1., 1.]
    if not limits:
        if not tissue:
            limits = [(0,25), (0, 25), (0, 25), (0, 12), (0, 100), (0, 100)]
        else:
            limits = [(0,40), (0, 30), (0, 20), (0, 18), (0, 100), (0, 100)]
    
    # normalize data
    for label, norm in zip(label_order, normalization):
        df[label] = df[label]*norm
        
    # pick only categories we want to plot
    if use_cols == 'All':
        use_cols = [i for i in range(len(limits))]
    plot_count = len(use_cols)
    
    label_order = [label_order[col] for col in use_cols]
    label_names = [label_names[col] for col in use_cols]
    normalization = [normalization[col] for col in use_cols]
    limits = [limits[col] for col in use_cols]
    kept = ['CellType'] + label_order
    df = df[kept].set_index('CellType')
    
    positions = np.arange(len(celltypes))+.5
    
    # do plotting
    axes, bot = define_axes(fig, top, plot_count, left=left, right=right, row_count=row_count)
    #ticklabels = [category.split('-')[0] for category in categories]
    if not ticklabels:
        ticklabels = [celltype_label.replace('-','-\n') for celltype_label in celltype_labels]
    for ax, label_name, label_ord, limit in zip(axes, label_names, label_order, limits):
        for celltype, position, color in zip(celltypes, positions, colors):
            data = df.loc[celltype, label_ord]
            plot_ax_violin(ax, data, position, color, width=.4, show_violin=show_violin, show_error=show_error)
        ax.set_ylabel(label_name, fontsize=8, labelpad=1, **hfont)
        ax.axis([0,len(celltypes),limit[0],limit[1]])
        ax.tick_params(size=1, labelsize=7, pad=1)
        ax.set_xticks(positions)
        ax.set_xticklabels(ticklabels, rotation=rotation, **hfont)

    return bot

def plot_generated_electrophys(fig, top, df, categories, colors=['red', 'blue'], left=.1, right=.9, use_cols='All', row_count=6, rotation=0, ticklabels=False, limits=None, show_violin=True, show_error=False):
    # do plotting of electrophys data in df_elec on figure
    
    # create sequencing data
    label_order = ['Voltage Resting (mV)', 'Max. Frequency (Hz)', 'Firing Threshold (pA)', 'Res. Input (MΩ)', 'Res. Series (MΩ)', 'Capacitance (pF)', 'Basewidth (ms)', 'Halfwidth (ms)', 'Symmetricity', 'Peak (mV)', 'Saq Potential (mV)', 'Attenuation']
    label_names = ['V resting (mV)', 'Max. AP firing (Hz)', 'Firing thresh (pA)', 'Input res. (MΩ)', 'Series res. (MΩ)', 'Capacitance (pF)', 'AP base-width (ms)', 'AP half-width (ms)', 'AP symmetricity', 'AP peak (mV)', 'Sag Potential (mV)', 'Attenuation']
    label_names = label_order
    
    if not limits:
        limits = [(-85,-40), (0, 125), (0, 180), (80, 500), (5, 30), (20, 120),
                 (.4,1.4), (.2,.8), (.28,.46), (60, 105), (0,15), (1.,1.7)]
    
    # pick only categories we want to plot
    if use_cols == 'All':
        use_cols = [i for i in range(len(limits))]
    plot_count = len(use_cols)
    
    label_order = [label_order[col] for col in use_cols]
    label_names = [label_names[col] for col in use_cols]
    limits = [limits[col] for col in use_cols]
    kept = ['CellTypes'] + label_order
    df = df[kept].set_index('CellTypes')
    
    positions = np.arange(len(categories))+.5
    
    # do plotting
    axes, bot = define_axes(fig, top, plot_count, left=left, right=right, row_count=row_count)
    #ticklabels = [category.split('-')[0] for category in categories]
    if not ticklabels:
        ticklabels = [category.replace('-','-\n') for category in categories]
    for ax, label_name, label_ord, limit in zip(axes, label_names, label_order, limits):
        for category, position, color in zip(categories, positions, colors):
            data = df.loc[category, label_ord]
            plot_ax_violin(ax, data, position, color, width=.4, show_violin=show_violin, show_error=show_error)
        ax.set_ylabel(label_name, fontsize=8, labelpad=1, **hfont)
        ax.axis([0,len(categories),limit[0],limit[1]])
        ax.tick_params(size=1, labelsize=7, pad=1)
        ax.set_xticks(positions)
        ax.set_xticklabels(ticklabels, rotation=rotation, **hfont)

    return bot

def plot_electrophys(fig, top, categories, refs, colors=['red', 'blue'], left=.1, right=.9, use_cols='All', row_count=6, height=.06, show_violin=True, show_error=False):
    # do plotting of electrophys on figure

    # read in data
    celltypes = {cell for category in categories for cell in refs[category]}
    dataframes = get_electrophys_data(categories)
    
    # create sequencing data
    label_order = ['Voltage Resting', 'Maximum Frequency', 'Firing Threshold', 'Resistance Input', 'Resistance Series', 'Capacitance', 'Basewidth', 'Halfwidth', 'Symmetricity', 'Peak', 'Saq Potential', 'Attenuation']
    label_names = ['V resting (mV)', 'Max. AP firing (Hz)', 'Firing thresh (pA)', 'Input res. (MΩ)', 'Series res. (MΩ)', 'Capacitance (pF)', 'AP base-width (ms)', 'AP half-width (ms)', 'AP symmetricity', 'AP peak (mV)', 'Sag Potential (mV)', 'Attenuation']
    
    limits = [(-85,-40), (0, 125), (0, 180), (80, 500), (5, 30), (20, 120),
             (.4,1.4), (.2,.8), (.28,.46), (60, 105), (0,15), (1.,1.7)]
    
    # pick only categories we want to plot
    if use_cols == 'All':
        use_cols = [i for i in range(len(limits))]
    plot_count = len(use_cols)
    
    label_order = [label_order[col] for col in use_cols]
    label_names = [label_names[col] for col in use_cols]
    limits = [limits[col] for col in use_cols]
    
    plot_data = [[dataframe[label] for dataframe in dataframes] for label in label_order]
    positions = np.arange(len(categories))+.5
    
    # do plotting
    axes, bot = define_axes(fig, top, plot_count, left=left, right=right, row_count=row_count, height=height, dw=0.045)
    ticklabels = [category.replace('-','-\n') for category in categories]
    for ax, label, data, limit in zip(axes, label_names, plot_data, limits):
        #ax.violinplot(data, positions=positions, showextrema=False)
        for cat_data, position, color in zip(data, positions, colors):
            plot_ax_violin(ax, cat_data, position, color, width=.25, show_violin=show_violin, show_error=show_error)
        ax.set_ylabel(label, fontsize=8, labelpad=1, **hfont)
        ax.axis([0,len(categories),limit[0],limit[1]])
        ax.tick_params(size=1, labelsize=7, pad=1)
        ax.set_xticks(positions)
        ax.set_xticklabels(ticklabels, **hfont)
    
    return bot