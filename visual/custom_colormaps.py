"""
The module contains custom colormaps for matplotlib,
as well as wrappers for fast generation of divergent colormaps
"""

import matplotlib.colors as mcolors
import matplotlib as mpl
import matplotlib.pyplot as plt

def _make_colormap(seq,name):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    print(cdict)

    return mcolors.LinearSegmentedColormap(name, cdict)


def get_rwb_colormap(vmin,vmax):
    """
    Create a  divergent colormap with
    with red (negative values), white (near 0) and blue (positive values)
    colors, linearly interpolated between vmin and vmax.
    The corresponding colormap is registered by name 'rwb'
    """
    c = mcolors.ColorConverter().to_rgb
    frac = 1.-abs(vmin)/vmax
    cred = (1, frac, frac)
    dv=vmax+abs(vmin)
    mid = abs(vmin)/dv
    rwb = _make_colormap([cred, c('white'),  mid, c('white'), c('blue')],'rwb')
    plt.register_cmap(name=rwb.name, cmap=rwb)
    return rwb

def test_colormap(cmap,vmin=-1,vmax=1):
    """
    The function plots a colorbar colored according to the cmap,
    with bounds vmin and vmax

    """
    fig, ax = plt.subplots(figsize=(6, 0.5))
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
    cax=ax, orientation='horizontal')
    return
