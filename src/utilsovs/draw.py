#draw.py
from utilsovs.globals import COLORSCHEME
import matplotlib.pyplot as plt
import pandas as pd
import logomaker

def draw_seqLogo(data,showplot,center_values):

    fig, ax = plt.subplots()

    ax.set_aspect(1.0/ax.get_data_ratio()*0.4)

    print (data.df)

    data.seqLogo = logomaker.Logo(data.df,
                          ax=ax,
                          shade_below=.5,
                          fade_below=.5,
                          center_values=center_values,
                          flip_below=False,
                          stack_order='small_on_top',
                          font_name='Arial Rounded MT Bold')

    for aa in COLORSCHEME.keys():
        data.seqLogo.style_single_glyph(p=0,c=aa,floor=0,ceiling=0)

    data.seqLogo.style_spines(visible=False)
    data.seqLogo.style_spines(spines=['left', 'bottom'], visible=True)
    data.seqLogo.style_xticks(rotation=90, fmt='%d', anchor=0)
    data.seqLogo.ax.set_ylabel(r'$\log_2\frac{observed}{abundance}$', labelpad=5,  fontsize=14)
    data.seqLogo.ax.set_xlabel("Position", labelpad=5,  fontsize=14)
    data.seqLogo.ax.xaxis.set_ticks_position('none')
    data.seqLogo.ax.yaxis.set_ticks_position('none')
    data.seqLogo.ax.xaxis.set_tick_params(pad=1)

    data.seqLogo.fig.tight_layout()

    if data.filepath != None:
        plt.savefig(data.filepath,orientation='landscape', bbox_inches = 'tight',transparent=True, dpi=300)
        print ('File %s successfully written to disk by %s routine' % (data.filepath,data.func))

    if showplot == True:
        plt.show()

    plt.close()
    return data
