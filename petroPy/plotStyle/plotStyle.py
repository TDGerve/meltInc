import matplotlib.pyplot as plt
    
colors = {
    'flatDesign': plt.cycler(color= ['#e27a3d', '#344d5c', '#df5a49', '#43b29d', '#efc94d']), 
    'firenze': plt.cycler(color= ['#8E2800','#468966', '#B64926','#FFF0A5', '#FFB03B']),
    'vitaminC': plt.cycler(color= ['#FD7400','#004358', '#FFE11A', '#1F8A70', '#BEDB39']),
    'bella': plt.cycler(color= ['#801637', '#047878', '#FFB733', '#F57336', '#C22121']),
    'buddha': plt.cycler(color= ['#192B33', '#FF8000', '#8FB359', '#FFD933', '#CCCC52']),
    'elemental': plt.cycler(color= ['#E64661', '#FFA644', '#998A2F', '#2C594F', '#002D40']),
    'carolina': plt.cycler(color= ['#73839C', '#2E4569', '#AECCCF', '#D5957D', '#9C7873']),
    '42': plt.cycler(color= ['#2469A6', '#C4E1F2', '#F2E205', '#F2D22E', '#D9653B']),
    'terrazaverde': plt.cycler(color= ['#DFE2F2', '#88ABF2', '#4384D9', '#56BFAC', '#D9B341'])
    }


def layout(fontSize= 16, axTitleSize= 16, axLabelSize= 16, tickLabelSize= 12, legendFontSize= 8, colors= colors['firenze']):

    plt.rcParams['figure.constrained_layout.use'] = True
    plt.rcParams['figure.figsize'] = (8, 7)
    plt.rcParams['savefig.dpi'] = 300

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.size'] = fontSize
    plt.rcParams['axes.titlesize'] = axTitleSize
    plt.rcParams['axes.labelsize'] = axLabelSize
    plt.rcParams['xtick.labelsize'] = tickLabelSize
    plt.rcParams['ytick.labelsize'] = tickLabelSize    
    plt.rcParams['legend.fontsize'] = legendFontSize

    plt.rcParams['axes.linewidth'] = 1.5
    plt.rcParams['lines.linewidth'] = 4

    plt.rcParams['lines.markeredgecolor'] = 'k'
    plt.rcParams['lines.markersize'] = 10

    plt.rcParams['axes.prop_cycle'] = colors
    plt.rcParams['axes.facecolor'] = 'whitesmoke' 

    plt.rcParams['axes.grid'] = True
    plt.rcParams['grid.color'] = 'snow'
    plt.rcParams['axes.axisbelow'] = True


   
