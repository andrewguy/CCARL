'''Plotting utility functions'''


def remove_top_right_borders(ax):
    '''Remove top and right borders from Matplotlib axis'''
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
