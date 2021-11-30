import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from ccarl.glycan_graph_methods import generate_digraph_from_glycan_string


from ccarl.glycan_plotting import draw_glycan_diagram


def render_features_pdf(model, pdfname):
    '''Renders features in a PDF document, with one feature per page.

    Args:
        model (ccarl.CCARLClassifier): Trained model
        pdfname (str): PDF output filename.
    Returns:
        None
    '''
    coefs = model._model.coef_[0]
    _sorted_coefs, sorted_features = zip(*sorted(zip(coefs, model.subtree_features), key=lambda x: -x[0]))
    pp = PdfPages(pdfname)
    for i, feature in enumerate(sorted_features):
        fig, ax = plt.subplots()
        fig.suptitle(f"Feature {i+1}")
        draw_glycan_diagram(feature, ax, draw_terminal_connection_labels=True)
        pp.savefig(fig)
        plt.close(fig)
    pp.close()
    return


def render_glycan_list_pdf(glycans, pdfname, format="CFG"):
    '''Renders glcyans in a PDF document, with one glycan per page.

    Args:
        glycans (list): List of glycan strings.
        pdfname (str): PDF output filename.
    Returns:
        None
    '''
    pp = PdfPages(pdfname)
    glycan_graphs = [generate_digraph_from_glycan_string(x, parse_linker=True, format=format) for x in glycans]
    for glycan, glycan_str in zip(glycan_graphs, glycans):
        fig, ax = plt.subplots()
        fig.suptitle(glycan_str)
        draw_glycan_diagram(glycan, ax, draw_terminal_connection_labels=True)
        pp.savefig(fig)
        plt.close(fig)
    pp.close()
    return
