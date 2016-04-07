import matplotlib as mpl
import numpy as np
mpl.use('Agg')
import matplotlib_venn as venn
from matplotlib import pyplot as plt

venn3_keys = ['100', '010' ,'001', '110', '011', '101', '111']
venn2_keys = ['10', '01', '11']


def venn3(subsets, title, unit_title, filename, set_labels=None, normalize=1.0, annotation=None):
    plt.figure()
    v = venn.venn3_unweighted(subsets=subsets, set_labels=set_labels)
    c = venn.venn3_circles(subsets=(1,1,1,1,1,1,1), linestyle='solid', linewidth=1.5, normalize_to=normalize)
    for i in range(len(venn3_keys)):
        label_id = venn3_keys[i]
        text = v.get_label_by_id(label_id)
        text.set_position(text.get_position() + np.array([0, 0.02]))
        # TEMPORALLY COUPLED WITH CREATION OF DIAGRAM
        subtitle = unit_title
        if text.get_text() != '1':
            subtitle += 's'
        text.set_text(text.get_text() + '\n' + subtitle)
        text.set_size(text.get_size() - 2)
    if annotation is not None:
        for a in annotation:
            text = v.get_label_by_id(a)
            xy= text.get_position() - np.array([0, 0.085])
            plt.annotate(annotation[a], xy=xy, xytext=xy, ha='center', textcoords='offset points', color='r', weight='bold')

    plt.title(title)
    plt.savefig(filename)
    plt.close()


def venn2(subsets, title, unit_title, filename, set_labels=None, normalize=1.0, annotation=None):
    plt.figure()
    v = venn.venn2(subsets=subsets, set_labels=set_labels)
    for i in range(len(venn2_keys)):
        label_id = venn2_keys[i]
        text = v.get_label_by_id(label_id)
        text.set_position(text.get_position() + np.array([0, 0.02]))
        # TEMPORALLY COUPLED WITH CREATION OF DIAGRAM
        subtitle = unit_title
        if text.get_text() != '1':
            subtitle += 's'
        text.set_text(text.get_text() + '\n' + subtitle)
        text.set_size(text.get_size() - 2)
    if annotation is not None:
        for a in annotation:
            text = v.get_label_by_id(a)
            xy = text.get_position() - np.array([0, 0.085])
            plt.annotate(annotation[a], xy=xy, xytext=xy, ha='center', textcoords='offset points', color='r', weight='bold')

    plt.title(title)
    plt.savefig(filename)
    plt.close()

class SimpleTable:
    def __init__(self, header_tuple):
        """
        initializes table as a list of tuples with a header. add rows as tuples to self.rows. sort them as you please
        :return: SimpleTable
        """
        self.header = header_tuple
        self.col = len(header_tuple)
        self.rows = []

    def add(self, row_tuple):
        if len(row_tuple) != self.col:
            raise ValueError
        self.rows.append(row_tuple)

    def markdown(self):
        result = "\n"
        result += SimpleTable._markdown_row(self.header) + "\n"
        separator = "| "
        for i in range(self.col):
            separator += "--- | "
        result += separator + "\n"
        for r in self.rows:
            result += SimpleTable._markdown_row(r) + "\n"
        return result

    @staticmethod
    def _markdown_row(tuple):
        result = "| "
        for r in tuple:
            result += str(r) + " | "
        return result




