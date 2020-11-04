import matplotlib
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from scipy import ndimage
import os as os
import math

from matplotlib.colors import LinearSegmentedColormap

class PlottingPackage(object):
    def __init__(self, cwd=None, save_dir_name="test"):
        """ A class for high-level customizaiton of figures
           
        The intention is to initialize with a save_dir_name as well as
        the super diretory where you want the save_dir_name to be.
        
        Then you select the dimensions you want such as set_X_dimensions()
        
        Standard considerations include a left axis with a label, a bottom
        axis with a label, a title at the top and no axis labels on the right.
        
        Standard font sizes are also given that were found to be reasonably
        legible in their respective formats. Feel free to adjust and modify
        as necessary.   
        
        """
        if cwd is None:
            self.cwd = os.getcwd()
        else:
            self.cwd = cwd

        self.set_save_dir(save_dir_name)

        self.colors = ["b", "r", "g", "m", "c", "y"]
        self.colors_fep = ["b", "r", "g", "m", "c"]
        self.colors_black = ["k"]

        self.line_types = ["-", "--", ":", "-."]

        self.alphabet = [chr(i) for i in range(ord('a'),ord('z')+1)]

        for col in self.colors:
            self.colors_black.append(col)

    def set_save_dir(self, save_dir_name):
        self.figures_dir = "%s/%s" % (self.cwd, save_dir_name) # directory for saving final figures
        self.png_figures_dir = "%s/%s_png" % (self.cwd, save_dir_name)

        ensure_dir(self.figures_dir)
        ensure_dir(self.png_figures_dir)

    def save_figure(self, fig, savefile):
        fig.savefig("%s/%s.pdf" % (self.figures_dir, savefile),  format="pdf", bbox_inches="tight", pad_inches=0.02, dpi=300)

        png_save_name = "%s/%s" % (self.png_figures_dir, savefile)

        fig.savefig("%s.png" % png_save_name, bbox_inches="tight", format="png", dpi=300)
        #fig.savefig("%s.png" % png_save_name, format="png", dpi=300)

    def set_presentation_dimensions(self):
        self.column_width_inches = 5 #width of column in a two-column article
        self.double_column_width_inches = 9. # width of a two-column spanning figure
        self.standard_box_height_inches = self.column_width_inches * (5. / 6.)
        self.padding_standard = 0.8
        self.padding_standard_width_inches = 0.75 # padding on the left for y-axis label
        self.padding_standard_height_inches = 0.75 # padding on the bottom for x-axis label
        self.padding_title_height_inches = 0.8 # padding on top for title
        self.padding_empty_inches = 0.16 # padding if no labels or tick marks exist
        self.padding_buffer_inches = 0.1 # extra padding between adjacent figures
        self.padding_buffer_xaxis = 0.3
        self.padding_buffer_yaxis = 0.3

        self.standard_font_size = 28
        self.standard_special_font_size = 38
        self.standard_thinline = 2
        self.standard_thickline = 4


        matplotlib.rcParams.update({"font.size":self.standard_font_size})
        matplotlib.rcParams.update({"axes.labelsize":self.standard_font_size})

        matplotlib.rcParams.update({"axes.labelpad":4.0})
        matplotlib.rcParams.update({"legend.fontsize":self.standard_font_size})

        matplotlib.rcParams.update({"xtick.labelsize":2*self.standard_font_size/3})
        matplotlib.rcParams.update({"xtick.major.size":3.5})
        matplotlib.rcParams.update({"xtick.minor.size":2.})
        matplotlib.rcParams.update({"xtick.major.pad":3.5})

        matplotlib.rcParams.update({"ytick.labelsize":2*self.standard_font_size/3})
        matplotlib.rcParams.update({"ytick.major.size":3.5})
        matplotlib.rcParams.update({"ytick.minor.size":2.})
        matplotlib.rcParams.update({"ytick.major.pad":3.5})

        matplotlib.rcParams.update({"lines.linewidth":2})
        matplotlib.rcParams.update({"patch.linewidth":1})
        matplotlib.rcParams.update({"axes.linewidth":1})

    def set_acs_dimensions(self):
        self.column_width_inches = 3.25 #width of column in a two-column article
        self.double_column_width_inches = 7. # width of a two-column spanning figure
        self.standard_box_height_inches = self.column_width_inches * (5. / 6.)
        self.padding_standard = 0.4
        self.padding_standard_width_inches = 0.32 # padding on the left for y-axis label
        self.padding_standard_width_label_inches = 0.2 # padding on the left for the y-axis when it is labeled but no tick-labels
        self.padding_standard_height_inches = 0.25 # padding on the bottom for x-axis label
        self.padding_title_height_inches = 0.20 # padding on top for title
        self.padding_empty_inches = 0.08 # padding if no labels or tick marks exist
        self.padding_buffer_inches = 0.1 # extra padding between adjacent figures
        self.padding_buffer_xaxis = 0.15 # extra padding along vertical dimension
        self.padding_buffer_yaxis = 0.15 # extra padding along horizontal dimension

        self.standard_font_size = 9
        self.standard_special_font_size = 13
        self.standard_thinline = 1
        self.standard_thickline = 2

        matplotlib.rcParams.update({"font.size":9})
        matplotlib.rcParams.update({"axes.labelsize":9})

        matplotlib.rcParams.update({"axes.labelpad":2.0})
        matplotlib.rcParams.update({"legend.fontsize":8})

        matplotlib.rcParams.update({"xtick.labelsize":6})
        matplotlib.rcParams.update({"xtick.major.size":1.75})
        matplotlib.rcParams.update({"xtick.minor.size":1.})
        matplotlib.rcParams.update({"xtick.major.pad":1.75})

        matplotlib.rcParams.update({"ytick.labelsize":6})
        matplotlib.rcParams.update({"ytick.major.size":1.75})
        matplotlib.rcParams.update({"ytick.minor.size":1.})
        matplotlib.rcParams.update({"ytick.major.pad":1.75})

        matplotlib.rcParams.update({"lines.linewidth":1})
        matplotlib.rcParams.update({"patch.linewidth":0.5})
        matplotlib.rcParams.update({"axes.linewidth":0.5})

    def set_pnas_dimensions(self):
        self.column_width_inches = 3.42 #width of column in a two-column article
        self.double_column_width_inches = 7. # width of a two-column spanning figure
        self.standard_box_height_inches = self.column_width_inches * (5. / 6.)
        self.padding_standard = 0.4
        self.padding_standard_width_inches = 0.32 # padding on the left for y-axis label
        self.padding_standard_width_label_inches = 0.2 # padding on the left for the y-axis when it is labeled but no tick-labels
        self.padding_standard_height_inches = 0.25 # padding on the bottom for x-axis label
        self.padding_title_height_inches = 0.20 # padding on top for title
        self.padding_empty_inches = 0.08 # padding if no labels or tick marks exist
        self.padding_buffer_inches = 0.1 # extra padding between adjacent figures
        self.padding_buffer_xaxis = 0.15 # extra padding along vertical dimension
        self.padding_buffer_yaxis = 0.15 # extra padding along horizontal dimension

        self.standard_font_size = 9
        self.standard_special_font_size = 13
        self.standard_thinline = 1
        self.standard_thickline = 2

        matplotlib.rcParams.update({"font.size":9})
        matplotlib.rcParams.update({"axes.labelsize":9})

        matplotlib.rcParams.update({"axes.labelpad":2.0})
        matplotlib.rcParams.update({"legend.fontsize":8})

        matplotlib.rcParams.update({"xtick.labelsize":6})
        matplotlib.rcParams.update({"xtick.major.size":1.75})
        matplotlib.rcParams.update({"xtick.minor.size":1.})
        matplotlib.rcParams.update({"xtick.major.pad":1.75})

        matplotlib.rcParams.update({"ytick.labelsize":6})
        matplotlib.rcParams.update({"ytick.major.size":1.75})
        matplotlib.rcParams.update({"ytick.minor.size":1.})
        matplotlib.rcParams.update({"ytick.major.pad":1.75})

        matplotlib.rcParams.update({"lines.linewidth":1})
        matplotlib.rcParams.update({"patch.linewidth":0.5})
        matplotlib.rcParams.update({"axes.linewidth":0.5})

    def set_thesis_dimensions(self):
        self.column_width_inches = 4 #width of column in a two-column article
        self.double_column_width_inches = 6. # width of a two-column spanning figure
        self.standard_box_height_inches = self.column_width_inches * (5. / 6.)
        self.padding_standard = 0.5
        self.padding_standard_width_inches = 0.42 # padding on the left for y-axis label
        self.padding_standard_width_label_inches = 0.3 # padding on the left for the y-axis when it is labeled but no tick-labels
        self.padding_standard_height_inches = 0.35 # padding on the bottom for x-axis label
        self.padding_title_height_inches = 0.30 # padding on top for title
        self.padding_empty_inches = 0.08 # padding if no labels or tick marks exist
        self.padding_buffer_inches = 0.1 # extra padding between adjacent figures
        self.padding_buffer_xaxis = 0.15 # extra padding along vertical dimension
        self.padding_buffer_yaxis = 0.15 # extra padding along horizontal dimension

        self.standard_font_size = 12
        self.standard_special_font_size = 16
        self.standard_thinline = 1
        self.standard_thickline = 2

        matplotlib.rcParams.update({"font.size":self.standard_font_size})
        matplotlib.rcParams.update({"axes.labelsize":self.standard_font_size})

        matplotlib.rcParams.update({"axes.labelpad":2.0})
        matplotlib.rcParams.update({"legend.fontsize":self.standard_font_size})

        matplotlib.rcParams.update({"xtick.labelsize":self.standard_font_size})
        matplotlib.rcParams.update({"xtick.major.size":1.75})
        matplotlib.rcParams.update({"xtick.minor.size":1.})
        matplotlib.rcParams.update({"xtick.major.pad":1.75})

        matplotlib.rcParams.update({"ytick.labelsize":self.standard_font_size})
        matplotlib.rcParams.update({"ytick.major.size":1.75})
        matplotlib.rcParams.update({"ytick.minor.size":1.})
        matplotlib.rcParams.update({"ytick.major.pad":1.75})

        matplotlib.rcParams.update({"lines.linewidth":1})
        matplotlib.rcParams.update({"patch.linewidth":0.5})
        matplotlib.rcParams.update({"axes.linewidth":0.5})

    def get_standard_axes(self, fig, fig_dimensions):
        """ Return an axes class that maximizes the use of a figure area per standard considerations"""
        fig_width_inches = fig_dimensions[0]
        fig_height_inches = fig_dimensions[1]
        left_padding_inches = self.padding_standard_width_inches
        right_padding_inches = self.padding_empty_inches
        top_padding_inches = self.padding_empty_inches
        lower_padding_inches = self.padding_standard_height_inches

        axes_width_inches = fig_width_inches - left_padding_inches - right_padding_inches
        axes_height_inches = fig_height_inches - top_padding_inches - lower_padding_inches

        left = left_padding_inches / fig_width_inches
        bottom = lower_padding_inches / fig_height_inches

        axes_width = axes_width_inches / fig_width_inches
        axes_height = axes_height_inches / fig_height_inches

        ax_plot = fig.add_axes([left, bottom, axes_width, axes_height])

        return ax_plot
        
    def plot_lines(self, x_list, y_list, savename, xlabel="x", ylabel="y", labels=None, disable_yticks=False, disable_xticks=False, colors=None, alphas=None,  double=False, width=None, height=None, legends_location="upper right", y_padding_factor=1.05, bigger_font=False):
        """ Plot multiple lines on a standard axes size."""
        if width is None:
            if double:
                fig_width_inches = self.double_column_width_inches
            else:
                fig_width_inches = self.column_width_inches
        else:
            fig_width_inches = width
        if height is None:
            fig_height_inches = 0.5 * fig_width_inches
        else:
            fig_height_inches = height

        if colors == None:
            use_colors = self.colors
        else:
            use_colors = colors

        if alphas == None:
            alphas = [1 for i in range(len(x_list))]


        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        ax = self.get_standard_axes(fig, (fig_width_inches, fig_height_inches))

        all_lines = []
        for idx, (x,y) in enumerate(zip(x_list, y_list)):
            line, = ax.plot(x, y, color=use_colors[idx], alpha=alphas[idx])
            all_lines.append(line)

        if labels is not None:
            ax.legend(all_lines, labels, loc=legends_location, fontsize=self.standard_font_size)

        if disable_xticks:
            ax.get_xaxis().set_ticks([])
        if disable_yticks:
            ax.get_yaxis().set_ticks([])

        if bigger_font:
            ax.set_xlabel(xlabel, fontsize=self.standard_special_font_size)
            ax.set_ylabel(ylabel, fontsize=self.standard_special_font_size)
        else:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        xmin, xmax = determine_min_max(x_list)
        ymin, ymax = determine_min_max(y_list)
        ymax *= y_padding_factor

        ax.axis([xmin, xmax, ymin, ymax])

        self.save_figure(fig, savename)
     
    def plot_spread(self, data, savename, bins=None, xname="data", yname="count"):
        """ Plot a histogram of data. Construct axes dimensions manually"""
        if bins is None:
            bins = np.sqrt(len(data))
            bins = int(bins)
        fig_width_inches = self.column_width_inches
        fig_height_inches = self.column_width_inches * (5./6.)
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))
        left_padding_inches = self.padding_standard_width_inches
        right_padding_inches = self.padding_empty_inches
        top_padding_inches = self.padding_empty_inches
        lower_padding_inches = self.padding_standard_height_inches

        axes_width_inches = fig_width_inches - left_padding_inches - right_padding_inches
        axes_height_inches = fig_height_inches - top_padding_inches - lower_padding_inches

        left = left_padding_inches / fig_width_inches
        bottom = lower_padding_inches / fig_height_inches

        axes_width = axes_width_inches / fig_width_inches
        axes_height = axes_height_inches / fig_height_inches

        ax_plot = fig.add_axes([left, bottom, axes_width, axes_height])

        ax_plot.hist(data, bins=bins)

        ax_plot.set_xlabel(xname)
        ax_plot.set_ylabel(yname)

        self.save_figure(fig, savename)

    def plot_many_images(self, image_list, savename, ncols=1, nrows=1):
        """ Give a list of images of the same sizes and arrange in a grid
           
        This seems highly specialized, but actually really useful if you
        have multiple protein structures you want to show in a single 
        figure with subfigures.
        
        """
        first_im = mpimg.imread("%s" % image_list[0])
        ratio = float(np.shape(first_im)[0]) / float(np.shape(first_im)[1])
        assert len(image_list) == ncols * nrows

        fig_width_inches = self.column_width_inches
        fig_height_padding_inches = self.padding_standard_width_inches

        image_height_inches = ratio * fig_width_inches / ncols

        fig_height_inches = (fig_height_padding_inches + image_height_inches) * nrows
        fig = plt.figure(figsize=(fig_width_inches, fig_height_inches))

        image_height = image_height_inches / fig_height_inches
        vertical_padding = fig_height_padding_inches / fig_height_inches
        vertical_diff = vertical_padding + image_height
        highest_vertical_start = 1.0 - vertical_diff
        for count,image in enumerate(image_list):
            horizontal_start = (1.0 / (ncols)) * (count % 2)
            vertical_start = highest_vertical_start - (vertical_diff * (np.floor(count/ncols)))
            ax1 = fig.add_axes([horizontal_start, vertical_start, 1.0/ncols, image_height])
            ax1.xaxis.set_visible(False)
            ax1.yaxis.set_visible(False)
            ax1.axis("off")
            this_img = mpimg.imread("%s" % image_list[count])
            ax1.imshow(this_img)
            ax1.set_title("(%s)" % self.alphabet[count], fontsize=self.standard_special_font_size)

        self.save_figure(fig, savename)

