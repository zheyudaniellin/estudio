""" plot_example_2.py
This does the exact same as plot_example_1.py, but this is an example of how to use this way of scripting and to practice using the python debugger tool. I highly recommend people who are new to python to use this tool. The rest of the plotting scripts in this directory are based in this form. If this form is confusing be sure to refer back to plot_example_1.py for comparison. 
All the settings will be under the function main() which is usually near the bottom of the file.

Once you changed the parameters for your data, simply run this from the terminal (make sure you're using python version 3):
python plot_example2.py

If you want to stop at some place to look at the variables interactively, simply go to that place in this script and write:
pdb.set_trace()

Now run the code from the terminal like before:
python plot_example2.py

And it will stop at the place. Simply type the variables you want to look at and it will show the values. If you want to let the script continue all the way to the end, then simply enter "c" 

See this webpage to learn more aboue pdb: 
https://docs.python.org/3/library/pdb.html

Here's another tutorial site:
https://realpython.com/python-debugging-pdb/
"""
import numpy as np
import matplotlib.pyplot as plt
import etools as eto
import pdb

def plot1(im):
    ax = plt.gca()

    pc = ax.pcolormesh(im.grid.x, im.grid.y, im.I[:,:,0].T, cmap=plt.cm.gist_heat)
    cb = plt.colorbar(pc, ax=ax)
    ax.set_aspect('equal')

    plt.show()

def main():
    # ==== settings ====
    # distance to source in pc
    dpc = 140.

    # determine the fits name
    fname = '/lustre/cv/projects/edisk/IRAS04302/export_SBLB/IRAS04302_SBLB_continuum_robust_0.0.pbcor.tt0.fits'

    # prepare the imageManager
    im = eto.image.imageManager()

    # read the fits image
    im.read_fits(fname, dpc, quant='I')

    # start plotting 
    plot1(im)

if __name__ == '__main__':
    main()

