import argparse as ap
import numpy as np
from time import sleep

parser = ap.ArgumentParser(description='Process triangulation and solution filenames')
parser.add_argument('-s', type=str,
                    help='name of the solution')
parser.add_argument('-t', type=str,
                    help='name of the triangulation')
args = parser.parse_args()

filename_c4n = "meshes/" + args.t + "_c4n.tsv";
c4n = np.loadtxt(filename_c4n, delimiter="\t")

filename_n4e = "meshes/" + args.t + "_n4e.tsv";
n4e = np.loadtxt(filename_n4e, delimiter="\t").astype(int)

filename_sol = "solutions/" + args.s + ".tsv";
u = np.loadtxt(filename_sol, delimiter="\t")

if len(u.shape) == 1:
    from matplotlib import pyplot as plt
    from matplotlib import tri
    
    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    t = tri.Triangulation(x=c4n[:,0], y=c4n[:,1], triangles=n4e, mask=None)

    tcf = ax.tricontourf(t, u)
    fig.colorbar(tcf)
    plt.show()

elif len(u.shape) == 2:
    from matplotlib import pyplot as plt
    from matplotlib import tri
    import matplotlib.animation as animation

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    t = tri.Triangulation(x=c4n[:,0], y=c4n[:,1], triangles=n4e, mask=None)

    lvls = np.arange(round(np.min(u), 1) - 0.1, round(np.max(u), 1) + 0.1, 0.05 * (np.max(u) - np.min(u)))
    tcf = ax.tricontourf(t, u[:,0], vmin=np.min(u), vmax=np.max(u), levels=lvls)

    fig.colorbar(tcf, ticks=lvls)

    

    def data_gen(framenumber, u, tcf):
        if framenumber == 1:
            sleep(1)
        #change variable for the next frame
        ax.clear()
        tcf = ax.tricontourf(t, u[:,framenumber], vmin=np.min(u), vmax=np.max(u), levels=lvls)
        return tcf,

    pam_ani = animation.FuncAnimation(fig, data_gen, fargs=(u, tcf),
                                  interval=60, blit=False, frames=u.shape[1] -1)


    plt.show();