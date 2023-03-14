import argparse as ap
import numpy as np

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
    
    import plotly.figure_factory as ff
    import plotly.graph_objects as go

    fig = ff.create_trisurf(x=c4n[:,0], y=c4n[:,1], z=u, simplices=n4e,
                            title="Solution", aspectratio=dict(x=1, y=1, z=1), gridcolor="rgb(50, 50, 50)")
    fig.show()
    '''
    from matplotlib import pyplot as plt
    from matplotlib import cm

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    t = ax.plot_trisurf(c4n[:,0], c4n[:,1], n4e, u, cmap=cm.coolwarm, linewidth=0.1);
    plt.show();
    '''

elif len(u.shape) == 2:
    from matplotlib import pyplot as plt
    import matplotlib.animation as animation
    from matplotlib import cm

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    t = ax.plot_trisurf(c4n[:,0], c4n[:,1], n4e, u[:,0], cmap=cm.coolwarm, linewidth=0.1);

    def data_gen(framenumber, u, t):
        #change variable for the next frame
        ax.clear()
        t = ax.plot_trisurf(c4n[:,0], c4n[:,1], n4e, u[:,framenumber], cmap=cm.coolwarm, linewidth=0.1);
        ax.axes.set_xlim3d(left=0, right=1) 
        ax.axes.set_ylim3d(bottom=0, top=1) 
        ax.axes.set_zlim3d(bottom=0, top=1) 
        return t,

    pam_ani = animation.FuncAnimation(fig, data_gen, fargs=(u, t),
                                  interval=60, blit=False, frames=u.shape[1] -1)


    plt.show();