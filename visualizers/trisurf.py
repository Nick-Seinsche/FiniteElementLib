import plotly.figure_factory as ff
import numpy as np

filename_c4n = "meshes/mytriang_c4n.tsv";
c4n = np.loadtxt(filename_c4n, delimiter="\t")

filename_n4e = "meshes/mytriang_n4e.tsv";
n4e = np.loadtxt(filename_n4e, delimiter="\t").astype(int)

filename_sol = "solutions/mysol.tsv";
u = np.loadtxt(filename_sol, delimiter="\t")

fig = ff.create_trisurf(x=c4n[:,0], y=c4n[:,1], z=u, simplices=n4e,
                        title="Solution", aspectratio=dict(x=1, y=1, z=1))
fig.show()