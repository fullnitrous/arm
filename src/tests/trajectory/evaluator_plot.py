import matplotlib
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

plt.style.use("dark_background")
plt.rcParams["grid.color"] = "#202020"

def create_plot(elev, azim):
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.w_xaxis.set_pane_color((0, 0, 0))
	ax.w_yaxis.set_pane_color((0, 0, 0))
	ax.w_zaxis.set_pane_color((0, 0, 0))
	ax.w_xaxis.line.set_color("#FF4A98")
	ax.tick_params(axis="x", colors="#FF4A98")
	ax.w_yaxis.line.set_color("#0AFAFA")
	ax.tick_params(axis="y", colors="#0AFAFA")
	ax.w_zaxis.line.set_color("#7F83FF")
	ax.tick_params(axis="z", colors="#7F83FF")
	ax.set_box_aspect([1, 1, 1])
	ax.elev = elev
	ax.azim = azim
	return fig, ax

def ldcsv(path):
    x = []
    y = []
    z = []
    with open(path, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line[:-1]
            xval, yval, zval = line.split(", ")
            x.append(float(xval))
            y.append(float(yval))
            z.append(float(zval))
    return (x, y, z)


px, py, pz = ldcsv("csv/points.csv");
gx, gy, gz = ldcsv("csv/linear.csv");
fig0, ax0 = create_plot(30, -130)
ax0.scatter(px, py, pz, c="w", alpha=1.0)
ax0.plot(gx, gy, gz, color="w");
fig0.savefig("img/linear.svg", bbox_inches="tight")
# will nuke
#fig0.show()
print("saved img/linear.svg")

gx, gy, gz = ldcsv("csv/cubic0.csv");
fig1, ax1 = create_plot(30, -130)
ax1.scatter(px, py, pz, c="w", alpha=1.0)
ax1.plot(gx, gy, gz, color="w");
fig1.savefig("img/cubic0.svg", bbox_inches="tight")
# will nuke
#fig1.show()
print("saved img/cubic0.svg")

gx, gy, gz = ldcsv("csv/cubic1.csv");
fig2, ax2 = create_plot(30, -130)
ax2.scatter(px, py, pz, c="w", alpha=1.0)
ax2.plot(gx, gy, gz, color="w");
fig2.savefig("img/cubic1.svg", bbox_inches="tight")
# will nuke
#fig2.show()
print("saved img/cubic1.svg")

gx, gy, gz = ldcsv("csv/quintic0.csv");
fig3, ax3 = create_plot(30, -130)
ax3.scatter(px, py, pz, c="w", alpha=1.0)
ax3.plot(gx, gy, gz, color="w");
fig3.savefig("img/quintic0.svg", bbox_inches="tight")
# will nuke
#fig3.show()
print("saved img/quintic0.svg")

gx, gy, gz = ldcsv("csv/quintic1.csv");
fig4, ax4 = create_plot(30, -130)
ax4.scatter(px, py, pz, c="w", alpha=1.0)
ax4.plot(gx, gy, gz, color="w");
fig4.savefig("img/quintic1.svg", bbox_inches="tight")
# will nuke
#fig4.show()
print("saved img/quintic1.svg")

px, py, pz = ldcsv("csv/arcp.csv");
gx, gy, gz = ldcsv("csv/arc.csv");
fig5, ax5 = create_plot(30, -130)
ax5.scatter(px, py, pz, c="w", alpha=1.0)
ax5.plot(gx, gy, gz, color="w");
fig5.savefig("img/arc.svg", bbox_inches="tight")
# will nuke
#fig5.show()
print("saved img/arc.svg")

px, py, pz = ldcsv("csv/bezierpoly.csv");
gx, gy, gz = ldcsv("csv/bezier.csv");
fig6, ax6 = create_plot(30, -130)
ax6.scatter(px, py, pz, c="w", alpha=1.0)
ax6.plot(px, py, pz, c="w", linestyle="dashed", alpha=1.0)
ax6.plot(gx, gy, gz, color="w");
fig6.savefig("img/bezier.svg", bbox_inches="tight")
# will nuke
#fig6.show()
print("saved img/bezier.svg")

px, py, pz = ldcsv("csv/bsplinepoly.csv");
gx, gy, gz = ldcsv("csv/bspline.csv");
fig7, ax7 = create_plot(24, 140)
ax7.scatter(px, py, pz, c="w", alpha=1.0)
ax7.plot(px, py, pz, c="w", linestyle="dashed", alpha=1.0)
ax7.plot(gx, gy, gz, color="w");
fig7.savefig("img/bspline.svg", bbox_inches="tight")
# will nuke
#fig7.show()
print("saved img/bspline.svg")

px, py, pz = ldcsv("csv/nurbspoly.csv");
gx, gy, gz = ldcsv("csv/nurbs.csv");
fig8, ax8 = create_plot(45, -115)
#ax8.scatter(px, py, pz, c="w", alpha=1.0)
#ax8.plot(px, py, pz, c="w", linestyle="dashed", alpha=1.0)
ax8.plot(gx, gy, gz, color="w");
fig8.savefig("img/nurbs.svg", bbox_inches="tight")
# will nuke
#fig8.show()
print("saved img/nurbs.svg")


# will nuke
#plt.show()
