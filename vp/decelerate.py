import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import math

def decelerate(x, a0, v0, v_min, a_max, j_max, special=False):
	if 2*a0**2 + 4*j_max*(v0 - v_min) < 0:
		if a0 <= 0:
			#print("impossible deceleration: cannot compensate j_max")
			return [0]*4, 0, 0, 0
		j_max = (a0**2) / (2*(v_min - v0))
		u0 = (2*a0) / (2*j_max)
		u1 = u0
		u2 = u0
	else:	
		x0 = (2*a0 + math.sqrt(2*a0**2 + 4*j_max*(v0 - v_min))) / (2*j_max)
		"""
		if x0 < 0:
			print("impossible deceleration: acceleration too low")
			return [0]*4, 0, 0, 0
		"""
		x0_lim = (a0 + a_max) / (j_max)
		s0     = min(x0, x0_lim)
		a1     = a0 - j_max*s0
		v1     = v0 + a0*s0 - 0.5*j_max*s0**2
		s2     = -a1 / j_max
		vr     = v1 + a1*s2 + 0.5*j_max*s2**2 - v_min
		s1     = vr / a_max
		v2     = v1 - a_max*s1
		u0     = s0
		u1     = u0+s1
		u2     = u1+s2

	vlines = [0, u0, u1, u2]
	
	if 0 <= x <= u0:
		v = v0 + (a0 - 0.5*j_max*x)*x
		a = a0 - j_max*x
		j = -j_max
	elif u0 <= x <= u1:
		v = v1 - a_max*(x-u0)
		a = -a_max
		j = 0
	elif u1 <= x <= u2:
		v = v2 + (a1 + 0.5*j_max*(x-u1))*(x-u1)
		a = a1 + j_max*(x-u1)
		j = j_max
	else:
		if special:
			v = v_min - 0.5*j_max*(x-u2)**2
			a = 0
			j = -j_max
		else:
			v = v_min
			a = 0
			j = 0
	
	return vlines, v, a, j

def gen_decelerate(n, a0, v0, v_min, a_max, j_max):	
	dx = 1 / (n-1)

	vlines, _, _, _ = decelerate(0, a0, v0, v_min, a_max, j_max)

	x  = []
	vy = []
	ay = []
	jy = []
	
	d = 0
	for i in range(2*n):
		x.append(d)
		_, v, a, j = decelerate(d, a0, v0, v_min, a_max, j_max)
		vy.append(v)
		ay.append(a)
		jy.append(j)
		d += dx
	
	return vlines, x, vy, ay, jy

if __name__ == "__main__":
	a0 = 0
	v0 = 3
	a_max = 2
	j_max = 5
	v_min = -2

	fig, ax = plt.subplots()
	fig.dpi = 200
	vlines, x, vy, ay, jy = gen_decelerate(300, a0, v0, v_min, a_max, j_max)

	vgraph, = ax.plot(x, vy, c="k")
	agraph, = ax.plot(x, ay, c="g")
	jgraph, = ax.plot(x, jy, c="r")

	vl0 = ax.axvline(x=vlines[0], color="k", linestyle="--")
	vl1 = ax.axvline(x=vlines[1], color="k", linestyle="--")
	vl2 = ax.axvline(x=vlines[2], color="k", linestyle="--")
	vl3 = ax.axvline(x=vlines[3], color="k", linestyle="--")

	fig.subplots_adjust(left=0.25, bottom=0.25)
	axa0 = fig.add_axes([0.25, 0.1, 0.65, 0.03])
	axa0_slider = Slider(
		ax=axa0,
		label="a0",
		valmin=-a_max,
		valmax=a_max,
		valinit=a0
	)

	axv0 = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
	axv0_slider = Slider(
		ax=axv0,
		label="v0",
		valmin=0,
		valmax=1,
		valinit=v0,
		orientation="vertical"
	)

	def update(val):
		vlines, x, vy, ay, jy = gen_decelerate(300, axa0_slider.val, axv0_slider.val, v_min, a_max, j_max)
		vgraph.set_ydata(vy)
		agraph.set_ydata(ay)
		jgraph.set_ydata(jy)
		vl0.set_xdata([vlines[0]])
		vl1.set_xdata([vlines[1]])
		vl2.set_xdata([vlines[2]])
		vl3.set_xdata([vlines[3]])
		fig.canvas.draw_idle()
	
	#plt.savefig("img/mtvp_decelerate.png", dpi=300, bbox_inches="tight")

	axa0_slider.on_changed(update)
	axv0_slider.on_changed(update)

	plt.show()
