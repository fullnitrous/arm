import matplotlib
import matplotlib.pyplot as plt
import math
from matplotlib.widgets import Slider, Button
import time
from copy import deepcopy

from decelerate import decelerate
from accelerate import accelerate

def v_lim(x):
	#return [0.2 + 0.2*x, 0.2, 0]
	#return [0.2, 0, 0]
	return [
		0.2*math.sin(20*x) + 0.3*x + 0.14,
		4*math.cos(20*x) + 0.3,
		-80 * math.sin(20*x)
	]

def ffol(x, params):
	x0, x1, a0, v0, v_max, v_min, v_min2, a_max, j_max = params	
	v1, a1, _ = v_lim(x)
	_, v2, a2, _ = decelerate(x1 - x, a1, v1, v_min, a_max, j_max, True)
	return v2, (v1, a1, v2, a2)

def facc(x, params):
	x0, x1, a0, v0, v_max, v_min, v_min2, a_max, j_max = params	
	_, v1, a1, _ = accelerate(x - x0, a0, v0, v_max, a_max, j_max)
	_, v2, a2, _ = decelerate(x1 - x, a1, v1, v_min, a_max, j_max, True)
	return v2, (v1, a1, v2, a2)

def fdec(x, params):
	x0, x1, a0, v0, v_max, v_min, v_min2, a_max, j_max = params
	_, v1, a1, _ = decelerate(x - x0, a0, v0, v_min2, a_max, j_max)
	_, v2, a2, _ = decelerate(x1 - x, a1, v1, v_min, a_max, j_max, True)
	return v2, (v1, a1, v2, a2)

def bsolve(a, b, n, y, f, params):
	tol = 1e-7
	for i in range(n):
		x = (a + b) / 2.0
		fx, ox = f(x, params)
		fa,  _ = f(a, params)
		fx -= y
		fa -= y
		if fx == 0.0 or (-tol < fx and fx < tol):
			return True, x, ox
		elif fx*fa < 0: b = x
		else:           a = x
	return False, 0, None

def bsearch(arr, idx, n):
	lo = 0
	hi = n-1
	mi = (lo + hi) // 2	
	if idx == arr[0]:   return 0	
	if idx == arr[n-1]: return n-2
	i = 0
	while idx < arr[mi] or idx >= arr[mi+1]:
		if idx < arr[mi]: hi = mi
		else:             lo = mi	
		mi = (lo + hi) // 2
		i += 1
	return mi

# sanity check this fucking thing lol
# maybe integrate into algorithm main loop 
def minimas(v, n, dx):
	mx = [0]
	mv = []
	mg = 1e10
	for i in range(1, n-1, 1):
		if v[i] < mg: mg = v[i]
		a0 = (v[i+0] - v[i-1]) / dx
		a1 = (v[i+1] - v[i+0]) / dx
		is_maxima = a0 > 0 and a1 < 0
		is_minima = a0 < 0 and a1 > 0
		if is_maxima and len(mx) > 1:
			mx[len(mx)-1] = i*dx
		elif is_minima:
			mx.append(i*dx)
			mv.append(v[i])
			lm = i
	if len(mx) > 1:
		print("hello")
		mx[len(mx)-1] = (n-2)*dx + 0.5*dx
	mx.append((n-1)*dx)
	mv.append(0)
	return mx, mv

def cut_op(u, ops, idx, x1, a1, v1):
	u[idx+1]       = x1
	ops[idx]["a1"] = a1
	ops[idx]["v1"] = v1
	return idx

def push_op(u, ops, opidx, op, x1, a0, v0, a1, v1, f, vmin):
	if ops[opidx]["op"] == op and ops[opidx]["vmin"] == vmin and op != "fol":
		u[opidx+1]       = x1
		ops[opidx]["v1"] = v1
		ops[opidx]["a1"] = a1
	elif ops[opidx]["op"] == None:
		u[opidx+1]       = x1
		ops[opidx]["op"] = op
		ops[opidx]["v1"] = v1
		ops[opidx]["a1"] = a1
	else:
		opidx += 1
		u[opidx+1]         = x1
		ops[opidx]["op"]   = op
		ops[opidx]["a0"]   = a0
		ops[opidx]["v0"]   = v0
		ops[opidx]["a1"]   = a1
		ops[opidx]["v1"]   = v1
		ops[opidx]["vmin"] = vmin
	return opidx

def op2coeff(a, b, op, v_max, a_max, j_max):
	u = []
	c = []
	i2j_acc = [j_max, 0, -j_max, 0]
	i2j_dec = [-j_max, 0, j_max, 0]
	a0 = op["a0"]
	v0 = op["v0"]
	if op["op"] == "acc" or op["op"] == "dec":
		px0 = a
		if op["op"] == "acc": segments, _, _, _ = accelerate(0, a0, v0, v_max, a_max, j_max)
		else:                 segments, _, _, _ = decelerate(0, a0, v0, op["vmin"], a_max, j_max)
		
		for i in range(len(segments)):
			x0 = a + segments[i]
			x1 = a + min(segments[i+1], (b-a)) if i < 3 else b
			j0 = i2j_acc[i] if op["op"] == "acc" else i2j_dec[i]
			c0 = 0.5*j0
			c1 = a0 - j0*x0
			c2 = v0 - a0*x0 + 0.5*j0*x0**2
			a0 = 2*c0*x1 + c1
			v0 = c0*x1**2 + c1*x1 + c2
			
			if x1 - px0 > 1e-15: px0 = x1
			else:                continue

			u.append(x1)
			c.append([c0, c1, c2])
		return u, c
	elif op["op"] == "fol":
		u.append(b)
		c.append([None, None, None, None])
		return u, c

def _eval(x, u, ops, opidx, v_max, v_min, a_max, j_max):
	k = bsearch(u, x, opidx+2)
	if ops[k]["op"] == "acc":
		vlines, v, a, j = accelerate(x - u[k], ops[k]["a0"], ops[k]["v0"], v_max, a_max, j_max)
		vlines[0] += u[k]
		vlines[1] += u[k]
		vlines[2] += u[k]
		vlines[3] += u[k]
	elif ops[k]["op"] == "dec":
		vlines, v, a, j = decelerate(x - u[k], ops[k]["a0"], ops[k]["v0"], v_min, a_max, j_max)
		vlines[0] += u[k]
		vlines[1] += u[k]
		vlines[2] += u[k]
		vlines[3] += u[k]
	elif ops[k]["op"] == "fol":
		v, a, j = v_lim(x)
		vlines = []
	return vlines, v, a, j

def aisec(x, params):
	u, ops, opidx, v_min, v_max, a_max, j_max = params	
	_, v, a, j = _eval(x, u, ops, opidx, v_max, v_min, a_max, j_max)
	_, la0, _ = v_lim(x)
	return a - la0, None

def bracket_acc_isec(u, i0, i1, dx, opidx, ops, v_max, v_min, a_max, j_max, da1):
	_probes, _, _, _ = _eval((i1)*dx, u, ops, opidx, v_max, v_min, a_max, j_max)
	_probes.append(u[opidx])	
	_probes = sorted(_probes)
	probes = [i0*dx]
	for p in _probes:
		if i0*dx < p < i1*dx:
			probes.append(p)
	probes.append((i1)*dx)

	for k in range(len(probes)-1):
		a = probes[k]
		b = probes[k+1]
		
		_, la0, _ = v_lim(a)
		_, la1, _ = v_lim(b)
		_, _, a0, _ = _eval(a, u, ops, opidx, v_max, v_min, a_max, j_max)
		_, _, a1, _ = _eval(b, u, ops, opidx, v_max, v_min, a_max, j_max)

		if (la0-a0)*(la1-a1) < 0:
			params = (u, ops, opidx, v_min, v_max, a_max, j_max)
			solved, x, _ = bsolve(a, b, 50, 0, aisec, params)
			return x, (a, b)

	return 0, None

def backtrack(TARGETX, FUCKOFF, TARGETV, u, ops, opidx, v_min, v_max, a_max, j_max):
	for k in range(opidx, -1, -1):
		_, dvk, dak, _ = decelerate(TARGETX-u[k], ops[k]["a0"], ops[k]["v0"], v_min, a_max, j_max)	
		if dvk <= TARGETV:
			params = (u[k], TARGETX, ops[k]["a0"], ops[k]["v0"], v_max, v_min, ops[k]["vmin"], a_max, j_max)
			if   ops[k]["op"] == "acc": f = facc
			elif ops[k]["op"] == "dec": f = fdec
			elif ops[k]["op"] == "fol": f = ffol
			solved, x, ox = bsolve(u[k], TARGETX, 50, TARGETV, f, params)
			if not solved: continue
			dv0, da0, dv1, da1 = ox
			dv1 = max(dv1, v_min)
			opidx = cut_op(u, ops, k, x, da0, dv0)
			opidx = push_op(u, ops, opidx, "dec", FUCKOFF, da0, dv0, da1, dv1, None, v_min)		
			return opidx, (dv0, da0, dv1, da1)

def algorithm(n, v, a_max, j_max, v_max):	
	def constrained(_v, _a, _j):
		cond0 = 1e-7 <= _v and _v <= v_max
		cond1 = -a_max <= _a and _a <= a_max
		cond2 = -j_max <= _j and _j <= j_max
		return cond0 and cond1

	n_max = int(v_max)
	v_max = 2

	i      = 0
	dx     = 1 / (n-1)
	v[0]   = 0
	v[n-1] = 0
	opidx  = 0
	ops    = []
	u      = [0]*(2*n)
	mx, mv = minimas(v, n, dx)

	for _ in range(2*n): ops.append(deepcopy({"op": None, "a0": 0, "v0": 0, "a1": 0, "v1": 0, "vmin": -1}))

	auxx = [[],[]]
	auxy = [[],[]]
	
	while i < n_max - 1:
		v_min = min(mv[bsearch(mx, (i+1)*dx, len(mx))], v[i+1])

		if ops[opidx]["op"] == "acc":
			aa0 = ops[opidx]["a1"]
			aa0 = ops[opidx]["a0"]
			av0 = ops[opidx]["v0"]
			au0 = u[opidx]
		else:
			aa0 = ops[opidx]["a1"]
			av0 = ops[opidx]["v1"]
			au0 = u[opidx+1]
		
		_, av1, aa1, _ = accelerate((i+1)*dx - au0, aa0, av0, v_max, a_max, j_max)
		opidx = push_op(u, ops, opidx, "acc", (i+1)*dx, aa0, av0, aa1, av1, None, -1)

		if av1 <= v[i+1]:
			i += 1
			continue
		
		tmp_ops = deepcopy(ops)
		tmp_u   = deepcopy(u)

		opidx_tmp, (dv0, da0, dv1, da1) = backtrack((i+1)*dx, (i+1)*dx, v[i+1], tmp_u, tmp_ops, opidx, v_min, v_max, a_max, j_max)	
		lv, la, lj = v_lim((i+1)*dx)
	
		tmp2_ops = deepcopy(tmp_ops)
		tmp2_u = deepcopy(tmp_u)
		tmp2_opidx = opidx_tmp

		kissed = False

		max_its = 50 if constrained(*v_lim((i+1)*dx)) else 0
		
		# create a function for this thing
		for k in range(max_its):
			isecx, bracketed = bracket_acc_isec(tmp2_u, i, i+1, dx, tmp2_opidx, tmp2_ops, v_max, v_min, a_max, j_max, da1)

			if bracketed != None:
				#print(isecx)
				#auxx[1] = [bracketed[0], bracketed[0], bracketed[1], bracketed[1]]
				#auxy[1] = [0.3, -0.1, 0.3, -0.1]

				tmp2_opidx, (dv0, da0, dv1, da1) = backtrack(isecx, (i+1)*dx, v_lim(isecx)[0], tmp2_u, tmp2_ops, tmp2_opidx, v_min, v_max, a_max, j_max)

				_, _, da1, _ = _eval((i+1)*dx, tmp2_u, tmp2_ops, tmp2_opidx, v_max, v_min, a_max, j_max)
					
				lv, la, lj = v_lim(isecx)
				
				v_diff = abs(lv-tmp2_ops[tmp2_opidx]["v1"])
				a_diff = abs(la-tmp2_ops[tmp2_opidx]["a1"])
				

				if v_diff < 1e-7 and a_diff < 1e-7 and -j_max <= lj and lj <= j_max:
					#auxx[0] = [isecx, isecx]
					#auxy[0] = [0.25, -0.1]
					opidx, (dv0, da0, dv1, da1) = backtrack(isecx, isecx, v_lim(isecx)[0], u, ops, opidx, v_min, v_max, a_max, j_max)
					kissed = True
					i += 1	
					i0 = i	
					while i < n_max - 2:
						lv, la, lj = v_lim((i+1)*dx)
						if constrained(lv, la, lj):
							lv1, la1 = lv, la
							i += 1
						else:
							break
					if i0 < i:
						opidx = push_op(u, ops, opidx, "fol", i*dx, da1, dv1, la1, lv1, None, -1)		
					break
			else:
				break

		if not kissed:
			opidx = opidx_tmp
			u = tmp_u
			ops = tmp_ops	
			i += 1
		

	
	xv = []
	vv = []
	av = []
	jv = []
	
	U = [0]	
	C = []
	for i in range(opidx+1):
		_u, _c = op2coeff(u[i], u[i+1], ops[i], v_max, a_max, j_max)
		U += _u
		C += _c

	max_l = U[len(U)-1]
	for i in range(5000):
		x = i/4999
		if x > max_l: break
		k = bsearch(U, x, len(U))
		xv.append(x)
		if C[k][0] != None:
			vv.append(C[k][0]*x**2 + C[k][1]*x + C[k][2])
			av.append((2*C[k][0]*x + C[k][1]) / 50)
			jv.append((2*C[k][0]) / 800)
		else:
			_v, _a, _j = v_lim(x)
			vv.append(_v)
			av.append(_a / 50)
			jv.append(_j / 800)
	
	# * right after the final velocity profile has been generated, it needs to
	#   be verified
	#
	# * make it so that the baseline velocity limit is defined only by its 0th
	#   derivative, then piecewise chebyshev approximate it, trivial to then
	#   integrate or differentiate
	#
	# * for the sections that are followed, chebyshev approximate only that
	#   section, do not use the first approximation
	#
	# * optimize away the deepcopies, copying the entire state on each iteration
	#   should not be needed


	return 0, auxx, auxy, xv, vv, av, jv





def gen_profile(limf, parms, n, a_max, j_max, v_max):
	vlx = []
	vly = []
	vldy = []
	vld2y = []

	for i in range(n):
		x = i / (n-1)
		vlx.append(x)
		vly.append(limf(x)[0])
		vldy.append(limf(x)[1] / 50)
		vld2y.append(limf(x)[2] / 800)

	success, auxx, auxy, x, v, a, j = algorithm(n, vly, a_max, j_max, v_max)
	return success, auxx, auxy, vlx, vly, vldy, vld2y, x, v, a, j

if __name__  == "__main__":
	n     = 100
	a_max = 20
	j_max = 100
	v_max = 13

	fig, ax = plt.subplots()
	fig.dpi = 250
	parms = [0.88888888888, 80.0, 8.888888, 1]
	success, auxx, auxy, vlx, vly, vldy, vld2y, x, v, a, j = gen_profile(v_lim, parms, n, a_max, j_max, v_max)

	vprofile, = ax.plot(x, v, c="k")
	aprofile, = ax.plot(x, a, c="g")
	jprofile, = ax.plot(x, j, c="r")

	aux0, = ax.plot(auxx[0], auxy[0], c="y")
	aux1, = ax.plot(auxx[1], auxy[1], c="m")

	ax.plot([0, 1], [0, 0], c="b")

	limit, = ax.plot(vlx, vly, c="k", linestyle="--")
	limit1, = ax.plot(vlx, vldy, c="g", linestyle="--")
	limit2, = ax.plot(vlx, vld2y, c="r", linestyle="--")

	#ax.scatter(vlx, vly, c="k")

	fig.subplots_adjust(left=0.3, bottom=0.25)

	axa_max = fig.add_axes([0.25, 0.1, 0.65, 0.03])
	axa_max_slider = Slider(
		ax=axa_max,
		label="a_max",
		valmin=0,
		valmax=20,
		valinit=a_max
	)   

	axj_max = fig.add_axes([0.15, 0.25, 0.0225, 0.63])
	axj_max_slider = Slider(
		ax=axj_max,
		label="j_max",
		valmin=0.01,
		valmax=100,
		valinit=j_max,
		orientation="vertical"
	)
	
	axv_max = fig.add_axes([0.05, 0.25, 0.0225, 0.63])
	axv_max_slider = Slider(
		ax=axv_max,
		label="v_max",
		valmin=0,
		valmax=n,
		valinit=v_max,
		orientation="vertical"
	)

	def update(val):
		success, auxx, auxy, vlx, vly, _, _, x, v, a, j = gen_profile(v_lim, parms, n, axa_max_slider.val, axj_max_slider.val, axv_max_slider.val)
		
		vprofile.set_xdata(x)
		vprofile.set_ydata(v)

		aprofile.set_xdata(x)
		aprofile.set_ydata(a)
		
		jprofile.set_xdata(x)
		jprofile.set_ydata(j)
		
		aux0.set_xdata(auxx[0])
		aux0.set_ydata(auxy[0])

		aux1.set_xdata(auxx[1])
		aux1.set_ydata(auxy[1])
		
		limit.set_ydata(vly)
		fig.canvas.draw_idle()
	
	axa_max_slider.on_changed(update)
	axj_max_slider.on_changed(update)
	axv_max_slider.on_changed(update)
	
	plt.show()
