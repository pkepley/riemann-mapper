import numpy as np

class RiemannMapper:
	def __init__(self, problem_ks):
		self.points_on_gamma = problem_ks.points_on_gamma
		self.problem_ks = problem_ks
		self.a = self.problem_ks.a
		self.S_a = self.problem_ks.S_a
		self.L_a = self.problem_ks.L_a
		self.n_zs = self.points_on_gamma.n_ps
		self.zs = self.points_on_gamma.ps.reshape((self.n_zs,1))

	def evaluate_riemann(self, ws):
		''' 
		   compute Riemann map at the points ws. i.e. f(ws) 
		'''
		self.ws = ws
		self.Sa_on_ws = self.evaluate_szego(self.ws)
		self.La_on_ws = self.evaluate_garabedian(self.ws)
		self.f_on_ws = self.Sa_on_ws /  self.La_on_ws
		return self.f_on_ws 

	def evaluate_szego(self, ws):
		return self.points_on_gamma.cauchy_integral(self.S_a, ws)

	def evaluate_garabedian(self, ws):
		m = (1.0/(2.0*np.pi))
		I1 = m / (ws - a)
		I2 = self.points_on_gamma.cauchy_integral(self.L_a - m / (self.zs - self.a), ws)
		return I1 + I2


if __name__ == '__main__':
	from svgpathtools import svg2paths, disvg, path
	from PointsOnGamma import PointsOnGamma
	from KernelsKS import KernelsKS
	import matplotlib.pyplot as plt
	paths, attrs = svg2paths('./imgs/indiana_map.svg')
	#paths, attrs = svg2paths('./imgs/circle.svg')
	debug_plots = True

	# Distribute n points around a piece-wise continuous path
	gamma = paths[0]
	n = 1000

	pog = PointsOnGamma(gamma, n)
	ps = pog.ps
	xs = ps.real
	ys = ps.imag

	# not guaranteed to fall inside of gamma, but set a the mean of gamma
	c = np.mean(ps)
	r = abs(ps[0] - c)
	a = c - 0.05 * (ps[0] - c)
	kks = KernelsKS(a, pog)	
	kks.solve_ks()

	rm = RiemannMapper(kks)
	for i in range(0,8):
		ws = a + ((i+1) / 11.0) * (ps - a)
		ws = ws[np.where(abs(ws - c) < r)]
		zz = rm.evaluate_riemann(ws)

		L_xx, L_yy = rm.La_on_ws.real, rm.La_on_ws.imag		
		plt.figure(1)
		plt.scatter(L_xx, L_yy, label="{0}".format(i))

		S_xx, S_yy = rm.Sa_on_ws.real, rm.Sa_on_ws.imag		
		plt.figure(2)
		plt.scatter(S_xx, S_yy, label="{0}".format(i))

		x_zz, y_zz = zz[0].real, zz[0].imag
		plt.figure(3)
		plt.plot(x_zz, y_zz, label="{0}".format(i))

	plt.legend()
	plt.show()
