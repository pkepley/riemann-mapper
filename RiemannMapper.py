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
		I1 = m / (ws - self.a)
		I2 = self.points_on_gamma.cauchy_integral(self.L_a - m / (self.zs - self.a), ws)
		return I1 + I2


if __name__ == '__main__':
	from svgpathtools import svg2paths, disvg, path
	from PointsOnGamma import PointsOnGamma
	from KernelsKS import KernelsKS
	import matplotlib.pyplot as plt

	# plot_type = 'indiana'
	plot_type = 'square'
	# plot_type = 'ellipse'
	# plot_type = 'circle'

	# Distribute n points around a piece-wise continuous path
	if plot_type == 'indiana':
		paths, attrs = svg2paths('./imgs/indiana_map.svg')
		gamma = paths[0]
		n = 4000
		
		pog = PointsOnGamma(gamma, n)
		ps = pog.ps
		xs = ps.real
		ys = ps.imag
		
		# not guaranteed to fall inside of gamma, but set a the mean of gamma
		c = np.mean(ps)
		v = ps[n/4] - c
		r = abs(ps[0] - c)
		#a = c #+ 0.5 * v
		a = 160+190.j
		print c, a, r
		kks = KernelsKS(a, pog)	
		kks.solve_ks()
		
		
		rm = RiemannMapper(kks)
		print rm.evaluate_riemann(a+.001)
		
		for i in range(0,200,5):
			# ws = a + ((i+1) / 300.0) * (ps - a)
			y = 50 + (i/200.0)*(300 - 50)
			ws = np.array([ x+y*1.0j for x in range(50,275)])
			zz = rm.evaluate_riemann(ws)
			
			ws_xx, ws_yy = ws.real, ws.imag		
			plt.figure(1)
			plt.plot(ps.real, ps.imag)
			plt.scatter(ws_xx, ws_yy, marker='.', label="{0}".format(i))
			
			L_xx, L_yy = rm.La_on_ws.real, rm.La_on_ws.imag		
			plt.figure(2)
			plt.scatter(L_xx, L_yy, label="{0}".format(i))
			
			S_xx, S_yy = rm.Sa_on_ws.real, rm.Sa_on_ws.imag		
			plt.figure(3)
			plt.scatter(S_xx, S_yy, label="{0}".format(i))


			tt = np.linspace(0,2*np.pi,100)
			circ_xs, circ_ys = np.cos(tt), np.sin(tt)
			
			x_zz, y_zz = zz[0].real, zz[0].imag
			plt.figure(4)
			if i == 0 :
				plt.plot(circ_xs, circ_ys, label="circle")
			plt.plot(x_zz, y_zz, label="{0}".format(i))
			
		plt.legend()
		plt.show()

	# Distribute n points around a piece-wise continuous path
	if plot_type == 'square':
		paths, attrs = svg2paths('./imgs/square.svg')
		gamma = paths[0]
		n = 4000
		
		pog = PointsOnGamma(gamma, n)
		ps = pog.ps
		xs = ps.real
		ys = ps.imag
		
		# not guaranteed to fall inside of gamma, but set a the mean of gamma
		c = np.mean(ps)
		v = ps[n/4] - c
		r = abs(ps[0] - c)
		a = 0.25+0.25j
		print c, a, r
		kks = KernelsKS(a, pog)	
		kks.solve_ks()
				
		rm = RiemannMapper(kks)
		
		for i in range(0,201,1):
			# ws = a + ((i+1) / 300.0) * (ps - a)
			y = 0.0125 + (i/200.0) * (1.0 - 2.0 * 0.0125)
			ws = np.array([ x+y*1.0j for x in np.linspace(0.0125, 1.0 - 0.0125,100)])
			zz = rm.evaluate_riemann(ws)
			
			ws_xx, ws_yy = ws.real, ws.imag		
			plt.figure(1)
			plt.plot(ps.real, ps.imag)
			plt.scatter(ws_xx, ws_yy, marker='.', label="{0}".format(i))
			
			L_xx, L_yy = rm.La_on_ws.real, rm.La_on_ws.imag		
			plt.figure(2)
			plt.scatter(L_xx, L_yy, label="{0}".format(i))
			
			S_xx, S_yy = rm.Sa_on_ws.real, rm.Sa_on_ws.imag		
			plt.figure(3)
			plt.scatter(S_xx, S_yy, label="{0}".format(i))


			tt = np.linspace(0,2*np.pi,100)
			circ_xs, circ_ys = np.cos(tt), np.sin(tt)
			
			x_zz, y_zz = zz[0].real, zz[0].imag
			plt.figure(4)
			if i == 0 :
				plt.plot(circ_xs, circ_ys, label="circle")
			plt.plot(x_zz, y_zz, label="{0}".format(i))
			
		plt.legend()
		plt.show()


	# Distribute n points around a piece-wise continuous path
	if plot_type in ['ellipse', 'circle']:
		paths, attrs = svg2paths('./imgs/{0}.svg'.format(plot_type))
		gamma = paths[0]
		n = 1000
		
		pog = PointsOnGamma(gamma, n)
		ps = pog.ps
		xs = ps.real
		ys = ps.imag
			
		c = np.mean(ps)
		v = ps[0] - c
		r = abs(ps[0] - c)
		a = c
		kks = KernelsKS(a, pog)	
		kks.solve_ks()
				
		rm = RiemannMapper(kks)
		print a, rm.evaluate_riemann(a+.001)
		
		for i in range(0,200,5):
			ws = a + ((i+1) / 300.0) * (ps - a)
			zz = rm.evaluate_riemann(ws)
			
			ws_xx, ws_yy = ws.real, ws.imag		
			plt.figure(1)
			plt.plot(ws_xx, ws_yy, label="{0}".format(i))
			plt.axis('equal')
			#plt.scatter(ws_xx, ws_yy, marker='.', label="{0}".format(i))
			
			L_xx, L_yy = rm.La_on_ws.real, rm.La_on_ws.imag		
			plt.figure(2)
			plt.scatter(L_xx, L_yy, label="{0}".format(i))
			plt.axis('equal')

			S_xx, S_yy = rm.Sa_on_ws.real, rm.Sa_on_ws.imag		
			plt.figure(3)
			plt.scatter(S_xx, S_yy, label="{0}".format(i))
			plt.axis('equal')
			
			x_zz, y_zz = zz[0].real, zz[0].imag
			plt.figure(4)
			plt.plot(x_zz, y_zz, label="{0}".format(i))
			plt.axis('equal')
			
		plt.legend()
		plt.show()

