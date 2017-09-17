import numpy as np

class ProblemKS:
	def __init__(self, a, points_on_gamma):
		# a must be inside of gamma, we should add error checking!
		self.points_on_gamma = points_on_gamma
		self.a = a
		
		# note, this will call set_a() to set the rhs for the KS problem.
		self.set_points_on_gamma(points_on_gamma)
		
	def set_points_on_gamma(self, points_on_gamma):
		self.points_on_gamma = points_on_gamma
		self.set_eval_points()
		self.set_cauchy_ca()
		self.set_KS_kernel()

	def set_a(self, a):
		self.a = a
		self.set_cauchy_ca()
		
	def set_eval_points(self):
		self.n_zs = self.points_on_gamma.n_ps
		self.zs = self.points_on_gamma.ps
		self.z_prime_s = self.points_on_gamma.ps_tangent
		self.Ts = self.points_on_gamma.ps_unit_tangent

	def set_cauchy_ca(self):
		# Cauchy Kernel, evaluated at z = ps and w = a
		self.C_a = 	(1.0/(2.0*np.pi*1.0j)) * (
			(self.Ts.conjugate()) / ((self.a - self.zs).conjugate()))

		self.c = np.sqrt(abs(self.z_prime_s)) * self.C_a

		# reshape
		self.C_a = self.C_a.reshape((self.n_zs,1))
		self.c = self.c.reshape((self.n_zs,1))
		
	def set_KS_kernel(self):
		w = self.zs.reshape((1,self.n_zs))
		Tw = self.Ts.reshape((1,self.n_zs))
		abs_z_prime_s = abs(self.z_prime_s).reshape((1,self.n_zs))

		z = self.zs.reshape((self.n_zs,1))
		Tz = self.Ts.reshape((self.n_zs,1))
		abs_w_prime_s = abs(self.z_prime_s).reshape((self.n_zs,1))

		# compute matrix for the KS kernel, A
		A0 = np.zeros((self.n_zs, self.n_zs))
		A0 = w - z		
		np.fill_diagonal(A0,1) # don't divide by zero
		A0 = 1.0 / A0
		np.fill_diagonal(A0,0) # K will vanish along diagonal, so set this to 0

		A1 = np.zeros((self.n_zs, self.n_zs))
		A1 = np.dot(A0, np.diag(Tw.flatten()))

		A2 = np.zeros((self.n_zs, self.n_zs))
		A2 = np.dot(np.diag(Tz.flatten()), A0)
		A2 = A2.conjugate()

		self.A = (1.0/(2.0*np.pi*1.0j)) *  (A1 - A2)	

		self.K = ((np.diag(abs_z_prime_s.flatten()) * self.A) * 
				  np.diag(abs_w_prime_s.flatten()))
		


if __name__ == '__main__':
	from svgpathtools import svg2paths, disvg, path
	from PointsOnGamma import PointsOnGamma
	import matplotlib.pyplot as plt
	paths, attrs = svg2paths('./imgs/indiana_map.svg')
	#paths, attrs = svg2paths('./imgs/circle2.svg')

	# Distribute n points around a piece-wise continuous path
	gamma = paths[0]
	n = 2000

	pog = PointsOnGamma(gamma, n)
	ps = pog.ps
	xs = ps.real
	ys = ps.imag

	# not guaranteed to fall inside of gamma, but set a the mean of gamma
	a = np.mean(xs) + np.mean(ys) * (1.0j)		
	kks = ProblemKS(a, pog)

	plt.figure()
	plt.plot(abs(kks.c))
	plt.plot(abs(kks.C_a))

	plt.figure()
	plt.imshow(abs(kks.A), interpolation='nearest')
	plt.colorbar()
	plt.show()
