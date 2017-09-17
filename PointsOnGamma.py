import numpy as np

class  PointsOnGamma:
	def __init__(self, gamma, n_ps):
		self.gamma = gamma
		self.n_ps = n_ps
		self.gamma_len = self.gamma.length()
		self.ps = np.zeros(self.n_ps, dtype=np.complex64)
		self.ps_seg_ind = np.zeros(self.n_ps, dtype=np.int32)
		self.ps_seg_param = np.zeros(self.n_ps, dtype=np.float64)
		self.ps_tangent = np.zeros(self.n_ps, dtype=np.complex64)
		self.ps_unit_tangent = np.zeros(self.n_ps, dtype=np.complex64)

		# For now, this is the only option to distribute points
		self.distribute_points_evenly()

	def distribute_points_evenly(self):		
		sep_len = self.gamma_len / self.n_ps
		accum_len = 0.0
		p_j_len =  0.0
		j = 0
		i = 0
		for i, gamma_i in enumerate(self.gamma):
			gamma_i_len = gamma_i.length()
			len_along_gamma_i = p_j_len - accum_len
			
			# this is a hack
			while(len_along_gamma_i < gamma_i_len - .000001):
				param_p_j = len_along_gamma_i / gamma_i_len			
				p_j = gamma_i.point(param_p_j)
				
				# save point information:
				self.ps[j] = p_j
				self.ps_seg_ind[j] = i
				self.ps_seg_param[j] = param_p_j
				self.ps_tangent[j] = gamma_i.derivative(param_p_j)
				self.ps_unit_tangent[j] = gamma_i.unit_tangent(param_p_j)
				
				# update:
				j = j + 1
				p_j_len = p_j_len + sep_len
				len_along_gamma_i = len_along_gamma_i + sep_len
				
			accum_len = accum_len + gamma_i_len
		
# Demo:
if __name__ == '__main__':
	from svgpathtools import svg2paths, disvg, path
	import matplotlib.pyplot as plt
	paths, attrs = svg2paths('./imgs/indiana_map.svg')
	#paths, attrs = svg2paths('./imgs/circle.svg')

	# Distribute n points around a piece-wise continuous path
	gamma = paths[0]
	n = 100

	pog = PointsOnGamma(gamma, n)
	ps = pog.ps

	# Plot
	xs = ps.real
	ys = -ps.imag
	plt.scatter(xs,ys)
	plt.axis('equal')
	plt.show()
