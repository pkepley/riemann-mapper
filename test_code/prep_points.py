import numpy as np

def distribute_points_evenly(gamma, n):
	gamma_len = gamma.length()
	sep_len = gamma_len / n
	ps = np.zeros(n, dtype=np.complex64)
	gamma_ind_ps = np.zeros(n, dtype=np.int32)
	accum_len = 0.0
	p_j_len =  0.0
	j = 0
	i = 0
	for i, gamma_i in enumerate(gamma):
		gamma_i_len = gamma_i.length()
		len_along_gamma_i = p_j_len - accum_len

		# this is a hack
		while(len_along_gamma_i < gamma_i_len - .000001):
			param_p_j = len_along_gamma_i / gamma_i_len			
			p_j = gamma_i.point(param_p_j)
			ps[j] = p_j
			gamma_ind_ps = i

			# update:
			j = j + 1
			p_j_len = p_j_len + sep_len
			len_along_gamma_i = len_along_gamma_i + sep_len

		accum_len = accum_len + gamma_i_len

	return ps, gamma_ind_ps


if __name__ == '__main__':
	from svgpathtools import svg2paths, disvg, path
	import matplotlib.pyplot as plt
	# paths, attrs = svg2paths('../imgs/indiana_map.svg')
	paths, attrs = svg2paths('../imgs/circle.svg')

	gamma = paths[0]
	n = 200
	(ps, gamma_ind_ps) = distribute_points_evenly(gamma, n)
	xs = ps.real
	ys = -ps.imag
	plt.scatter(xs,ys)
	plt.axis('equal')
	plt.show()
