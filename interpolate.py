from numpy import interp 

def interpolate(prop, time, inter_pts, **kwargs):
	inter_vals = interp(inter_pts, time, prop, kwargs)
	return inter_pts, inter_vals

