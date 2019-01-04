#!/usr/env/bin python

class TangosInputHandler:
	''' Handler for extracting data from Tangos databases. At each timestep of interest, extracts the halo number of halos
		within some radius of the center halo. This script will default to following the main progenitor of the given central halo.
	Inputs:
		sim_name - Name of the simulation recognizable by Tangos. Type(str)
		center_halo - Number of the halo to be centered on, default halo ordering by number of particles. Type(int)
		latest_timestep - Lowest redshift timestep of interest. Type(int)
		earliest_timestep - Highest redshift timestep of interest. Type(int)
		region_size - Size of the surrounding halo search radius in comoving units. Type(float)
		boxsize - Simulation box size in comoving coordinates. Default 25e3, Type(float)
		periodic - Flag for simulation periodic boundary conditions. Default True
	'''

	def wrap(relpos, box_size=25e3):
		''' If simulation has periodic boundary conditions, periodically wrap coordinates relative to the center point.
		Inputs:
			relpos - (N,M) array of M-dimensional coordinates with a pre-defined center.
			boxsize - Size of simulation box in relpos units.
		Outputs:
			Replaces relpos in-place.
		'''
	    bphys = box_size
		bad = np.where(np.abs(relpos) > bphys/2.)
		if type(bphys) == np.ndarray:
		    relpos[bad] = -1.0 * (relpos[bad] / np.abs(relpos[bad])) * np.abs(bphys[bad] - np.abs(relpos[bad]))
		else:
		    relpos[bad] = -1.0 * (relpos[bad] / np.abs(relpos[bad])) * np.abs(bphys - np.abs(relpos[bad]))

	def __init__(self, sim_name, center_halo, latest_timestep, earliest_timestep, region_size, box_size=25e3, periodic=True):
		import numpy as np
		import tangos
		from scipy.spatial import KDTree

		self.simulation = tangos.get_simulation(str(sim_name))
		self.region_size = 5*float(region_size)
		self.box_size = float(box_size)
		self.timestep, self.num_timesteps = self.find_timesteps(int(latest_timestep), int(earliest_timestep))
		self.center_halo = self.timestep[center_halo]

		self.relevant_halos = self.gather_relevant_halos()
		
	def surrounding_halos(self):
		''' Finds all halos within a given radius of the central halo.
		Outputs:
			Tangos halo_number() of halos within a radius of region_size from the central halo.
		'''
		coords_center = self.center_halo.calculate('shrink_center')
		halo_num, coords = self.timestep.calculate_all('halo_num', 'shrink_center')
		coords_rel = coords - coords_center

		if self.periodic == True:
			wrap(coords_rel, self.box_size)

		tree = KDTree(coords_rel)
		query = tree.query_ball_point([0,0,0], self.region_size)
		return halo_num[query]

	def find_timesteps(self, earliest_timestep, latest_timestep):
		''' Gather all available timesteps from the simulation database and trim the list according to given bounds.
		Inputs:
			earliest_timestep - Highest redshift timestep number.
			latest_timestep - Lowest redshift timestep number.
		Outputs:
			Object containing highest redshift timestep desired.
			Number of steps visible to Tangos between earliest and latest desired timesteps.
		'''
		all_steps = [str(step) for step in self.simulation.timesteps]
		early = [i for i,step in enumerate(all_steps) if earliest_timestep in step][0]
		late = [i for i,step in enumerate(all_steps) if latest_timestep in step][0]
		return all_steps[late], late - early

	def gather_relevant_halos(self):
		''' Generate a list of halo numbers for all halos surrounding the central halo, for each timestep.
		Outputs:
			List of length [Number of timesteps], where each entry is a list of halos surrounding the central halo or central halo progenitor.
		'''
		relevant_halos = []
		for step in range(len(num_timesteps)):
			relevant_halos += [self.surrounding_halos()]
			self.timestep = self.timestep.previous
			self.center_halo = self.center_halo.previous
		return relevant_halos

	def major_progenitors(self, halo):
		return tangos.relation_finding.MultiHopMajorProgenitorsStrategy(halo).all()