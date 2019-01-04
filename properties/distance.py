from properties import PropertyCalculation

class Distance(PropertyCalculation):
	name = 'distance'

	def calculate(self, handler='tangos'):
		X, Y, Z = self.center_halo.calculate('shrink_center')
		return np.sqrt(X**2 + Y**2 + Z**2)