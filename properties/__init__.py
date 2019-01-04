class PropertyCalculation:
	def __init__(self, handler_name='tangos'):
		self.handler_name = handler_name
		if self.handler_name == 'tangos':
			from handlers import TangosInputHandler as handler

		if self.handler_name == 'pynbody'
			from handlers import PynbodyInputHandler as handler