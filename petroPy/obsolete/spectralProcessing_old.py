	# def baselineCorrect(self, birs, method = 'gcvspline', s = 0.001):
	# 	# Baseline correction with general cross-validated splines
	# 	# For n interpolation regions birs has to be a numpy array with shape (n,2), where each row is [lower_limit, upper_limit]
	# 	if self.norm == 2:
	# 		warn('run normalisation again to normalise baseline corrected spectrum')
	# 	spectrum = self.intensities[self.spectrumSelect]
	# 	self.intensities['BC'], self.baseline = blcorrect(x_input = self.x, y_input = spectrum, bir = birs, method = method, s = s)
	# 	self.intensities['BC'], self.baseline = self.intensities['BC'].reshape(-1), self.baseline.reshape(-1)
	# 	self.BC = 3
	# 	self.spectrumSelect = intensityDict[self.BC + self.norm]