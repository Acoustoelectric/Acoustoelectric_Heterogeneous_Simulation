# -*- coding: utf-8 -*-

# Read in the saved file. 
inputs = []
trivial_producer = analysis.extractors.DataCacheImporter(inputs=inputs)
trivial_producer.Name = "AE Complex Source"
trivial_producer.FileName = u"c:\\Users\\jeantoul\Desktop\\sim4life_tests\\ATI\\scripts\\E_C_source.cache"
trivial_producer.UpdateAttributes()
document.AllAlgorithms.Add(trivial_producer)

# Now interpolate it to be on the same grid as the transient thermal diffusion solver. 

# now update this to be the source. 