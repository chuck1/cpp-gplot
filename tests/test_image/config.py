import pbs2

e = pbs2.Executable(self, 'test_image', __file__)

e._test = True

e.add_dep('gplot')

#e.args.append('-lglpk')

self.parts.append(e)


