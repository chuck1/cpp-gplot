import pbs

e = pbs.Executable(self, 'test', __file__)

e._test = True

e.add_dep('gplot')

#e.args.append('-lglpk')

self.parts.append(e)


