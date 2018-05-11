import pbs

e = pbs.Executable(self, 'test_image', __file__)

e._test = True

e.add_dep('gplot')

e.args.args.append('-lboost_program_options')
e.args.args.append('-lncurses')

self.parts.append(e)


