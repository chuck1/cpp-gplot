import os
import pbs

l = pbs.Library(self, 'gplot', __file__)

l.doc_out_dir = "/media/sf_P_DRIVE/html/gplot"

self.parts.append(l)

self.execfile(os.path.join(__dir__, 'tests/config.py'))


