#	LanHEP output produced at Mon Jul 15 17:09:46 2024
#	from the file '/home/belyaev/packages/lanhep/lanhep402/mdl/qcd.mdl'
#	Model named 'QCD'


from __future__ import division
from object_library import all_particles, Particle
import propagators as Prop

Ggh =  Particle(  pdg_code = 0,
              name = 'G.c',
              antiname = 'G.C',
              spin = -1,
              color = 8,
              mass = 'ZERO',
              width = 'ZERO',
              texname = 'G.c',
              antitexname = 'G.C',
              charge = 0,
              LeptonNumber = 0,
              GhostNumber = 1)

Ggh__tilde__ = Ggh.anti()

G =  Particle(  pdg_code = 21,
              name = 'G',
              antiname = 'G',
              spin = 3,
              color = 8,
              mass = 'ZERO',
              width = 'ZERO',
              texname = 'G',
              antitexname = 'G',
              charge = 0,
              LeptonNumber = 0,
              GhostNumber = 0)

q =  Particle(  pdg_code = 1000501,
              name = 'q',
              antiname = 'Q',
              spin = 2,
              color = 3,
              mass = 'Mq',
              width = 'ZERO',
              texname = 'q',
              antitexname = 'Q',
              charge = 0,
              LeptonNumber = 0,
              GhostNumber = 0)

Q = q.anti()

