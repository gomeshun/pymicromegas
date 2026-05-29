#     LanHEP output produced at Mon Jul 15 17:09:46 2024
#     from the file '/home/belyaev/packages/lanhep/lanhep402/mdl/qcd.mdl'
#     Model named 'QCD'

from object_library import all_couplings, Coupling
from function_library import complexconjugate, re, im, csc, sec, acsc, asec

GC_1 = Coupling(name = 'GC_1',
               value = 'gg',
               order = { })

GC_2 = Coupling(name = 'GC_2',
               value = '-complex(0,1)*gg',
               order = { })

GC_3 = Coupling(name = 'GC_3',
               value = '-gg',
               order = { })

GC_4 = Coupling(name = 'GC_4',
               value = '-complex(0,1)*(gg*gg)',
               order = { })

GC_5 = Coupling(name = 'GC_5',
               value = 'complex(0,1)*(gg*gg)',
               order = { })

