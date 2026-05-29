#     LanHEP output produced at Mon Jul 15 17:09:46 2024
#     from the file '/home/belyaev/packages/lanhep/lanhep402/mdl/qcd.mdl'
#     Model named 'QCD'

from object_library import all_parameters, Parameter
from function_library import complexconjugate, re, im, csc, sec, acsc, asec

ZERO = Parameter(name = 'ZERO',
                  nature = 'internal',
                  type = 'real',
                  value = '0.0',
                  texname = '0')

Sqrt2 = Parameter(name = 'Sqrt2',
                  nature = 'internal',
                  type = 'real',
                  value = 'cmath.sqrt(2.0)',
                  texname = '\\sqrt{2}')

gg = Parameter(name = 'gg',
               nature = 'external',
               type = 'real',
               value = 1.13,
               texname = '\\text{gg}',
               lhablock = 'USERDEF',
               lhacode = [ 1 ] )

Mq = Parameter(name = 'Mq',
               nature = 'external',
               type = 'real',
               value = 0.02,
               texname = '\\text{Mq}',
               lhablock = 'USERDEF',
               lhacode = [ 2 ] )

