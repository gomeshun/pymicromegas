#	LanHEP output produced at Mon Jul 15 17:09:46 2024
#	from the file '/home/belyaev/packages/lanhep/lanhep402/mdl/qcd.mdl'
#	Model named 'QCD'


from object_library import all_lorentz, Lorentz
from function_library import complexconjugate, re, im, csc, sec, acsc, asec

UUV_1 = Lorentz(name = 'UUV_1',
                  spins = [-1, -1, 3],
                  structure = ' P(3,2)')

FFV_1 = Lorentz(name = 'FFV_1',
                  spins = [2, 2, 3],
                  structure = 'Gamma(3,1,2)')

VVV_1 = Lorentz(name = 'VVV_1',
                  spins = [3, 3, 3],
                  structure = ' Metric(1,3)* P(2,3)')

VVV_2 = Lorentz(name = 'VVV_2',
                  spins = [3, 3, 3],
                  structure = ' Metric(2,3)* P(1,3)')

VVV_3 = Lorentz(name = 'VVV_3',
                  spins = [3, 3, 3],
                  structure = ' Metric(1,2)* P(3,1)')

VVV_4 = Lorentz(name = 'VVV_4',
                  spins = [3, 3, 3],
                  structure = ' Metric(1,3)* P(2,1)')

VVV_5 = Lorentz(name = 'VVV_5',
                  spins = [3, 3, 3],
                  structure = ' Metric(1,2)* P(3,2)')

VVV_6 = Lorentz(name = 'VVV_6',
                  spins = [3, 3, 3],
                  structure = ' Metric(2,3)* P(1,2)')

VVVV_1 = Lorentz(name = 'VVVV_1',
                  spins = [3, 3, 3, 3],
                  structure = ' Metric(1,3)* Metric(2,4)')

VVVV_2 = Lorentz(name = 'VVVV_2',
                  spins = [3, 3, 3, 3],
                  structure = ' Metric(1,4)* Metric(2,3)')

VVVV_3 = Lorentz(name = 'VVVV_3',
                  spins = [3, 3, 3, 3],
                  structure = ' Metric(1,2)* Metric(3,4)')

