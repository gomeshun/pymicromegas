#	LanHEP output produced at Mon Jul 15 17:09:46 2024
#	from the file '/home/belyaev/packages/lanhep/lanhep402/mdl/qcd.mdl'
#	Model named 'QCD'


from object_library import all_vertices, Vertex
import particles as P
import couplings as C
import lorentz as L

V_1 = Vertex(name='V_1',
              particles = [  P.Ggh__tilde__,  P.Ggh,  P.G ],
              lorentz = [ L.UUV_1  ],
              color = [ 'f(1,2,3)' ],
              couplings = {(0,0):C.GC_1})

V_2 = Vertex(name='V_2',
              particles = [  P.Q,  P.q,  P.G ],
              lorentz = [ L.FFV_1  ],
              color = [ 'T(3,2,1)' ],
              couplings = {(0,0):C.GC_2})

V_3 = Vertex(name='V_3',
              particles = [  P.G,  P.G,  P.G ],
              lorentz = [ L.VVV_1, L.VVV_2, L.VVV_3, L.VVV_4, L.VVV_5, L.VVV_6  ],
              color = [ 'f(1,2,3)' ],
              couplings = {(0,0):C.GC_3,(0,1):C.GC_1,(0,2):C.GC_3,(0,3):C.GC_1,(0,4):C.GC_1,(0,5):C.GC_3})

V_4 = Vertex(name='V_4',
              particles = [  P.G,  P.G,  P.G,  P.G ],
              lorentz = [ L.VVVV_1, L.VVVV_2, L.VVVV_3  ],
              color = [ 'f(1,2,-1)*f(3,4,-1)', 'f(1,3,-1)*f(2,4,-1)', 'f(1,4,-1)*f(2,3,-1)' ],
              couplings = {(0,0):C.GC_4,(2,0):C.GC_5,(0,1):C.GC_5,(1,1):C.GC_5,(1,2):C.GC_4,(2,2):C.GC_4})

