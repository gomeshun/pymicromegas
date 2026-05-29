
#include"../include/micromegas.h"
#include"../include/micromegas_aux.h"

extern REAL MassArray(int id,  int i);
extern REAL MixMatrix(int id, int i,int j);


static void multMC(int dim,COMPLEX*M1,COMPLEX*M2,COMPLEX*M3) 
{
  for(int i=0;i<dim;i++) for(int j=0;j<dim;j++)
  { COMPLEX R=0;
    for(int k=0;k<dim;k++) R+=M1[i*dim+k]*M2[k*dim+j];
    M3[i*dim+j]=R;
  }  
}



int  MtoSim(int dim, int i, int j) 
{ if(j<i) { int jm=j;j=i;i=jm;}  
  return   ((i)*dim-((i)*((i)+1))/2 +(j));
}  



int main(void)
{
   REAL  MTestSym[6] = {0., 1., 1.,
                            0., 1.,
                                0.};
printf("Diagonalization of symmetric  matrix:\n"
"          0  1  1\n"
"          1  0  1\n"
"          1  1  0\n" );

   int dim=3;

   printf("\n  For   func*.mdl\n");
   REAL m[3][3];
   m[0][0]=MTestSym[0];
   m[0][1]=MTestSym[1];
   m[0][2]=MTestSym[2];
   m[1][1]=MTestSym[3];
   m[1][2]=MTestSym[4];
   m[2][2]=MTestSym[5];

   initDiagonal();
   int id=rDiagonal(dim,m[0][0],m[0][1],m[0][2],m[1][1],m[1][2],m[2][2]);

   printf("Result: masses"); for(int i=1;i<=dim;i++)  printf(" %E ",(double)MassArray(id,i)); printf("\n");
   printf("Rotation matrix:\n");
   for(int i=1;i<=3;i++)
   { for(int j=1;j<=3;j++) printf("%12.4E  ",(double)MixMatrix(id,i,j)); printf("\n");}

   printf("\n Restoration of initial matrix\n");

   for(int i=1;i<=dim;i++) for(int j=1; j<=dim;j++)
   { REAL mm=0;
     for(int k=1;k<=dim;k++) mm+=MixMatrix(id,k,j)*MassArray(id,k)*MixMatrix(id,k, i);
      printf(" %12.4E ", (double) mm); if(j==dim) printf("\n");
   }
// Below is an example of the function to diagonalize the mass matric MTestSym
// when not working with a matrix defined in a model file func*.mdl
printf("\n  When not using model files :\n");

   REAL E[dim],  V[dim*dim];
   rJacobi(MTestSym,dim, E , V);
   printf("Result  Masses: "); for(int k=0;k<dim;k++) printf("%E ",(double)E[k]); printf("\n");
    printf("Rotation matrix:\n");

    for(int i=0;i<3;i++)
    { for(int j=0;j<3;j++) printf("%12.4E  ",(double)V[i*dim+j]); printf("\n");}

printf("//========    Neutrino 6x6 mass matrix ======\n");

REAL  phi1=1,   phi2=2,   phi3=3;      // rotation angles 
REAL  Acp=1,    A1=2,    A2=3;         // phases
REAL  m1=1E-12, m2=1E-11, m3=1E-10;    // left  masses
REAL  M1=1E7 ,  M2=1E8,   M3=1E9;      // masses of sterile neutrino

COMPLEX  R12[9]=  { Cos(phi3), Sin(phi3),  0,
                   -Sin(phi3), Cos(phi3),  0,
                    0,        0,           1};

COMPLEX  R13[9]=  { Cos(phi2),                       0,  Sin(phi2)*(Cos(Acp)+I*Sin(Acp)),
                     0,                              1,        0,
                   -Sin(phi2)*(Cos(Acp)-I*Sin(Acp)), 0,  Cos(phi2)  
               };
COMPLEX  R23[9]=  { 1,       0,        0,
                    0,   Cos(phi1),  Sin(phi1),
                    0,  -Sin(phi1),  Cos(phi1)  
               };
COMPLEX FF[9]  = { Cos(A1)+I*Sin(A1), 0,                 0,
                       0,             Cos(A2)+I*Sin(A2), 0,
                       0,                    0,          1};
                       
COMPLEX  MM[9] =  { Sqrt(m1*M1),   0,         0,
                    0,        Sqrt(m2*M2),    0, 
                    0,          0,         Sqrt(m3*M3)
               };  
               
                  
                    
COMPLEX MUL[9],MUL_[9];
 
     multMC(3,R12,R13,MUL);
     multMC(3,MUL,R23,MUL_);
     multMC(3,MUL_,FF,MUL);
     multMC(3,MUL,MM,MUL_);


COMPLEX Msum[36]; 
     for(int i=0;i<3;i++) for(int j=0;j<3;j++) 
     {  Msum[6*i+j]=0; Msum[6*(i+3)+j+3]=0;
        Msum[6*i+j+3]=MUL_[i*3+j]; Msum[6*(i+3)+j]=MUL_[j*3+i];
     }
     
     Msum[7*3]=M1; Msum[7*4]=M2;Msum[7*5]=M3;

COMPLEX MsumS[21];


for(int i=0; i<6;i++) for(int j=i;j<6;j++)   MsumS[MtoSim(6,i,j)]=Msum[i*6+j]; 

COMPLEX VV[36];
REAL EV[6];


int err=cJacobiS(MsumS, 6 , EV, VV);     
printf("err=%d  masses: ",err); for(int i=0;i<6;i++) printf(" %E ",(double)EV[i]); printf("\n"); 

}
