      SUBROUTINE COUPLINGFITS(FLAG)

*   Subroutine to compute the LHC constraints on the SM Higgs couplings in
*   the kappa framework from ATLAS 2207.00092 Fig.6, CMS 2207.00043 Fig.4a
*   flag=1: B_invisible = 0,  kappa_V possibly > 1 (default)
*   flag=2: B_invisible free, B_undetected > 0, kappa_V <= 1
*   Formulas for combined uncertainties:
*   See http://pdg.lbl.gov/2018/reviews/rpp2018-rev-rpp-intro.pdf page 6

      IMPLICIT NONE

      INTEGER I,FLAG

      DOUBLE PRECISION k(2),wm(2),sigmam(2),wp(2),sigmap(2)
      DOUBLE PRECISION ksum,wmsum,wpsum,CMIN(7),CMAX(7)

      COMMON/CMIMAX/CMIN,CMAX


      IF(FLAG.EQ.0)RETURN


************************************************************************
*   k_V

      k(1)=0.99252d0		! k_z from ATLAS
      k(2)=1.04d0		! k_z from CMS

*   sigma-

      sigmam(1)=0.05652d0	! sigam-(k_z) from ATLAS
      sigmam(2)=0.07d0		! sigam-(k_z) from CMS

*   sigma+ is not necessary since k_V > 1 is impossible in the NMSSM

*   Combination

      wmsum=0d0
      wpsum=0d0
      ksum=0d0
      DO I=1,2
        wm(i)=1d0/sigmam(i)**2
        wmsum=wmsum+wm(i)
        ksum=ksum+wm(i)*k(i)
      ENDDO
      ksum=ksum/wmsum
      CMIN(1)=ksum-2d0/dsqrt(wmsum)
      CMAX(1)=1d0

************************************************************************
*   k_top

      IF(FLAG.EQ.1)THEN
       k(1)=0.94394d0		! from ATLAS
       k(2)=1.01d0		! from CMS
      ELSE
       k(1)=0.93723d0		! from ATLAS
       k(2)=1.01d0		! from CMS
      ENDIF

c   sigma+

      IF(FLAG.EQ.1)THEN
       sigmap(1)=0.11139d0	! from ATLAS
       sigmap(2)=0.11d0		! from CMS
      ELSE
       sigmap(1)=0.11108d0	! from ATLAS
       sigmap(2)=0.10d0		! from CMS
      ENDIF

c   sigma-

      IF(FLAG.EQ.1)THEN
       sigmam(1)=0.10907d0	! from ATLAS
       sigmam(2)=0.10d0 	! from CMS
      ELSE
       sigmam(1)=0.10904d0	! from ATLAS
       sigmam(2)=0.10d0 	! from CMS
      ENDIF

*   Combination

      wmsum=0d0
      wpsum=0d0
      ksum=0d0
      DO I=1,2
       wp(i)=1d0/sigmap(i)**2
       wpsum=wpsum+wp(i)
       wm(i)=1d0/sigmam(i)**2
       wmsum=wmsum+wm(i)
      ENDDO
      ksum=(wm(1)+wp(1))/2d0*k(1)+(wm(2)+wp(2))*k(2)/2d0
      ksum=ksum/((wm(1)+wp(1))/2d0+(wp(2)+wm(2))/2d0)
      CMIN(2)=ksum-2d0/dsqrt(wmsum)
      CMAX(2)=ksum+2d0/dsqrt(wpsum)

************************************************************************
*   k_bottom

      IF(FLAG.EQ.1)THEN
       k(1)=0.89158d0		! from ATLAS
       k(2)=0.99d0		! from CMS
      ELSE
       k(1)=0.81941d0		! from ATLAS
       k(2)=0.90d0		! from CMS
      ENDIF

c   sigma+

      IF(FLAG.EQ.1)THEN
       sigmap(1)=0.11399d0	! from ATLAS
       sigmap(2)=0.17d0		! from CMS
      ELSE
       sigmap(1)=0.08544d0	! from ATLAS
       sigmap(2)=0.10d0		! from CMS
      ENDIF

c   sigma-

      IF(FLAG.EQ.1)THEN
       sigmam(1)=0.10709d0	! from ATLAS
       sigmam(2)=0.16d0 	! from CMS
      ELSE
       sigmam(1)=0.08095d0	! from ATLAS
       sigmam(2)=0.12d0 	! from CMS
      ENDIF

*   Combination

      wmsum=0d0
      wpsum=0d0
      ksum=0d0
      DO I=1,2
       wp(i)=1d0/sigmap(i)**2
       wpsum=wpsum+wp(i)
       wm(i)=1d0/sigmam(i)**2
       wmsum=wmsum+wm(i)
      ENDDO
      ksum=(wm(1)+wp(1))/2d0*k(1)+(wm(2)+wp(2))*k(2)/2d0
      ksum=ksum/((wm(1)+wp(1))/2d0+(wp(2)+wm(2))/2d0)
      CMIN(3)=ksum-2d0/dsqrt(wmsum)
      CMAX(3)=ksum+2d0/dsqrt(wpsum)

************************************************************************
*   k_gluon

      IF(FLAG.EQ.1)THEN
       k(1)=0.94913d0		! from ATLAS
       k(2)=0.92d0		! from CMS
      ELSE
       k(1)=0.94092d0		! from ATLAS
       k(2)=0.93d0		! from CMS
      ENDIF

c   sigma+

      IF(FLAG.EQ.1)THEN
       sigmap(1)=0.07273d0	! from ATLAS
       sigmap(2)=0.08d0		! from CMS
      ELSE
       sigmap(1)=0.06919d0	! from ATLAS
       sigmap(2)=0.07d0		! from CMS
      ENDIF

c   sigma-

      IF(FLAG.EQ.1)THEN
       sigmam(1)=0.0668d0	! from ATLAS
       sigmam(2)=0.08d0 	! from CMS
      ELSE
       sigmam(1)=0.064d0	! from ATLAS
       sigmam(2)=0.07d0 	! from CMS
      ENDIF

*   Combination

      wmsum=0d0
      wpsum=0d0
      ksum=0d0
      DO I=1,2
       wp(i)=1d0/sigmap(i)**2
       wpsum=wpsum+wp(i)
       wm(i)=1d0/sigmam(i)**2
       wmsum=wmsum+wm(i)
      ENDDO
      ksum=(wm(1)+wp(1))/2d0*k(1)+(wm(2)+wp(2))*k(2)/2d0
      ksum=ksum/((wm(1)+wp(1))/2d0+(wp(2)+wm(2))/2d0)
      CMIN(4)=ksum-2d0/dsqrt(wmsum)
      CMAX(4)=ksum+2d0/dsqrt(wpsum)

************************************************************************
*   k_gamma

      IF(FLAG.EQ.1)THEN
       k(1)=1.0094d0		! from ATLAS
       k(2)=1.10d0		! from CMS
      ELSE
       k(1)=0.98367d0		! from ATLAS
       k(2)=1.07d0		! from CMS
      ENDIF

c   sigma+

      IF(FLAG.EQ.1)THEN
       sigmap(1)=0.06324d0	! from ATLAS
       sigmap(2)=0.08d0		! from CMS
      ELSE
       sigmap(1)=0.05281d0	! from ATLAS
       sigmap(2)=0.05d0		! from CMS
      ENDIF

c   sigma-

      IF(FLAG.EQ.1)THEN
       sigmam(1)=0.06103d0	! from ATLAS
       sigmam(2)=0.08d0 	! from CMS
      ELSE
       sigmam(1)=0.05011d0	! from ATLAS
       sigmam(2)=0.06d0 	! from CMS
      ENDIF

*   Combination

      wmsum=0d0
      wpsum=0d0
      ksum=0d0
      DO I=1,2
       wp(i)=1d0/sigmap(i)**2
       wpsum=wpsum+wp(i)
       wm(i)=1d0/sigmam(i)**2
       wmsum=wmsum+wm(i)
      ENDDO
      ksum=(wm(1)+wp(1))/2d0*k(1)+(wm(2)+wp(2))*k(2)/2d0
      ksum=ksum/((wm(1)+wp(1))/2d0+(wp(2)+wm(2))/2d0)
      CMIN(5)=ksum-2d0/dsqrt(wmsum)
      CMAX(5)=ksum+2d0/dsqrt(wpsum)

************************************************************************
*   k_tau

      IF(FLAG.EQ.1)THEN
       k(1)=0.92927d0		! from ATLAS
       k(2)=0.92d0		! from CMS
      ELSE
       k(1)=0.91166d0		! from ATLAS
       k(2)=0.91d0		! from CMS
      ENDIF

c   sigma+

      IF(FLAG.EQ.1)THEN
       sigmap(1)=0.07302d0	! from ATLAS
       sigmap(2)=0.08d0		! from CMS
      ELSE
       sigmap(1)=0.06569d0	! from ATLAS
       sigmap(2)=0.07d0		! from CMS
      ENDIF

c   sigma-

      IF(FLAG.EQ.1)THEN
       sigmam(1)=0.07012d0	! from ATLAS
       sigmam(2)=0.08d0 	! from CMS
      ELSE
       sigmam(1)=0.06296d0	! from ATLAS
       sigmam(2)=0.07d0 	! from CMS
      ENDIF

*   Combination

      wmsum=0d0
      wpsum=0d0
      ksum=0d0
      DO I=1,2
       wp(i)=1d0/sigmap(i)**2
       wpsum=wpsum+wp(i)
       wm(i)=1d0/sigmam(i)**2
       wmsum=wmsum+wm(i)
      ENDDO
      ksum=(wm(1)+wp(1))/2d0*k(1)+(wm(2)+wp(2))*k(2)/2d0
      ksum=ksum/((wm(1)+wp(1))/2d0+(wp(2)+wm(2))/2d0)
      CMIN(6)=ksum-2d0/dsqrt(wmsum)
      CMAX(6)=ksum+2d0/dsqrt(wpsum)

************************************************************************
*   BR_BSM

      CMIN(7)=0d0
      CMAX(7)=0.16d0


      END

      
