      SUBROUTINE LHCHIG(PAR,PROB)

*   Subroutine to check LHC Higgs constraints
*      PROB(45) =/= 0  excluded by t -> bH+ (LHC)

      IMPLICIT NONE

      INTEGER I,J

      DOUBLE PRECISION PAR(*),PROB(*),SIG(5,8),S,ggF13,SMEM
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJSM,BREESM,BRMMSM,BRLLSM,BRCCSM
      DOUBLE PRECISION BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM
      DOUBLE PRECISION brtopbw,brtopbh,brtopneutrstop(5,2),TBH
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH

      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BR_top2body/brtopbw,brtopbh,brtopneutrstop
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     .       HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     .       HCBRSUSY,HCWIDTH
      COMMON/LHCSIG/SIG
      COMMON/SMEM/SMEM


* Loop over H1, H2, H3

      DO I=1,3

       DO J=1,8
        SIG(I,J)=0d0
       ENDDO

       CALL HDECAY(SMASS(I),BRJJSM,BREESM,BRMMSM,BRLLSM,
     .      BRCCSM,BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   H -> tautau
* VBF/VH
       IF(BRLLSM.NE.0d0)SIG(I,1)=CV(I)**2*BRLL(I)/BRLLSM
* ggF
       IF(BRLLSM.NE.0d0)SIG(I,2)=CJ(I)**2*BRLL(I)/BRLLSM
       
*   H -> bb
* VBF/VH
       IF(BRBBSM.NE.0d0)SIG(I,3)=CV(I)**2*BRBB(I)/BRBBSM
* ttH
       IF(BRBBSM.NE.0d0)SIG(I,4)=CU(I)**2*BRBB(I)/BRBBSM

*   H -> ZZ/WW
* VBF/VH
       IF(BRZZSM.NE.0d0)SIG(I,5)=CV(I)**2*BRZZ(I)/BRZZSM
* ggF
       IF(BRZZSM.NE.0d0)SIG(I,6)=CJ(I)**2*BRZZ(I)/BRZZSM
       
*   H -> gammagamma
* VBF/VH
       IF(BRGGSM.NE.0d0)SIG(I,7)=CV(I)**2*BRGG(I)/BRGGSM
* ggF
       IF(BRGGSM.NE.0d0)SIG(I,8)=CJ(I)**2*BRGG(I)/BRGGSM

      ENDDO

* Loop over A1, A2

      DO I=1,2

       DO J=1,8
        SIG(I+3,J)=0d0
       ENDDO

       CALL HDECAY(PMASS(I),BRJJSM,BREESM,BRMMSM,BRLLSM,
     .      BRCCSM,BRBBSM,BRTTSM,BRWWSM,BRZZSM,BRGGSM,BRZGSM)

*   A -> tautau
* ggF
       IF(BRLLSM.NE.0d0)SIG(I+3,2)=CJ(I+3)**2*BRLL(I+3)/BRLLSM
       
*   A -> bb
* ttH
       IF(BRBBSM.NE.0d0)SIG(I+3,4)=CU(I+3)**2*BRBB(I+3)/BRBBSM

*   A -> ZZ/WW
* ggF
       SIG(I+3,6)=0d0
       
*   A -> gammagamma
* ggF
       IF(BRGGSM.NE.0d0)SIG(I+3,8)=CJ(I+3)**2*BRGG(I+3)/BRGGSM

      ENDDO

* Bound on Br(t->bH+)*BR(H+->tau nu)

      PROB(45)=DDIM(brtopbh*HCBRL/TBH(CMASS),1d0)

* Other LHC bounds

      CALL HIGGS_KAPPAS(PROB)
      CALL HEAVYHA_TAUTAU(PROB)
      CALL HEAVYHA_TOPTOP(PROB)
      CALL HEAVYH_VV(PROB)
      CALL LIGHTHA_GAMGAM(PROB)
      CALL HSM_AA_LEPTONS(PROB)
      CALL HSM_AA_PHOTONS(PROB)
      CALL HSM_ZA(PROB)
      CALL HEAVYA_ZHSM(PROB)
      CALL HEAVYHA_ZAH(PROB)
      CALL HEAVYHA_HSMHA(PROB)
      CALL HEAVYH_HSMHSM(PROB)

      END


      DOUBLE PRECISION FUNCTION TBH(M)

* ATLAS constraints on BR(t->bH+)*BR(H+->taunu), ATLAS-CONF-2011-151 tab.5

      IMPLICIT NONE
      INTEGER I,N
      PARAMETER(N=8)
      DOUBLE PRECISION X(N),Y(N),M

      DATA X/90d0,100d0,110d0,120d0,130d0,140d0,150d0,160d0/ 
      DATA Y/.104d0,.098d0,.095d0,.077d0,.066d0,.071d0,.052d0,.141d0/ 

      TBH=1d9
      DO I=1,N-1
       IF((M.GE.X(I)).AND.(M.LE.X(I+1)))THEN
        TBH=(Y(I)+(Y(I+1)-Y(I))*(M-X(I))/(X(I+1)-X(I)))
        RETURN
       ENDIF
      ENDDO

      END


      SUBROUTINE HIGGS_KAPPAS(PROB)

*      If HFLAG=/=0:
*      If HFLAG = 1: H1 is SM-like
*      If HFLAG = 2: H2 is SM-like
*      If SMASS(X) is in the window MHMIN...MHMAX: DX=0
*      If SMASS(X) < MHMIN: DX = SMASS(X)/MHMIN - 1 < 0
*      If SMASS(X) > MHMAX: DX = SMASS(X)/MHMAX - 1 > 0
*
*      CV=C(1), CU=C(2), CB=C(3), CJ=C(4) (2gluons), CG=C(5) (2photons),
*      CL=C(6) (2taus), BR_bsm=C(7) are constrained by the combination
*      of kappas from subroutine couplingfits
*
*      PROB(46) =/= 0  no Higgs in the MHmin-MHmax GeV range
*      PROB(67) =/= 0 k_WZ(H_SM) 2 sigma away from LHC measured value
*      PROB(68) =/= 0 k_top(H_SM) 2 sigma away from LHC measured value
*      PROB(69) =/= 0 k_bot(H_SM) 2 sigma away from LHC measured value
*      PROB(70) =/= 0 k_glu(H_SM) 2 sigma away from LHC measured value
*      PROB(71) =/= 0 k_gam(H_SM) 2 sigma away from LHC measured value
*      PROB(72) =/= 0 k_tau(H_SM) 2 sigma away from LHC measured value
*      PROB(73) =/= 0 BR_bsm(H_SM) 2 sigma away from LHC measured value

      IMPLICIT NONE
      INTEGER I,GMUFLAG,HFLAG,PDGH1,PDGH2
      DOUBLE PRECISION PROB(*),D1,D2,C(7),CMIN(7),CMAX(7)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION MHmin,MHmax
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5),CB(5),CL(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)

      COMMON/GMUFLAG/GMUFLAG,HFLAG
      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/PDGHIGGS/PDGH1,PDGH2
      COMMON/CMIMAX/CMIN,CMAX

      IF(HFLAG.EQ.2)THEN
       D1=1d99
      ELSE
       D1=DDIM(SMASS(1)/MHMAX,1d0)-DDIM(1d0,SMASS(1)/MHMIN)
      ENDIF
      IF(HFLAG.EQ.1)THEN
       D2=1d99
      ELSE
       D2=DDIM(SMASS(2)/MHMAX,1d0)-DDIM(1d0,SMASS(2)/MHMIN)
      ENDIF

      IF(D1.EQ.0d0 .AND. D2.EQ.0d0)THEN
       IF(CV(1)**2.GT.CV(2)**2)THEN
        PDGH1=25
        PDGH2=35
        C(1)=CV(1)
        C(2)=CU(1)
        C(3)=CB(1)
        C(4)=CJ(1)
        C(5)=CG(1)
        C(6)=CL(1)
        C(7)=BRNEU(1,1,1)+BRHAA(1,1)*BRNEU(4,1,1)**2
       ELSE
        PDGH2=25
        PDGH1=35
        C(1)=CV(2)
        C(2)=CU(2)
        C(3)=CB(2)
        C(4)=CJ(2)
        C(5)=CG(2)
        C(6)=CL(2)
        C(7)=BRNEU(2,1,1)+BRHAA(2,1)*BRNEU(4,1,1)**2
     .      +BRHHH(1)*BRNEU(1,1,1)**2
       ENDIF
      ELSEIF(D1.EQ.0d0)THEN
       PDGH1=25
       PDGH2=35
       C(1)=CV(1)
       C(2)=CU(1)
       C(3)=CB(1)
       C(4)=CJ(1)
       C(5)=CG(1)
       C(6)=CL(1)
       C(7)=BRNEU(1,1,1)+BRHAA(1,1)*BRNEU(4,1,1)**2
      ELSEIF(D2.EQ.0d0)THEN
       PDGH2=25
       PDGH1=35
       C(1)=CV(2)
       C(2)=CU(2)
       C(3)=CB(2)
       C(4)=CJ(2)
       C(5)=CG(2)
       C(6)=CL(2)
        C(7)=BRNEU(2,1,1)+BRHAA(2,1)*BRNEU(4,1,1)**2
     .      +BRHHH(1)*BRNEU(1,1,1)**2
      ELSE
       IF(DABS(D1).LT.DABS(D2))THEN
        PROB(46)=D1
        PDGH1=25
        PDGH2=35
        C(1)=CV(1)
        C(2)=CU(1)
        C(3)=CB(1)
        C(4)=CJ(1)
        C(5)=CG(1)
        C(6)=CL(1)
        C(7)=BRNEU(1,1,1)+BRHAA(1,1)*BRNEU(4,1,1)**2
       ELSE
        PROB(46)=D2
        PDGH2=25
        PDGH1=35
        C(1)=CV(2)
        C(2)=CU(2)
        C(3)=CB(2)
        C(4)=CJ(2)
        C(5)=CG(2)
        C(6)=CL(2)
        C(7)=BRNEU(2,1,1)+BRHAA(2,1)*BRNEU(4,1,1)**2
     .      +BRHHH(1)*BRNEU(1,1,1)**2
       ENDIF
      ENDIF

      PROB(67)=-DDIM(1d0,C(1)/CMIN(1))
      DO I=2,6
       PROB(66+I)=DDIM(C(I)/CMAX(I),1d0)-DDIM(1d0,C(I)/CMIN(I))
      ENDDO
      PROB(73)=DDIM(C(7)/CMAX(7),1d0)

      END


      SUBROUTINE HEAVYHA_TAUTAU(PROB)

* Constraints from ggF/bb->H/A->tautau
* PROB(51) =/= 0: excluded by H/A->tautau

      IMPLICIT NONE

      INTEGER I,I1,J,J1,JBAR,JBARbb,NX,NA13,NC13,NC22
      PARAMETER(NX=18,NA13=14,NC13=28,NC22=32)

      DOUBLE PRECISION PROB(*),HMAS(NX),LCMS(NX),LCMSbb(NX)
      DOUBLE PRECISION LATLAS(NX),LATLASbb(NX)
      DOUBLE PRECISION HCMS13(NC13),LCMS13(NC13),LCMSbb13(NC13)
      DOUBLE PRECISION HATLAS13(NA13),LATLAS13(NA13),LATLASbb13(NA13)
      DOUBLE PRECISION HCMS22(NC22),LCMS22(NC22),LCMSbb22(NC22)
      DOUBLE PRECISION MH(5),SIG(5),SIGbb(5),LCMSH(5),LATLASH(5)
      DOUBLE PRECISION LCMSHbb(5),LATLASHbb(5)
      DOUBLE PRECISION DEL,SIGTOT,MBAR,LCMSMB,LATLASMB
      DOUBLE PRECISION SIGTOTbb,MBARbb,LCMSMBbb,LATLASMBbb
      DOUBLE PRECISION ggF8,ggF13,bbH8,bbH13
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION MHmin,MHmax,D1,D2,CBSM(5)

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/HIGGSFIT/MHmin,MHmax

* 8 TeV limits
      DATA HMAS/90d0,100d0,120d0,140d0,160d0,180d0,200d0,250d0,300d0,
     . 350d0,400d0,450d0,500d0,600d0,700d0,800d0,900d0,1000d0/

* Upper limit on ggF->H->tautau (8 TeV) from CMS-PAS-HIG-13-021, Tab. 7
      DATA LCMS/50.2d0,31.3d0,7.38d0,2.27d0,.845d0,.549d0,.517d0,.315d0,
     . .15d0,.112d0,.103d0,.607d-1,.385d-1,.193d-1,.143d-1,.115d-1,
     . .923d-2,.865d-2/

* Upper limit on bbH->tautau (8 TeV) from CMS-PAS-HIG-13-021, Tab. 8
      DATA LCMSbb/6.03d0,4.14d0,1.76d0,1.25d0,.814d0,.659d0,.553d0,
     . .217d0,
     . .975d-1,.638d-1,.613d-1,.431d-1,.320d-1,.203d-1,.173d-1,.166d-1,
     . .146d-1,.133d-1/

* Upper limit on ggF->H->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLAS/29.1d0,24.0d0,5.25d0,2.02d0,1.39d0,1.00d0,.794d0,
     .  .281d0,.127d0,.112d0,.773d-1,.400d-1,.240d-1,.177d-1,.127d-1,
     .  .993d-2,.840d-2,.735d-2/

* Upper limit on bbH->tautau (8 TeV) from ATLAS-CONF-2014-049, Fig. 7
      DATA LATLASbb/6.32d0,6.32d0,2.73d0,1.27d0,.966d0,.606d0,.393d0,
     .  .305d0,.116d0,.101d0,.656d-1,.363d-1,.238d-1,.159d-1,.117d-1,
     .  .943d-2,.785d-2,.716d-2/


* 13 TeV limits
      DATA HATLAS13/2d2,2.5d2,3d2,3.5d2,4d2,5d2,6d2,7d2,8d2,
     . 1d3,1.2d3,1.5d3,2d3,2.25d3/
      DATA HCMS13/9d1,1d2,1.1d2,1.2d2,1.3d2,1.4d2,1.6d2,1.8d2,2d2,
     . 2.5d2,3.5d2,4d2,4.5d2,5d2,6d2,7d2,8d2,9d2,1d3,1.2d3,1.4d3,
     . 1.6d3,1.8d3,2d3,2.3d3,2.6d3,2.9d3,3.2d3/

* Upper limit on ggF->H->tautau (13 TeV) from ATLAS 2002.12223, Fig. 2(a)
      DATA LATLAS13/2.28d-1,1.68d-1,1.65d-1,1.14d-1,8.29d-2,4.95d-2,
     . 1.48d-2,5.07d-3,2.93d-3,1.79d-3,1.58d-3,1.54d-3,1.27d-3,1.16d-3/

* Upper limit on bbH->tautau (13 TeV) from ATLAS 2002.12223, Fig. 2(b)
      DATA LATLASbb13/2.03d-1,1.07d-1,1.18d-1,1.19d-1,7.72d-2,2.85d-2,
     . 7.34d-3,2.45d-3,1.64d-3,1.11d-3,9.58d-4,1.01d-3,1.03d-3,1.02d-3/

* Upper limit on ggF->H->tautau (13 TeV) from CMS 1803.06553, Fig. 7(a)
      DATA LCMS13/1.74d1,2.59d1,1.42d1,8.61d0,4.73d0,2.82d0,1.86d0,
     . 1.30d0,8.29d-1,5.56d-1,1.09d-1,7.43d-2,6.25d-2,5.91d-2,5.01d-2,
     . 3.51d-2,2.47d-2,1.73d-2,1.37d-2,8.93d-3,6.42d-3,4.99d-3,4.43d-3,
     . 4.00d-3,3.74d-3,3.64d-3,3.53d-3,3.40d-3/

* Upper limit on bbH->tautau (13 TeV) from CMS 1803.06553, Fig. 7(b)
      DATA LCMSbb13/1.39d1,7.13d0,7.13d0,5.45d0,3.15d0,2.71d0,1.45d0,
     . 1.05d0,7.86d-1,3.61d-1,5.82d-2,4.38d-2,3.59d-2,3.00d-2,1.44d-2,
     . 8.86d-3,6.03d-3,4.74d-3,3.70d-3,2.90d-3,2.84d-3,2.46d-3,2.41d-3,
     . 2.40d-3,2.16d-3,2.35d-3,2.45d-3,2.48d-3/


* CMS 2208.02717 (13 TeV) limits
      DATA HCMS22/60d0,80d0,95d0,100d0,120d0,125d0,130d0,140d0,
     .160d0,180d0,200d0,250d0,300d0,350d0,400d0,450d0,500d0,600d0,
     .700d0,800d0,900d0,1000d0,1200d0,1400d0,1600d0,1800d0,2000d0,
     .2300d0,2600d0,2900d0,3200d0,3500d0/

* Upper limit on ggF->H->tautau from Fig. 10(a)
      DATA LCMS22/0.4731d1,0.3357d1,0.1528d2,0.1031d2,0.1829d1,
     .0.1094d1,0.7822d0,0.4895d0,0.3523d0,0.2251d0,0.2902d0,0.9341d-1,
     .0.5349d-1,0.3550d-1,0.3155d-1,0.2977d-1,0.2044d-1,0.1031d-1,
     .0.6764d-2,0.5813d-2,0.6194d-2,0.6469d-2,0.4737d-2,0.3016d-2,
     .0.1745d-2,0.1125d-2,0.7931d-3,0.5724d-3,0.4280d-3,0.3527d-3,
     .0.3050d-3,0.2760d-3/

* Upper limit on bbH->tautau from Fig. 10(b)
      DATA LCMSbb22/0.2435d2,0.1128d2,0.3126d1,0.2212d1,0.6480d0,
     .0.5643d0,0.4237d0,0.3337d0,0.2519d0,0.2018d0,0.2274d0,0.7975d-1,
     .0.4472d-1,0.3544d-1,0.3162d-1,0.2992d-1,0.2242d-1,0.1154d-1,
     .0.6579d-2,0.4048d-2,0.3097d-2,0.2313d-2,0.1462d-2,0.1058d-2,
     .0.9034d-3,0.7669d-3,0.6303d-3,0.5352d-3,0.4137d-3,0.3576d-3,
     .0.3442d-3,0.3123d-3/


      PROB(51)=0d0

      DO I=1,5
       CBSM(I)=1d0
      ENDDO
      D1=DDIM(SMASS(1)/MHMAX,1d0)-DDIM(1d0,SMASS(1)/MHMIN)
      D2=DDIM(SMASS(2)/MHMAX,1d0)-DDIM(1d0,SMASS(2)/MHMIN)
      IF(D1.EQ.0d0 .AND. D2.EQ.0d0)THEN
       IF(CV(1)**2.GT.CV(2)**2)THEN
        CBSM(1)=0d0
       ELSE
        CBSM(2)=0d0
       ENDIF
      ELSEIF(D1.EQ.0d0)THEN
        CBSM(1)=0d0
      ELSEIF(D2.EQ.0d0)THEN
        CBSM(2)=0d0
      ENDIF

* Loop over 5 Higgses
      DO I=1,5
        IF(I.LE.3) MH(I)=SMASS(I)
        IF(I.GE.4) MH(I)=PMASS(I-3)


* 8 TeV limits
        J=1
        DOWHILE(HMAS(J).LE.MH(I) .AND. J.LT.NX)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LE.HMAS(NX)) THEN

* ggF Signal cross section*BR:
          SIG(I)=CBSM(I)*CJ(I)**2*BRLL(I)*ggF8(MH(I))

* CMS ggF limit:
          LCMSH(I)=LCMS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMS(J)-LCMS(J-1))
          PROB(51)=PROB(51)+DDIM(SIG(I)/LCMSH(I),1d0)

* ATLAS ggF limit:
          LATLASH(I)=LATLAS(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLAS(J)-LATLAS(J-1))
* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLAS(J-1)=LATLAS(7)=.96d0 and not .794d0:
          IF(J.EQ.8) THEN
            LATLASH(I)=.96d0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLAS(J)-.96d0)
          ENDIF
          PROB(51)=PROB(51)+DDIM(SIG(I)/LATLASH(I),1d0)

* bbH Signal cross section*BR:
          SIGbb(I)=CBSM(I)*CB(I)**2*BRLL(I)*bbH8(MH(I))

* CMS Hbb limit:
          LCMSHbb(I)=LCMSbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LCMSbb(J)-LCMSbb(J-1))
          PROB(51)=PROB(51)+DDIM(SIGbb(I)/LCMSHbb(I),1d0)

* ATLAS Hbb limit:
          LATLASHbb(I)=LATLASbb(J-1)+(MH(I)-HMAS(J-1))/
     .      (HMAS(J)-HMAS(J-1))*(LATLASbb(J)-LATLASbb(J-1))
* Correct for jump in Fig.7 at MA=200 GeV: J=8, 
* modif. LATLASbb(J-1)=LATLASbb(7)=.858d0 and not .393d0:
          IF(J.EQ.8) THEN
            LATLASHbb(I)=.858d0+(MH(I)-HMAS(J-1))/
     .        (HMAS(J)-HMAS(J-1))*(LATLASbb(J)-.858d0)
          ENDIF
          PROB(51)=PROB(51)+DDIM(SIGbb(I)/LATLASHbb(I),1d0)

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HMAS(J1).LE.MH(I1) .AND. J1.LT.NX)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HMAS(1) .AND. MH(I1).LE.HMAS(NX)) THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG(I)*MH(I)+SIG(I1)*MH(I1))/(SIG(I)+SIG(I1))
             JBAR=1
             DOWHILE(HMAS(JBAR).LE.MBAR.AND.JBAR.LT.NX)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HMAS(JBARbb).LE.MBARbb.AND.JBARbb.LT.NX)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG(I)+SIG(I1)
     .         -SIG(I)*SIG(I1)*DEL/(SIG(I)+SIG(I1))
* Continue only if SIGTOT > 0 and 90<MBAR<1000
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0d0.AND.MBAR.GE.HMAS(1).AND.
     .         MBAR.LE.HMAS(NX).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* CMS ggF limit at MBAR:
               LCMSMB=LCMS(JBAR-1)+(MBAR-HMAS(JBAR-1))/
     .          (HMAS(JBAR)-HMAS(JBAR-1))*(LCMS(JBAR)-LCMS(JBAR-1))
               PROB(51)=PROB(51)+DDIM(SIGTOT/LCMSMB,1d0)

* ATLAS ggF limit at MBAR:
               LATLASMB=LATLAS(JBAR-1)+(MBAR-HMAS(JBAR-1))/
     .          (HMAS(JBAR)-HMAS(JBAR-1))*(LATLAS(JBAR)-LATLAS(JBAR-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBAR=8, 
* modif. LATLAS(JBAR-1)=LATLAS(7)=.96d0 and not .794d0:
               IF(J.EQ.8) THEN
                 LATLASMB=.96d0+(MBAR-HMAS(JBAR-1))/
     .            (HMAS(JBAR)-HMAS(JBAR-1))*(LATLAS(JBAR)-.96d0)
               ENDIF
               PROB(51)=PROB(51)+DDIM(SIGTOT/LATLASMB,1d0)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 90<MBARbb<1000
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0d0.AND.MBARbb.GE.HMAS(1).AND.
     .         MBARbb.LE.HMAS(NX).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* CMS bbH limit at MBARbb:
              LCMSMBbb=LCMSbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .         (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .         (LCMSbb(JBARbb)-LCMSbb(JBARbb-1))
              PROB(51)=PROB(51)+DDIM(SIGTOTbb/LCMSMBbb,1d0)

* ATLAS bbH limit at MBARbb:
              LATLASMBbb=LATLASbb(JBARbb-1)+(MBARbb-HMAS(JBARbb-1))/
     .         (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .         (LATLASbb(JBARbb)-LATLASbb(JBARbb-1))
* Correct for jump in Fig.7 at MA=200 GeV: JBARbb=8, 
* modif. LATLASbb(JBARbb-1)=LATLASbb(7)=.858d0 and not .393d0:
              IF(J.EQ.8) THEN
                LATLASMBbb=.858d0+(MBARbb-HMAS(JBARbb-1))/
     .           (HMAS(JBARbb)-HMAS(JBARbb-1))*
     .           (LATLASbb(JBARbb)-.858d0)
              ENDIF
              PROB(51)=PROB(51)+DDIM(SIGTOTbb/LATLASMBbb,1d0)
             ENDIF
            ENDIF
          ENDDO
        ENDIF



* ATLAS 13 TeV limits
        J=1
        DOWHILE(HATLAS13(J).LE.MH(I) .AND. J.LT.NA13)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LE.HATLAS13(NA13)) THEN

* ggF Signal cross section*BR:
          SIG(I)=CBSM(I)*CJ(I)**2*BRLL(I)*ggF13(MH(I))

* ATLAS ggF limit:
          LATLASH(I)=LATLAS13(J-1)+(MH(I)-HATLAS13(J-1))/
     .      (HATLAS13(J)-HATLAS13(J-1))*(LATLAS13(J)-LATLAS13(J-1))
          PROB(51)=PROB(51)+DDIM(SIG(I)/LATLASH(I),1d0)

* bbH Signal cross section*BR:
          SIGbb(I)=CBSM(I)*CB(I)**2*BRLL(I)*bbH13(MH(I))

* ATLAS Hbb limit:
          LATLASHbb(I)=LATLASbb13(J-1)+(MH(I)-HATLAS13(J-1))/
     .      (HATLAS13(J)-HATLAS13(J-1))*(LATLASbb13(J)-LATLASbb13(J-1))
          PROB(51)=PROB(51)+DDIM(SIGbb(I)/LATLASHbb(I),1d0)

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HATLAS13(J1).LE.MH(I1) .AND. J1.LT.NA13)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HATLAS13(1) .AND. MH(I1).LE.HATLAS13(NA13))THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG(I)*MH(I)+SIG(I1)*MH(I1))/(SIG(I)+SIG(I1))
             JBAR=1
             DOWHILE(HATLAS13(JBAR).LE.MBAR .AND. JBAR.LT.NA13)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HATLAS13(JBARbb).LE.MBARbb .AND. JBARbb.LT.NA13)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG(I)+SIG(I1)
     .         -SIG(I)*SIG(I1)*DEL/(SIG(I)+SIG(I1))
* Continue only if SIGTOT > 0 and 200<MBAR<2250
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0d0.AND.MBAR.GE.HATLAS13(1).AND.
     .         MBAR.LE.HATLAS13(NA13).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* ATLAS ggF limit at MBAR:
               LATLASMB=LATLAS13(JBAR-1)+(MBAR-HATLAS13(JBAR-1))/
     .          (HATLAS13(JBAR)-HATLAS13(JBAR-1))*
     .          (LATLAS13(JBAR)-LATLAS13(JBAR-1))
               PROB(51)=PROB(51)+DDIM(SIGTOT/LATLASMB,1d0)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 200<MBARbb<2250
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0d0.AND.MBARbb.GE.HATLAS13(1).AND.
     .         MBARbb.LE.HATLAS13(NA13).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* ATLAS bbH limit at MBARbb:
               LATLASMBbb=LATLASbb13(JBARbb-1)+
     .           (MBARbb-HATLAS13(JBARbb-1))/
     .           (HATLAS13(JBARbb)-HATLAS13(JBARbb-1))*
     .           (LATLASbb13(JBARbb)-LATLASbb13(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(SIGTOTbb/LATLASMBbb,1d0)
             ENDIF
            ENDIF
          ENDDO
        ENDIF


* CMS 13 TeV limits
        J=1
        DOWHILE(HCMS13(J).LE.MH(I) .AND. J.LT.NC13)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LE.HCMS13(NC13)) THEN

* ggF Signal cross section*BR:
          SIG(I)=CBSM(I)*CJ(I)**2*BRLL(I)*ggF13(MH(I))

* CMS ggF limit:
          LCMSH(I)=LCMS13(J-1)+(MH(I)-HCMS13(J-1))/
     .      (HCMS13(J)-HCMS13(J-1))*(LCMS13(J)-LCMS13(J-1))
          PROB(51)=PROB(51)+DDIM(SIG(I)/LCMSH(I),1d0)

* bbH Signal cross section*BR:
          SIGbb(I)=CBSM(I)*CB(I)**2*BRLL(I)*bbH13(MH(I))

* CMS Hbb limit:
          LCMSHbb(I)=LCMSbb13(J-1)+(MH(I)-HCMS13(J-1))/
     .      (HCMS13(J)-HCMS13(J-1))*(LCMSbb13(J)-LCMSbb13(J-1))
          PROB(51)=PROB(51)+DDIM(SIGbb(I)/LCMSHbb(I),1d0)

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HCMS13(J1).LE.MH(I1) .AND. J1.LT.NC13)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HCMS13(1) .AND. MH(I1).LE.HCMS13(NC13))THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG(I)*MH(I)+SIG(I1)*MH(I1))/(SIG(I)+SIG(I1))
             JBAR=1
             DOWHILE(HCMS13(JBAR).LE.MBAR .AND. JBAR.LT.NC13)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HCMS13(JBARbb).LE.MBARbb .AND. JBARbb.LT.NC13)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG(I)+SIG(I1)
     .         -SIG(I)*SIG(I1)*DEL/(SIG(I)+SIG(I1))
* Continue only if SIGTOT > 0 and 90<MBAR<3200
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0d0.AND.MBAR.GE.HCMS13(1).AND.
     .         MBAR.LE.HCMS13(NC13).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* CMS ggF limit at MBAR:
               LCMSMB=LCMS13(JBAR-1)+(MBAR-HCMS13(JBAR-1))/
     .          (HCMS13(JBAR)-HCMS13(JBAR-1))*
     .          (LCMS13(JBAR)-LCMS13(JBAR-1))
               PROB(51)=PROB(51)+DDIM(SIGTOT/LCMSMB,1d0)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 90<MBARbb<3200
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0d0.AND.MBARbb.GE.HCMS13(1).AND.
     .         MBARbb.LE.HCMS13(NC13).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* CMS bbH limit at MBARbb:
               LCMSMBbb=LCMSbb13(JBARbb-1)+
     .           (MBARbb-HCMS13(JBARbb-1))/
     .           (HCMS13(JBARbb)-HCMS13(JBARbb-1))*
     .           (LCMSbb13(JBARbb)-LCMSbb13(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(SIGTOTbb/LCMSMBbb,1d0)
             ENDIF
            ENDIF
          ENDDO
        ENDIF


* CMS 2208.02717 (13 TeV) limits
        J=1
        DOWHILE(HCMS22(J).LE.MH(I) .AND. J.LT.NC22)
          J=J+1
        ENDDO
        IF(J.GE.2 .AND. MH(I).LE.HCMS22(NC22)) THEN

* ggF Signal cross section*BR:
          SIG(I)=CBSM(I)*CJ(I)**2*BRLL(I)*ggF13(MH(I))

* CMS ggF limit:
          LCMSH(I)=LCMS22(J-1)+(MH(I)-HCMS22(J-1))/
     .      (HCMS22(J)-HCMS22(J-1))*(LCMS22(J)-LCMS22(J-1))
          PROB(51)=PROB(51)+DDIM(SIG(I)/LCMSH(I),1d0)

* bbH Signal cross section*BR:
          SIGbb(I)=CBSM(I)*CB(I)**2*BRLL(I)*bbH13(MH(I))

* CMS Hbb limit:
          LCMSHbb(I)=LCMSbb22(J-1)+(MH(I)-HCMS22(J-1))/
     .      (HCMS22(J)-HCMS22(J-1))*(LCMSbb22(J)-LCMSbb22(J-1))
          PROB(51)=PROB(51)+DDIM(SIGbb(I)/LCMSHbb(I),1d0)

* Combine signal rates of 2 Higgses
          DO I1=1,I-1
            J1=1
            DOWHILE(HCMS22(J1).LE.MH(I1) .AND. J1.LT.NC22)
              J1=J1+1
            ENDDO
            IF(MH(I1).GE.HCMS22(1) .AND. MH(I1).LE.HCMS22(NC22))THEN

* Average masses weighted by the signal rates (MBAR for ggF, MBARbb for bbH):
             MBAR=(SIG(I)*MH(I)+SIG(I1)*MH(I1))/(SIG(I)+SIG(I1))
             JBAR=1
             DOWHILE(HCMS22(JBAR).LE.MBAR .AND. JBAR.LT.NC22)
               JBAR=JBAR+1
             ENDDO
             MBARbb=(SIGbb(I)*MH(I)
     .         +SIGbb(I1)*MH(I1))/(SIGbb(I)+SIGbb(I1))
             JBARbb=1
             DOWHILE(HCMS22(JBARbb).LE.MBARbb .AND. JBARbb.LT.NC22)
               JBARbb=JBARbb+1
             ENDDO
* DEL=mass difference divided by a (small) resolution squared:
* [DEL < 1 only if |MH(I)-MH(I1)| < (MH(I)+MH(I1))/15d0;
*  otherwise the combined signal rate is small]
             DEL=((MH(I)-MH(I1))/(MH(I)+MH(I1))*15d0)**2

* Estimate of the combined ggF signal rates:
             SIGTOT=SIG(I)+SIG(I1)
     .         -SIG(I)*SIG(I1)*DEL/(SIG(I)+SIG(I1))
* Continue only if SIGTOT > 0 and 60<MBAR<3500
*      and |MH(I)-MH(I1)|/MBAR<0.20:
             IF(SIGTOT.GT.0d0.AND.MBAR.GE.HCMS22(1).AND.
     .         MBAR.LE.HCMS22(NC22).AND.
     .         dabs(MH(I)-MH(I1))/MBAR.LE.0.20d0) THEN

* CMS ggF limit at MBAR:
               LCMSMB=LCMS22(JBAR-1)+(MBAR-HCMS22(JBAR-1))/
     .          (HCMS22(JBAR)-HCMS22(JBAR-1))*
     .          (LCMS22(JBAR)-LCMS22(JBAR-1))
               PROB(51)=PROB(51)+DDIM(SIGTOT/LCMSMB,1d0)
             ENDIF

* Estimate of the combined bbH signal rates:
             SIGTOTbb=SIGbb(I)+SIGbb(I1)
     .         -SIGbb(I)*SIGbb(I1)*DEL/(SIGbb(I)+SIGbb(I1))
* Continue only if SIGTOTbb > 0 and 90<MBARbb<3200
*      and |MH(I)-MH(I1)|/MBARbb<0.20:
             IF(SIGTOTbb.GT.0d0.AND.MBARbb.GE.HCMS22(1).AND.
     .         MBARbb.LE.HCMS22(NC22).AND.
     .         dabs(MH(I)-MH(I1))/MBARbb.LE.0.20d0) THEN

* CMS bbH limit at MBARbb:
               LCMSMBbb=LCMSbb22(JBARbb-1)+
     .           (MBARbb-HCMS22(JBARbb-1))/
     .           (HCMS22(JBARbb)-HCMS22(JBARbb-1))*
     .           (LCMSbb22(JBARbb)-LCMSbb22(JBARbb-1))
               PROB(51)=PROB(51)+DDIM(SIGTOTbb/LCMSMBbb,1d0)
             ENDIF
            ENDIF
          ENDDO
        ENDIF


      ENDDO

      END


      SUBROUTINE HEAVYHA_TOPTOP(PROB)

* Constraints from H/A->toptop from 1908.01115 (CMS)
* PROB(77) =/= 0: excluded by H/A->toptop

      IMPLICIT NONE

      INTEGER I,J,N1,N2
      PARAMETER(N1=15,N2=18)

      DOUBLE PRECISION PROB(*),HM1(N1),HM2(N2),L,L1,L2,MH,RW
      DOUBLE PRECISION LH1(N1),LH2(N1),LH3(N1),LH4(N1),LH5(N2),LH6(N2)
      DOUBLE PRECISION LA1(N1),LA2(N1),LA3(N1),LA4(N1),LA5(N2),LA6(N2)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH

      DATA HM1/400d0,425d0,450d0,475d0,500d0,525d0,550d0,575d0,600d0,
     .625d0,650d0,675d0,700d0,725d0,750d0/
      DATA HM2/400d0,425d0,450d0,475d0,500d0,525d0,550d0,575d0,600d0,
     .625d0,650d0,675d0,700d0,725d0,735d0,740d0,745d0,750d0/
      DATA LH1/0.845d0,0.715d0,0.695d0,0.655d0,0.655d0,0.705d0,0.715d0,
     .0.705d0,0.635d0,0.615d0,0.625d0,0.705d0,0.785d0,0.875d0,0.935d0/
      DATA LH2/1.055d0,0.885d0,0.855d0,0.795d0,0.785d0,0.845d0,0.855d0,
     .0.865d0,0.765d0,0.735d0,0.735d0,0.825d0,0.935d0,1.045d0,1.115d0/
      DATA LH3/1.375d0,1.155d0,1.085d0,0.985d0,0.965d0,1.045d0,1.075d0,
     .1.125d0,0.965d0,0.905d0,0.865d0,0.955d0,1.095d0,1.305d0,1.425d0/
      DATA LH4/1.655d0,1.365d0,1.185d0,1.165d0,1.105d0,1.185d0,1.255d0,
     .1.275d0,1.285d0,1.225d0,1.135d0,1.135d0,1.345d0,1.625d0,1.755d0/
      DATA LH5/2.135d0,1.655d0,1.325d0,1.185d0,1.135d0,1.195d0,1.295d0,
     .1.335d0,1.285d0,1.225d0,1.115d0,1.115d0,1.045d0,1.055d0,0.995d0,
     .0.985d0,0.975d0,1.005d0/
      DATA LH6/2.585d0,2.085d0,1.665d0,1.655d0,1.545d0,1.595d0,1.655d0,
     .1.585d0,1.495d0,1.405d0,1.415d0,1.455d0,1.425d0,1.345d0,1.295d0,
     .1.285d0,1.245d0,1.265d0/
      DATA LA1/0.595d0,0.565d0,0.575d0,0.565d0,0.545d0,0.535d0,0.585d0,
     .0.585d0,0.585d0,0.555d0,0.575d0,0.605d0,0.685d0,0.775d0,0.825d0/
      DATA LA2/0.735d0,0.695d0,0.695d0,0.695d0,0.655d0,0.635d0,0.725d0,
     .0.715d0,0.725d0,0.675d0,0.705d0,0.725d0,0.815d0,0.935d0,1.015d0/
      DATA LA3/0.925d0,0.865d0,0.865d0,0.865d0,0.815d0,0.765d0,0.915d0,
     .0.935d0,0.955d0,0.885d0,0.885d0,0.875d0,0.995d0,1.205d0,1.325d0/
      DATA LA4/1.095d0,1.025d0,1.035d0,1.005d0,0.875d0,0.975d0,1.115d0,
     .1.175d0,1.215d0,1.085d0,1.055d0,1.045d0,1.245d0,1.545d0,1.615d0/
      DATA LA5/1.325d0,1.195d0,1.135d0,1.055d0,0.975d0,0.995d0,1.175d0,
     .1.325d0,1.385d0,1.325d0,1.355d0,1.515d0,1.055d0,0.985d0,0.995d0,
     .0.975d0,0.955d0,1.025d0/
      DATA LA6/1.825d0,1.785d0,1.625d0,1.675d0,1.525d0,1.495d0,1.605d0,
     .1.685d0,1.745d0,1.675d0,1.535d0,1.585d0,1.645d0,1.545d0,1.505d0,
     .1.435d0,1.475d0,1.335d0/

      PROB(77)=0d0

      DO J=1,3
       MH=SMASS(J)
       RW=WIDTH(J)/MH

       IF(MH.GE.HM1(1).AND.MH.LT.HM1(N1))THEN

        IF(RW.LE.0.5d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L=LH1(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LH1(I)-LH1(I-1))
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.1d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L1=LH1(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LH1(I)-LH1(I-1))
         L2=LH2(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LH2(I)-LH2(I-1))
         L=L1+(RW-0.5d-2)/0.5d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.2.5d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L1=LH2(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LH2(I)-LH2(I-1))
         L2=LH3(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LH3(I)-LH3(I-1))
         L=L1+(RW-1d-2)/1.5d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.5d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L1=LH3(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LH3(I)-LH3(I-1))
         L2=LH4(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LH4(I)-LH4(I-1))
         L=L1+(RW-2.5d-2)/2.5d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.10d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L1=LH4(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LH4(I)-LH4(I-1))
         I=1
         DOWHILE(HM2(I).LE.MH)
          I=I+1
         ENDDO
         L2=LH5(I-1)+(MH-HM2(I-1))/(HM2(I)-HM2(I-1))*(LH5(I)-LH5(I-1))
         L=L1+(RW-5d-2)/5d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.25d-2)THEN

         I=1
         DOWHILE(HM2(I).LE.MH)
          I=I+1
         ENDDO
         L1=LH5(I-1)+(MH-HM2(I-1))/(HM2(I)-HM2(I-1))*(LH5(I)-LH5(I-1))
         L2=LH6(I-1)+(MH-HM2(I-1))/(HM2(I)-HM2(I-1))*(LH6(I)-LH6(I-1))
         L=L1+(RW-10d-2)/15d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSE

         I=1
         DOWHILE(HM2(I).LE.MH)
          I=I+1
         ENDDO
         L=LH6(I-1)+(MH-HM2(I-1))/(HM2(I)-HM2(I-1))*(LH6(I)-LH6(I-1))
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ENDIF

       ENDIF
      ENDDO

      DO J=4,5
       MH=PMASS(J-3)
       RW=WIDTH(J)/MH

       IF(MH.GE.HM1(1).AND.MH.LT.HM1(N1))THEN

        IF(RW.LE.0.5d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L=LA1(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LA1(I)-LA1(I-1))
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.1d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L1=LA1(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LA1(I)-LA1(I-1))
         L2=LA2(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LA2(I)-LA2(I-1))
         L=L1+(RW-0.5d-2)/0.5d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.2.5d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L1=LA2(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LA2(I)-LA2(I-1))
         L2=LA3(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LA3(I)-LA3(I-1))
         L=L1+(RW-1d-2)/1.5d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.5d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L1=LA3(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LA3(I)-LA3(I-1))
         L2=LA4(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LA4(I)-LA4(I-1))
         L=L1+(RW-2.5d-2)/2.5d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.10d-2)THEN

         I=1
         DOWHILE(HM1(I).LE.MH)
          I=I+1
         ENDDO
         L1=LA4(I-1)+(MH-HM1(I-1))/(HM1(I)-HM1(I-1))*(LA4(I)-LA4(I-1))
         I=1
         DOWHILE(HM2(I).LE.MH)
          I=I+1
         ENDDO
         L2=LA5(I-1)+(MH-HM2(I-1))/(HM2(I)-HM2(I-1))*(LA5(I)-LA5(I-1))
         L=L1+(RW-5d-2)/5d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSEIF(RW.LE.25d-2)THEN

         I=1
         DOWHILE(HM2(I).LE.MH)
          I=I+1
         ENDDO
         L1=LA5(I-1)+(MH-HM2(I-1))/(HM2(I)-HM2(I-1))*(LA5(I)-LA5(I-1))
         L2=LA6(I-1)+(MH-HM2(I-1))/(HM2(I)-HM2(I-1))*(LA6(I)-LA6(I-1))
         L=L1+(RW-10d-2)/15d-2*(L2-L1)
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ELSE

         I=1
         DOWHILE(HM2(I).LE.MH)
          I=I+1
         ENDDO
         L=LA6(I-1)+(MH-HM2(I-1))/(HM2(I)-HM2(I-1))*(LA6(I)-LA6(I-1))
         PROB(77)=PROB(77)+DDIM(DABS(CU(J))/L,1d0)

        ENDIF

       ENDIF
      ENDDO

      END


      SUBROUTINE HEAVYH_VV(PROB)

* Constraints from H->WW from ATLAS: 1710.01123 + CMS: 2109.06055
* Constraints from H->ZZ from ATLAS: 1712.06386 + CMS: 2109.08268
* Constraints from H->VV from ATLAS: 1808.02380 + CMS: 2210.00043
* PROB(50) =/= 0: excluded by H->VV

      IMPLICIT NONE

      INTEGER I,J,N1G,N1V,N2G,N2V,N3G,N3V,N4,N5,N6
      PARAMETER (N1G=26,N1V=21,N2G=57,N2V=54,N3G=18,N3V=13,
     .           N4=36,N5=11,N6=48)

      DOUBLE PRECISION PROB(*)
      DOUBLE PRECISION M1G(N1G),L1G(N1G),M1V(N1V),L1V(N1V)
      DOUBLE PRECISION M2G(N2G),L2G(N2G),M2V(N2V),L2V(N2V)
      DOUBLE PRECISION M3G(N3G),L3G(N3G),M3V(N3V),L3V(N3V)
      DOUBLE PRECISION M4(N4),L4G(N4),L4V(N4)
      DOUBLE PRECISION M5(N5),L5G(N5),L5V(N5)
      DOUBLE PRECISION M6(N6),L6(N6)
      DOUBLE PRECISION MH,SIG,L,ggF13,VBF13
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH

*   ATLAS 1710.01123 H->WW
*   fig. 5a ggF
      DATA M1G/200d0,250d0,300d0,400d0,500d0,600d0,700d0,750d0,800d0,
     .900d0,1000d0,1200d0,1400d0,1600d0,1800d0,2000d0,2200d0,2400d0,
     .2600d0,2800d0,3000d0,3200d0,3400d0,3600d0,3800d0,4000d0/
      DATA L1G/6.388d0,4.632d0,2.930d0,1.280d0,0.6347d0,0.2939d0,
     .0.2216d0,0.2152d0,0.2050d0,0.1607d0,0.1361d0,0.1007d0,0.7298d-1,
     .0.5040d-1,0.3283d-1,0.2289d-1,0.1812d-1,0.1434d-1,0.1239d-1,
     .0.1113d-1,0.9903d-2,0.9160d-2,0.8640d-2,0.7915d-2,0.7687d-2,
     .0.7612d-2/
*   fig. 5b VBF
      DATA M1V/200d0,250d0,300d0,400d0,500d0,600d0,700d0,750d0,800d0,
     .900d0,1000d0,1200d0,1400d0,1600d0,1800d0,2000d0,2200d0,2400d0,
     .2600d0,2800d0,3000d0/
      DATA L1V/1.340d0,0.9523d0,0.6835d0,0.4321d0,0.27054d0,0.14637d0,
     .0.8895d-1,0.7038d-1,0.6446d-1,0.4581d-1,0.3486d-1,0.2291d-1,
     .0.1612d-1,0.1227d-1,0.9806d-2,0.8723d-2,0.7912d-2,0.7318d-2,
     .0.6703d-2,0.6509d-2,0.6260d-2/

*   ATLAS 1712.06386 H->ZZ
*   fig. 6a ggF
      DATA M2G/200d0,205d0,215d0,220d0,230d0,240d0,250d0,260d0,270d0,
     .280d0,290d0,300d0,320d0,340d0,360d0,380d0,400d0,420d0,440d0,460d0,
     .480d0,500d0,520d0,540d0,560d0,580d0,600d0,620d0,640d0,660d0,680d0,
     .700d0,720d0,740d0,760d0,780d0,800d0,820d0,840d0,860d0,880d0,900d0,
     .920d0,940d0,960d0,980d0,1000d0,1020d0,1040d0,1060d0,1080d0,1100d0,
     .1120d0,1140d0,1160d0,1180d0,1200d0/
      DATA L2G/0.3622d0,0.6746d0,0.1476d0,0.3070d0,0.2424d0,0.6799d0,
     .0.2007d0,0.3143d0,0.1728d0,0.2367d0,0.1635d0,0.1212d0,0.1222d0,
     .0.8851d-1,0.1129d0,0.9651d-1,0.9882d-1,0.8245d-1,0.6359d-1,
     .0.7385d-1,0.8509d-1,0.5348d-1,0.7502d-1,0.4497d-1,0.4641d-1,
     .0.3257d-1,0.3551d-1,0.3387d-1,0.4190d-1,0.3812d-1,0.3010d-1,
     .0.3723d-1,0.3361d-1,0.2340d-1,0.1862d-1,0.1391d-1,0.1306d-1,
     .0.1748d-1,0.1804d-1,0.1735d-1,0.1668d-1,0.1391d-1,0.1359d-1,
     .0.1402d-1,0.1359d-1,0.1402d-1,0.1348d-1,0.1425d-1,0.1370d-1,
     .0.1459d-1,0.1494d-1,0.1494d-1,0.1459d-1,0.1380d-1,0.1236d-1,
     .0.1198d-1,0.1090d-1/
*   fig. 6b VBF
      DATA M2V/200d0,205d0,220d0,225d0,235d0,250d0,260d0,275d0,290d0,
     .300d0,340d0,360d0,380d0,400d0,420d0,440d0,460d0,480d0,500d0,520d0,
     .540d0,560d0,580d0,600d0,620d0,640d0,660d0,680d0,700d0,720d0,740d0,
     .760d0,780d0,800d0,820d0,840d0,860d0,880d0,900d0,920d0,940d0,960d0,
     .980d0,1000d0,1020d0,1040d0,1060d0,1080d0,1100d0,1120d0,1140d0,
     .1160d0,1180d0,1200d0/
      DATA L2V/0.1414d0,0.1792d0,0.9537d-1,0.2535d0,0.4134d0,0.3040d0,
     .0.1152d0,0.6957d-1,0.9097d-1,0.7352d-1,0.6481d-1,0.6903d-1,
     .0.7237d-1,0.7769d-1,0.5942d-1,0.5035d-1,0.5579d-1,0.5803d-1,
     .0.4438d-1,0.5363d-1,0.4006d-1,0.4580d-1,0.4727d-1,0.5321d-1,
     .0.4134d-1,0.4267d-1,0.4580d-1,0.4617d-1,0.5406d-1,0.4653d-1,
     .0.3315d-1,0.3237d-1,0.3016d-1,0.2809d-1,0.3064d-1,0.2596d-1,
     .0.2288d-1,0.2148d-1,0.1806d-1,0.1723d-1,0.1778d-1,0.1723d-1,
     .0.1736d-1,0.1682d-1,0.1736d-1,0.1669d-1,0.1778d-1,0.1806d-1,
     .0.1806d-1,0.1764d-1,0.1682d-1,0.1518d-1,0.1471d-1,0.1360d-1/

*   ATLAS 1808.02380 H->VV
*   fig. 5a ggF
      DATA M3G/300d0,400d0,500d0,600d0,700d0,800d0,900d0,1000d0,1200d0,
     .1400d0,1600d0,1800d0,2000d0,2200d0,2400d0,2600d0,2800d0,3000d0/
      DATA L3G/0.3674d0,0.3230d0,0.1187d0,0.4686d-1,0.5406d-1,0.2073d-1,
     .0.1849d-1,0.2133d-1,0.1239d-1,0.9047d-2,0.5890d-2,0.4822d-2,
     .0.2533d-2,0.1797d-2,0.1876d-2,0.1293d-2,0.1187d-2,0.1293d-2/
*   fig. 5b VBF
      DATA M3V/500d0,600d0,700d0,800d0,1200d0,1600d0,1800d0,2000d0,
     .2200d0,2400d0,2600d0,2800d0,3000d0/
      DATA L3V/0.1331d0,0.8667d-1,0.3184d-1,0.1721d-1,0.1089d-1,
     .0.1014d-1,0.7841d-2,0.3571d-2,0.2462d-2,0.2533d-2,0.2840d-2,
     .0.2965d-2,0.3139d-2/

*   CMS 2109.06055 H->WW
      DATA M4/1000d0,1100d0,1200d0,1300d0,1400d0,1500d0,1600d0,1700d0,
     .1800d0,1900d0,2000d0,2100d0,2200d0,2300d0,2400d0,2500d0,2600d0,
     .2700d0,2800d0,2900d0,3000d0,3100d0,3200d0,3300d0,3400d0,3500d0,
     .3600d0,3700d0,3800d0,3900d0,4000d0,4100d0,4200d0,4300d0,4400d0,
     .4500d0/
*   fig. 7a ggF
      DATA L4G/0.18042d-1,0.15276d-1,0.88625d-2,0.38633d-2,0.24359d-2,
     .0.15276d-2,0.13739d-2,0.13906d-2,0.14461d-2,0.12212d-2,0.95475d-3,
     .0.77436d-3,0.66599d-3,0.76401d-3,0.7903d-3,0.69625d-3,0.5706d-3,
     .0.48255d-3,0.46988d-3,0.50412d-3,0.50732d-3,0.50487d-3,0.47801d-3,
     .0.4466d-3,0.39591d-3,0.3359d-3,0.29024d-3,0.25351d-3,0.23168d-3,
     .0.22228d-3,0.20523d-3,0.18265d-3,0.14839d-3,0.1231d-3,0.11303d-3,
     .0.11388d-3/
*   fig. 7b VBF
      DATA L4V/0.8595d-2,0.5847d-2,0.21881d-2,0.10402d-2,0.77083d-3,
     .0.63397d-3,0.54256d-3,0.51527d-3,0.61421d-3,0.68039d-3,0.73358d-3,
     .0.6557d-3,0.50355d-3,0.36764d-3,0.33695d-3,0.30533d-3,0.26351d-3,
     .0.29013d-3,0.34464d-3,0.35875d-3,0.36379d-3,0.36429d-3,0.35523d-3,
     .0.32926d-3,0.29349d-3,0.26006d-3,0.23557d-3,0.22171d-3,0.20718d-3,
     .0.20186d-3,0.1977d-3,0.19035d-3,0.17117d-3,0.14697d-3,0.14633d-3,
     .0.15863d-3/

*   CMS 2109.08268 H->ZZ
      DATA M5/1000d0,1200d0,1400d0,1600d0,1800d0,2000d0,2500d0,3000d0,
     .3500d0,4000d0,4500d0/
*   fig. 8a ggF
      DATA L5G/17.276d-3,7.7785d-3,4.8242d-3,2.8842d-3,1.9181d-3,
     .1.3698d-3,0.50885d-3,0.32886d-3,0.23925d-3,0.20666d-3,0.18949d-3/
*   fig. 8b VBF
      DATA L5V/7.5522d-3,2.5987d-3,2.1118d-3,1.6235d-3,0.98156d-3,
     .0.66831d-3,0.34251d-3,0.25988d-3,0.20052d-3,0.17619d-3,0.16417d-3/

*   CMS 2210.00043 fig 7a ggF
      DATA M6/1.3d3,1.4d3,1.5d3,1.6d3,1.7d3,1.8d3,1.9d3,2.0d3,2.1d3,
     .2.2d3,2.3d3,2.4d3,2.5d3,2.6d3,2.7d3,2.8d3,2.9d3,3.0d3,3.1d3,
     .3.2d3,3.3d3,3.4d3,3.5d3,3.6d3,3.7d3,3.8d3,3.9d3,4.0d3,4.1d3,
     .4.2d3,4.3d3,4.4d3,4.5d3,4.6d3,4.7d3,4.8d3,4.9d3,5.0d3,5.1d3,
     .5.2d3,5.3d3,5.4d3,5.5d3,5.6d3,5.7d3,5.8d3,5.9d3,6.0d3/
      DATA L6/0.15491d-1,0.16078d-1,0.14531d-1,0.12432d-1,0.1029d-1,
     .0.81815d-2,0.77797d-2,0.85144d-2,0.73271d-2,0.53995d-2,0.16809d-2,
     .0.22025d-2,0.25633d-2,0.23658d-2,0.21641d-2,0.20099d-2,0.19128d-2,
     .0.19728d-2,0.19226d-2,0.16198d-2,0.10852d-2,0.7885d-3,0.68754d-3,
     .0.56331d-3,0.41948d-3,0.32354d-3,0.30818d-3,0.31221d-3,0.30505d-3,
     .0.28427d-3,0.25948d-3,0.2396d-3,0.23552d-3,0.23696d-3,0.24895d-3,
     .0.25245d-3,0.23965d-3,0.21617d-3,0.19027d-3,0.16926d-3,0.15557d-3,
     .0.14782d-3,0.14488d-3,0.1425d-3,0.14051d-3,0.14021d-3,0.14077d-3,
     .0.14178d-3/

      PROB(50)=0d0

      DO I=1,3
       MH=SMASS(I)

       J=1
       DOWHILE(M1G(J).LE.MH .AND. J.LT.N1G)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M1G(N1G)) THEN
          SIG=CJ(I)**2*ggF13(MH)*BRWW(I)
          L=L1G(J-1)+(MH-M1G(J-1))/(M1G(J)-M1G(J-1))*(L1G(J)-L1G(J-1))
          PROB(50)=PROB(50)+DDIM(SIG/L,1d0)
       ENDIF

       J=1
       DOWHILE(M1V(J).LE.MH .AND. J.LT.N1V)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M1V(N1V)) THEN
          SIG=CV(I)**2*VBF13(MH)*BRWW(I)
          L=L1V(J-1)+(MH-M1V(J-1))/(M1V(J)-M1V(J-1))*(L1V(J)-L1V(J-1))
          PROB(50)=PROB(50)+DDIM(SIG/L,1d0)
       ENDIF

       J=1
       DOWHILE(M2G(J).LE.MH .AND. J.LT.N2G)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M2G(N2G)) THEN
          SIG=CJ(I)**2*ggF13(MH)*BRWW(I)
          L=L2G(J-1)+(MH-M2G(J-1))/(M2G(J)-M2G(J-1))*(L2G(J)-L2G(J-1))
          PROB(50)=PROB(50)+DDIM(SIG/L,1d0)
       ENDIF

       J=1
       DOWHILE(M2V(J).LE.MH .AND. J.LT.N2V)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M2V(N2V)) THEN
          SIG=CV(I)**2*VBF13(MH)*BRWW(I)
          L=L2V(J-1)+(MH-M2V(J-1))/(M2V(J)-M2V(J-1))*(L2V(J)-L2V(J-1))
          PROB(50)=PROB(50)+DDIM(SIG/L,1d0)
       ENDIF

       J=1
       DOWHILE(M3G(J).LE.MH .AND. J.LT.N3G)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M3G(N3G)) THEN
          SIG=CJ(I)**2*ggF13(MH)*BRWW(I)
          L=L3G(J-1)+(MH-M3G(J-1))/(M3G(J)-M3G(J-1))*(L3G(J)-L3G(J-1))
          PROB(50)=PROB(50)+DDIM(SIG/L,1d0)
       ENDIF

       J=1
       DOWHILE(M3V(J).LE.MH .AND. J.LT.N3V)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M3V(N3V)) THEN
          SIG=CV(I)**2*VBF13(MH)*BRWW(I)
          L=L3V(J-1)+(MH-M3V(J-1))/(M3V(J)-M3V(J-1))*(L3V(J)-L3V(J-1))
          PROB(50)=PROB(50)+DDIM(SIG/L,1d0)
       ENDIF

       J=1
       DOWHILE(M4(J).LE.MH .AND. J.LT.N4)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M4(N4)) THEN
          SIG=BRWW(I)
          L=L4G(J-1)+(MH-M4(J-1))/(M4(J)-M4(J-1))*(L4G(J)-L4G(J-1))
          PROB(50)=PROB(50)+DDIM(CJ(I)**2*ggF13(MH)*SIG/L,1d0)
          L=L4V(J-1)+(MH-M4(J-1))/(M4(J)-M4(J-1))*(L4V(J)-L4V(J-1))
          PROB(50)=PROB(50)+DDIM(CV(I)**2*VBF13(MH)*SIG/L,1d0)
       ENDIF

       J=1
       DOWHILE(M5(J).LE.MH .AND. J.LT.N5)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M5(N5)) THEN
          SIG=BRZZ(I)
          L=L5G(J-1)+(MH-M5(J-1))/(M5(J)-M5(J-1))*(L5G(J)-L5G(J-1))
          PROB(50)=PROB(50)+DDIM(CJ(I)**2*ggF13(MH)*SIG/L,1d0)
          L=L5V(J-1)+(MH-M5(J-1))/(M5(J)-M5(J-1))*(L5V(J)-L5V(J-1))
          PROB(50)=PROB(50)+DDIM(CV(I)**2*VBF13(MH)*SIG/L,1d0)
       ENDIF

       J=1
       DOWHILE(M6(J).LE.MH .AND. J.LT.N6)
         J=J+1
       ENDDO
       IF(J.GE.2 .AND. MH.LE.M6(N6)) THEN
          SIG=CJ(I)**2*(BRWW(I)+BRZZ(I))*ggF13(MH)
          L=L6(J-1)+(MH-M6(J-1))/(M6(J)-M6(J-1))*(L6(J)-L6(J-1))
          PROB(50)=PROB(50)+DDIM(SIG/L,1d0)
       ENDIF

      ENDDO

      END


      SUBROUTINE LIGHTHA_GAMGAM(PROB)

* Constraints from ggF->H/A->gamgam
* from ATLAS-CONF-2018-025 + 2407.07546
* and CMS-PAS-HIG-17-013 + 2405.18149
* PROB(53) =/= 0: excluded by ggF->H/A->gamgam

      IMPLICIT NONE

      INTEGER I,J,NHGG1,NHGG2,NHGG3,NHGG4,NHGG5,NHGG6

      DOUBLE PRECISION PROB(*),MH,ggF8,ggF13,VBF13,SHGG,LHGG
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION HGG1(300,2),HGG2(300,2),HGG3(300,2)
      DOUBLE PRECISION HGG4(400,2),HGG5(400,2),HGG6(100,2)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/LHCHGG/HGG1,HGG2,HGG3,HGG4,HGG5,HGG6,
     .      NHGG1,NHGG2,NHGG3,NHGG4,NHGG5,NHGG6

      PROB(53)=0d0

* Loop over 5 Higgses
      DO I=1,5
        IF(I.LE.3) MH=SMASS(I)
        IF(I.GE.4) MH=PMASS(I-3)

* CMS 8 TeV limit:
        J=1
        DOWHILE(HGG1(J,1).LE.MH .AND. J.LT.NHGG1)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG1(NHGG1,1)) THEN
          LHGG=HGG1(J-1,2)+(MH-HGG1(J-1,1))/
     .     (HGG1(J,1)-HGG1(J-1,1))*(HGG1(J,2)-HGG1(J-1,2))
          SHGG=CJ(I)**2*BRGG(I)*ggF8(MH)
          PROB(53)=PROB(53)+DDIM(SHGG/LHGG,1d0)
        ENDIF

* CMS 13 TeV limit:
        J=1
        DOWHILE(HGG2(J,1).LE.MH .AND. J.LT.NHGG2)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG2(NHGG2,1)) THEN
          LHGG=HGG2(J-1,2)+(MH-HGG2(J-1,1))/
     .     (HGG2(J,1)-HGG2(J-1,1))*(HGG2(J,2)-HGG2(J-1,2))
          SHGG=CJ(I)**2*BRGG(I)*ggF13(MH)
          PROB(53)=PROB(53)+DDIM(SHGG/LHGG,1d0)
        ENDIF

* ATLAS 13 TeV limit:
        J=1
        DOWHILE(HGG3(J,1).LE.MH .AND. J.LT.NHGG3)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG3(NHGG3,1)) THEN
          LHGG=HGG3(J-1,2)+(MH-HGG3(J-1,1))/
     .     (HGG3(J,1)-HGG3(J-1,1))*(HGG3(J,2)-HGG3(J-1,2))
          SHGG=CJ(I)**2*BRGG(I)*ggF13(MH)
          PROB(53)=PROB(53)+DDIM(SHGG/LHGG,1d0)
        ENDIF

* New CMS 13 TeV limit (2405.18149):
        J=1
        DOWHILE(HGG4(J,1).LE.MH .AND. J.LT.NHGG4)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG4(NHGG4,1)) THEN
          LHGG=HGG4(J-1,2)+(MH-HGG4(J-1,1))/
     .     (HGG4(J,1)-HGG4(J-1,1))*(HGG4(J,2)-HGG4(J-1,2))
          SHGG=CJ(I)**2*BRGG(I)*ggF13(MH)
          PROB(53)=PROB(53)+DDIM(SHGG/LHGG,1d0)
        ENDIF

        IF(I.LE.3)THEN
        J=1
        DOWHILE(HGG5(J,1).LE.MH .AND. J.LT.NHGG5)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG5(NHGG5,1)) THEN
          LHGG=HGG5(J-1,2)+(MH-HGG5(J-1,1))/
     .     (HGG5(J,1)-HGG5(J-1,1))*(HGG5(J,2)-HGG5(J-1,2))
          SHGG=CV(I)**2*BRGG(I)*VBF13(MH)
          PROB(53)=PROB(53)+DDIM(SHGG/LHGG,1d0)
        ENDIF
        ENDIF

* New ATLAS 13 TeV limit (2407.07546):
        J=1
        DOWHILE(HGG6(J,1).LE.MH .AND. J.LT.NHGG6)
          J=J+1
        ENDDO
        IF(J.GT.1 .AND. MH.LE.HGG6(NHGG6,1)) THEN
          LHGG=HGG6(J-1,2)+(MH-HGG6(J-1,1))/
     .     (HGG6(J,1)-HGG6(J-1,1))*(HGG6(J,2)-HGG6(J-1,2))
          SHGG=CJ(I)**2*BRGG(I)*ggF13(MH)
          PROB(53)=PROB(53)+DDIM(SHGG/LHGG,1d0)
        ENDIF

      ENDDO

      END


      SUBROUTINE HSM_AA_LEPTONS(PROB)

*     PROB(52) =/= 0: excluded 

      IMPLICIT NONE

      INTEGER NHAATAUS0,NHAATAUS1,NHAATAUS2,NHAATAUS3,NHAAMUS0
      INTEGER NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4,NHAAMUS5
      INTEGER NHAAMUTAU1,NHAAMUTAU2,NHAAMUB1,NHAAMUB2
      INTEGER NHAATAUB,NHAAGJ,NHAABS,NHAABS2,NHAAMUTAU3
      INTEGER NHAAMUB3,NHAAMUB4,NHAAMUS6,NHAAMUS7,NHAAMUS8
      INTEGER NHAATAUS4,NHAABS3
      INTEGER I,J,K

      DOUBLE PRECISION PROB(*),D1(3),D2(3),C2BRHAA,MH,MA
      DOUBLE PRECISION LOWBOUND,HIGHBOUND,LIMIT,ggF8,ggF13
      DOUBLE PRECISION HAATAUS0(100,2),HAATAUS1(300,2),HAATAUS2(100,2)
      DOUBLE PRECISION HAATAUS3(100,2),HAAMUS0(100,2),HAAMUS1(100,2)
      DOUBLE PRECISION HAAMUS2(100,2),HAAMUS3(100,2),HAAMUS4(100,4)
      DOUBLE PRECISION HAAMUS5(100,4),HAAMUTAU1(100,2),HAAMUTAU2(100,2)
      DOUBLE PRECISION HAAMUB1(100,2),HAAMUB2(200,2),HAATAUB(100,2)
      DOUBLE PRECISION HAAGJ(100,2),HAABS(100,2),HAABS2(100,2)
      DOUBLE PRECISION HAAMUTAU3(200,2),HAAMUB3(100,2),HAAMUB4(100,2)
      DOUBLE PRECISION HAAMUS6(200,2),HAAMUS7(200,2),HAAMUS8(400,2)
      DOUBLE PRECISION HAATAUS4(100,2),HAABS3(100,2)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION MHmin,MHmax

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/LHCHAA/HAATAUS0,HAATAUS1,HAATAUS2,HAATAUS3,
     .      HAAMUS0,HAAMUS1,HAAMUS2,HAAMUS3,HAAMUS4,HAAMUS5,
     .      HAAMUTAU1,HAAMUTAU2,HAAMUB1,HAAMUB2,HAATAUB,HAAGJ,
     .      HAABS,HAABS2,HAAMUTAU3,HAAMUB3,HAAMUB4,HAAMUS6,
     .      HAAMUS7,HAAMUS8,HAATAUS4,HAABS3,
     .      NHAATAUS0,NHAATAUS1,NHAATAUS2,NHAATAUS3,NHAAMUS0,
     .      NHAAMUS1,NHAAMUS2,NHAAMUS3,NHAAMUS4,NHAAMUS5,
     .      NHAAMUTAU1,NHAAMUTAU2,NHAAMUB1,NHAAMUB2,
     .      NHAATAUB,NHAAGJ,NHAABS,NHAABS2,NHAAMUTAU3,
     .      NHAAMUB3,NHAAMUB4,NHAAMUS6,NHAAMUS7,NHAAMUS8,
     .      NHAATAUS4,NHAABS3

      PROB(52)=0d0

* HSM -> a1a1 or h1h1
* Determining the SM-like Higgs and its coupling*BR

      DO I=1,3
       D1(I)=DDIM(SMASS(I)/MHMAX,1d0)-DDIM(1d0,SMASS(I)/MHMIN)
      ENDDO

      DO K=1,2

      MH=0d0
      C2BRHAA=0d0
      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CJ(I)**2*BRHAA(I,1)
         MH=MH+SMASS(I)*CJ(I)**2*BRHAA(I,1)
        ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CJ(I)**2*BRHHH(I-1)
         MH=MH+SMASS(I)*CJ(I)**2*BRHHH(I-1)
        ENDIF
       ENDDO
      ENDIF
      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

* ggF->HSM->AA->4mu from CMS-PAS-HIG-13-010, 0.25GeV < M_A < 3.55GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS0(1,1)
     . .AND. HAAMUS0(I,1).LE.MA .AND. I.LT.NHAAMUS0)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS0(NHAAMUS0,1))THEN
       LIMIT=HAAMUS0(I-1,2)
     .    +(MA-HAAMUS0(I-1,1))/(HAAMUS0(I,1)-HAAMUS0(I-1,1))
     .         *(HAAMUS0(I,2)-HAAMUS0(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4mu from 2110.13673 (ATLAS), 1.5GeV < M_A < 60GeV
* normalized to SM XS

      IF((MA.GT.1.5d0).AND.(MA.LT.2d0).OR.
     .   (MA.GT.4.5d0).AND.(MA.LT.8.5d0).OR.
     .   (MA.GT.8.5d0).AND.(MA.LT.60d0))THEN
      I=1
      DOWHILE(MA.GE.HAAMUS6(1,1)
     . .AND. HAAMUS6(I,1).LE.MA .AND. I.LT.NHAAMUS6)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS6(NHAAMUS6,1))THEN
       LIMIT=HAAMUS6(I-1,2)
     .    +(MA-HAAMUS6(I-1,1))/(HAAMUS6(I,1)-HAAMUS6(I-1,1))
     .         *(HAAMUS6(I,2)-HAAMUS6(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF
      ENDIF

* ggF->HSM->AA->4mu from 2111.01299 (CMS), 4.2GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS7(1,1)
     . .AND. HAAMUS7(I,1).LE.MA .AND. I.LT.NHAAMUS7)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS7(NHAAMUS7,1))THEN
       LIMIT=HAAMUS7(I-1,2)
     .    +(MA-HAAMUS7(I-1,1))/(HAAMUS7(I,1)-HAAMUS7(I-1,1))
     .         *(HAAMUS7(I,2)-HAAMUS7(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

      I=1
      DOWHILE(MA.GE.HAAMUS8(1,1)
     . .AND. HAAMUS8(I,1).LE.MA .AND. I.LT.NHAAMUS8)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS8(NHAAMUS8,1))THEN
       LIMIT=HAAMUS8(I-1,2)
     .    +(MA-HAAMUS8(I-1,1))/(HAAMUS8(I,1)-HAAMUS8(I-1,1))
     .         *(HAAMUS8(I,2)-HAAMUS8(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from CMS-HIG-PAS-15-011, 20GeV < M_A < 62.5GeV
* BR(A->2mu) converted in BR(A->2tau), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS0(1,1)
     .  .AND. HAATAUS0(I,1).LE.MA .AND. I.LT.NHAATAUS0)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS0(NHAATAUS0,1))THEN
       LIMIT=HAATAUS0(I-1,2)
     .    +(MA-HAATAUS0(I-1,1))/(HAATAUS0(I,1)-HAATAUS0(I-1,1))
     .         *(HAATAUS0(I,2)-HAATAUS0(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1701.02032 (CMS), 15GeV < M_A < 62.5GeV
* BR(A->2tau) converted in BR(A->2mu), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS1(1,1)
     .  .AND. HAAMUS1(I,1).LE.MA .AND. I.LT.NHAAMUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS1(NHAAMUS1,1))THEN
       LIMIT=HAAMUS1(I-1,2)
     .    +(MA-HAAMUS1(I-1,1))/(HAAMUS1(I,1)-HAAMUS1(I-1,1))
     .         *(HAAMUS1(I,2)-HAAMUS1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1701.02032 (CMS), 15GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUTAU1(1,1)
     .  .AND. HAAMUTAU1(I,1).LE.MA .AND. I.LT.NHAAMUTAU1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUTAU1(NHAAMUTAU1,1))THEN
       LIMIT=HAAMUTAU1(I-1,2)
     .    +(MA-HAAMUTAU1(I-1,1))/(HAAMUTAU1(I,1)-HAAMUTAU1(I-1,1))
     .         *(HAAMUTAU1(I,2)-HAAMUTAU1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRLL(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1805.04865 (CMS), 15GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUTAU2(1,1)
     .  .AND. HAAMUTAU2(I,1).LE.MA .AND. I.LT.NHAAMUTAU2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUTAU2(NHAAMUTAU2,1))THEN
       LIMIT=HAAMUTAU2(I-1,2)
     .    +(MA-HAAMUTAU2(I-1,1))/(HAAMUTAU2(I,1)-HAAMUTAU2(I-1,1))
     .         *(HAAMUTAU2(I,2)-HAAMUTAU2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRLL(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 2005.08694 (CMS), 3.5GeV < M_A < 21GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUTAU3(1,1)
     .  .AND. HAAMUTAU3(I,1).LE.MA .AND. I.LT.NHAAMUTAU3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUTAU3(NHAAMUTAU3,1))THEN
       LIMIT=HAAMUTAU3(I-1,2)
     .    +(MA-HAAMUTAU3(I-1,1))/(HAAMUTAU3(I,1)-HAAMUTAU3(I-1,1))
     .         *(HAAMUTAU3(I,2)-HAAMUTAU3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRLL(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2tau from 1505.01609 (ATLAS), 3.7GeV < M_A < 50GeV
* BR(A->2mu) converted in BR(A->2tau), normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS1(1,1)
     .  .AND. HAATAUS1(I,1).LE.MA .AND. I.LT.NHAATAUS1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS1(NHAATAUS1,1))THEN
       LIMIT=HAATAUS1(I-1,2)
     .    +(MA-HAATAUS1(I-1,1))/(HAATAUS1(I,1)-HAATAUS1(I-1,1))
     .         *(HAATAUS1(I,2)-HAATAUS1(I-1,2))
       PROB(52)=PROB(52)+DDIM(C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from 1510.06534 (CMS), 4GeV < M_A < 8GeV

      I=1
      DOWHILE(MA.GE.HAATAUS2(1,1)
     .  .AND. HAATAUS2(I,1).LE.MA .AND. I.LT.NHAATAUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS2(NHAATAUS2,1))THEN
       LIMIT=HAATAUS2(I-1,2)
     .    +(MA-HAATAUS2(I-1,1))/(HAATAUS2(I,1)-HAATAUS2(I-1,1))
     .         *(HAATAUS2(I,2)-HAATAUS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(ggF8(MH)*C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from 1510.06534 (CMS), 4GeV < M_A < 8GeV
* BR(A->2tau) converted in BR(A->2mu), Fig. 7 of 1701.02032
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS3(1,1)
     .  .AND. HAAMUS3(I,1).LE.MA .AND. I.LT.NHAAMUS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS3(NHAAMUS3,1))THEN
       LIMIT=HAAMUS3(I-1,2)
     .    +(MA-HAAMUS3(I-1,1))/(HAAMUS3(I,1)-HAAMUS3(I-1,1))
     .         *(HAAMUS3(I,2)-HAAMUS3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from CMS-PAS-HIG-14-022, 5GeV < M_A < 15GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS3(1,1)
     .  .AND. HAATAUS3(I,1).LE.MA .AND. I.LT.NHAATAUS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS3(NHAATAUS3,1))THEN
       LIMIT=HAATAUS3(I-1,2)
     .    +(MA-HAATAUS3(I-1,1))/(HAATAUS3(I,1)-HAATAUS3(I-1,1))
     .         *(HAATAUS3(I,2)-HAATAUS3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from 1907.07235 (CMS), 4GeV < M_A < 15GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUS4(1,1)
     .  .AND. HAATAUS4(I,1).LE.MA .AND. I.LT.NHAATAUS4)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUS4(NHAATAUS4,1))THEN
       LIMIT=HAATAUS4(I-1,2)
     .    +(MA-HAATAUS4(I-1,1))/(HAATAUS4(I,1)-HAATAUS4(I-1,1))
     .         *(HAATAUS4(I,2)-HAATAUS4(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRLL(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->4tau from CMS-PAS-HIG-14-022, 5GeV < M_A < 15GeV
* BR(A->2tau) converted in BR(A->2mu), Fig. 7 of 1701.02032
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUS2(1,1)
     .  .AND. HAAMUS2(I,1).LE.MA .AND. I.LT.NHAAMUS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUS2(NHAAMUS2,1))THEN
       LIMIT=HAAMUS2(I-1,2)
     .    +(MA-HAAMUS2(I-1,1))/(HAAMUS2(I,1)-HAAMUS2(I-1,1))
     .         *(HAAMUS2(I,2)-HAAMUS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2b from CMS-PAS-HIG-14-041, 25GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUB1(1,1)
     .  .AND. HAAMUB1(I,1).LE.MA .AND. I.LT.NHAAMUB1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUB1(NHAAMUB1,1))THEN
       LIMIT=HAAMUB1(I-1,2)
     .    +(MA-HAAMUB1(I-1,1))/(HAAMUB1(I,1)-HAAMUB1(I-1,1))
     .         *(HAAMUB1(I,2)-HAAMUB1(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRBB(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2b from 1807.00539 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUB2(1,1)
     .  .AND. HAAMUB2(I,1).LE.MA .AND. I.LT.NHAAMUB2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUB2(NHAAMUB2,1))THEN
       LIMIT=HAAMUB2(I-1,2)
     .    +(MA-HAAMUB2(I-1,1))/(HAAMUB2(I,1)-HAAMUB2(I-1,1))
     .         *(HAAMUB2(I,2)-HAAMUB2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRBB(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2b from 2110.00313 (ATLAS), 16GeV < M_A < 62GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUB3(1,1)
     .  .AND. HAAMUB3(I,1).LE.MA .AND. I.LT.NHAAMUB3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUB3(NHAAMUB3,1))THEN
       LIMIT=HAAMUB3(I-1,2)
     .    +(MA-HAAMUB3(I-1,1))/(HAAMUB3(I,1)-HAAMUB3(I-1,1))
     .         *(HAAMUB3(I,2)-HAAMUB3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRBB(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2mu2b from 1812.06359 (CMS), 20GeV < M_A < 62.5GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAMUB4(1,1)
     .  .AND. HAAMUB4(I,1).LE.MA .AND. I.LT.NHAAMUB4)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAMUB4(NHAAMUB4,1))THEN
       LIMIT=HAAMUB4(I-1,2)
     .    +(MA-HAAMUB4(I-1,1))/(HAAMUB4(I,1)-HAAMUB4(I-1,1))
     .         *(HAAMUB4(I,2)-HAAMUB4(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRMM(7-3*K)*BRBB(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2tau2b from 1805.10191 (CMS), 15GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAATAUB(1,1)
     .  .AND. HAATAUB(I,1).LE.MA .AND. I.LT.NHAATAUB)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAATAUB(NHAATAUB,1))THEN
       LIMIT=HAATAUB(I-1,2)
     .    +(MA-HAATAUB(I-1,1))/(HAATAUB(I,1)-HAATAUB(I-1,1))
     .         *(HAATAUB(I,2)-HAATAUB(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRLL(7-3*K)*BRBB(7-3*K)/LIMIT,1d0)
      ENDIF

* ggF->HSM->AA->2gamma2jet from 1803.11145 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      I=1
      DOWHILE(MA.GE.HAAGJ(1,1)
     .  .AND. HAAGJ(I,1).LE.MA .AND. I.LT.NHAAGJ)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAAGJ(NHAAGJ,1))THEN
       LIMIT=HAAGJ(I-1,2)
     .    +(MA-HAAGJ(I-1,1))/(HAAGJ(I,1)-HAAGJ(I-1,1))
     .         *(HAAGJ(I,2)-HAAGJ(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRGG(7-3*K)*BRJJ(7-3*K)/LIMIT,1d0)
      ENDIF

      ENDIF

      ENDDO

* pp->VHSM->AA->4b from 1806.07355 (ATLAS), 20GeV < M_A < 60GeV
* normalized to SM XS

      DO K=1,2

      C2BRHAA=0d0
      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CV(I)**2*BRHAA(I,1)
        ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CV(I)**2*BRHHH(I-1)
        ENDIF
       ENDDO
      ENDIF

      IF(C2BRHAA.NE.0d0)THEN

      I=1
      DOWHILE(MA.GE.HAABS(1,1)
     .  .AND. HAABS(I,1).LE.MA .AND. I.LT.NHAABS)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAABS(NHAABS,1))THEN
       LIMIT=HAABS(I-1,2)
     .    +(MA-HAABS(I-1,1))/(HAABS(I,1)-HAABS(I-1,1))
     .         *(HAABS(I,2)-HAABS(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRBB(7-3*K)**2/LIMIT,1d0)
      ENDIF

      ENDIF
      ENDDO

* pp->VHSM->AA->4b from 2005.12236 (ATLAS), 15GeV < M_A < 30GeV
* normalized to SM XS

      DO K=1,2

      C2BRHAA=0d0
      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CV(I)**2*BRHAA(I,1)
        ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CV(I)**2*BRHHH(I-1)
        ENDIF
       ENDDO
      ENDIF

      IF(C2BRHAA.NE.0d0)THEN

      I=1
      DOWHILE(MA.GE.HAABS2(1,1)
     .  .AND. HAABS2(I,1).LE.MA .AND. I.LT.NHAABS2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAABS2(NHAABS2,1))THEN
       LIMIT=HAABS2(I-1,2)
     .    +(MA-HAABS2(I-1,1))/(HAABS2(I,1)-HAABS2(I-1,1))
     .         *(HAABS2(I,2)-HAABS2(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRBB(7-3*K)**2/LIMIT,1d0)
      ENDIF

      ENDIF
      ENDDO

* pp->VHSM->AA->4b from 2403.10341 (CMS), 15GeV < M_A < 60GeV
* normalized to SM XS

      DO K=1,2

      C2BRHAA=0d0
      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CV(I)**2*BRHAA(I,1)
        ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         C2BRHAA=C2BRHAA+CV(I)**2*BRHHH(I-1)
        ENDIF
       ENDDO
      ENDIF

      IF(C2BRHAA.NE.0d0)THEN

      I=1
      DOWHILE(MA.GE.HAABS3(1,1)
     .  .AND. HAABS3(I,1).LE.MA .AND. I.LT.NHAABS3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. MA.LE.HAABS3(NHAABS3,1))THEN
       LIMIT=HAABS3(I-1,2)
     .    +(MA-HAABS3(I-1,1))/(HAABS3(I,1)-HAABS3(I-1,1))
     .         *(HAABS3(I,2)-HAABS3(I-1,2))
       PROB(52)=PROB(52)
     .          +DDIM(C2BRHAA*BRBB(7-3*K)**2/LIMIT,1d0)
      ENDIF

      ENDIF
      ENDDO

* ggF -> H -> a1a1 or h1h1 -> 4mu from 1506.00424 (CMS 8TeV)
* for 0.25GeV < M_A < 3.55GeV, 86GeV < M_H < 150GeV

      DO I=1,3
       D1(I)=DDIM(SMASS(I)/HAAMUS4(NHAAMUS4,1),1d0)
     .      -DDIM(1d0,SMASS(I)/HAAMUS4(1,1))
      ENDDO

        DO K=1,2

      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
         D2(I)=CJ(I)**2*BRHAA(I,1)
        ELSE
         D2(I)=0d0
        ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         D2(I)=CJ(I)**2*BRHHH(I-1)
        ELSE
         D2(I)=0d0
        ENDIF
       ENDDO
      ENDIF

       DO J=K,3

      MH=0d0
      C2BRHAA=0d0
      DO I=K,3
       IF(DABS(SMASS(I)-SMASS(J)).LT.3d0)THEN
        C2BRHAA=C2BRHAA+D2(I)
        MH=MH+SMASS(I)*D2(I)
       ENDIF
      ENDDO
      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

       I=1
       DOWHILE(MH.GE.HAAMUS4(1,1)
     .   .AND. HAAMUS4(I,1).LE.MH .AND. I.LT.NHAAMUS4)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.HAAMUS4(NHAAMUS4,1))THEN
        IF(MA.GE.0.25d0 .AND. MA.LE.3.55d0)THEN
         IF(MA.LT.2d0)THEN
          LOWBOUND=HAAMUS4(I-1,2)+(MA-0.25d0)/(2d0-0.25d0)
     .            *(HAAMUS4(I-1,3)-HAAMUS4(I-1,2))
          HIGHBOUND=HAAMUS4(I,2)+(MA-0.25d0)/(2d0-0.25d0)
     .             *(HAAMUS4(I,3)-HAAMUS4(I,2))
         ELSE
          LOWBOUND=HAAMUS4(I-1,3)+(MA-2d0)/(3.55d0-2d0)
     .            *(HAAMUS4(I-1,4)-HAAMUS4(I-1,3))
          HIGHBOUND=HAAMUS4(I,3)+(MA-2d0)/(3.55d0-2d0)
     .             *(HAAMUS4(I,4)-HAAMUS4(I,3))
         ENDIF
         LIMIT=(LOWBOUND
     .        +(MH-HAAMUS4(I-1,1))/(HAAMUS4(I,1)-HAAMUS4(I-1,1))
     .        *(HIGHBOUND-LOWBOUND))/1d3 ! translated in pb
         PROB(52)=PROB(52)
     .           +DDIM(ggF8(MH)*C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
        ENDIF
       ENDIF

      ENDIF
       ENDDO
        ENDDO

* ggF -> H -> a1a1 or h1h1 -> 4mu from CMS-PAS-HIG-18-003 (CMS 13TeV)
* for 0.25GeV < M_A < 3.55GeV, 90GeV < M_H < 150GeV

      DO I=1,3
       D1(I)=DDIM(SMASS(I)/150d0,1d0)-DDIM(1d0,SMASS(I)/90d0)
      ENDDO

        DO K=1,2

      IF(K.EQ.1)THEN
       MA=PMASS(1)
       DO I=1,3
        IF(D1(I).EQ.0d0)THEN
         D2(I)=CJ(I)**2*BRHAA(I,1)
        ELSE
         D2(I)=0d0
        ENDIF
       ENDDO
      ELSE
       MA=SMASS(1)
       DO I=2,3
        IF(D1(I).EQ.0d0)THEN
         D2(I)=CJ(I)**2*BRHHH(I-1)
        ELSE
         D2(I)=0d0
        ENDIF
       ENDDO
      ENDIF

       DO J=K,3

      MH=0d0
      C2BRHAA=0d0
      DO I=K,3
       IF(DABS(SMASS(I)-SMASS(J)).LT.3d0)THEN
        C2BRHAA=C2BRHAA+D2(I)
        MH=MH+SMASS(I)*D2(I)
       ENDIF
      ENDDO
      IF(C2BRHAA.NE.0d0)THEN
      MH=MH/C2BRHAA

       I=1
       DOWHILE(MA.GE.HAAMUS5(1,1)
     .   .AND. HAAMUS5(I,1).LE.MA .AND. I.LT.NHAAMUS5)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MA.LE.HAAMUS5(NHAAMUS5,1))THEN
        IF(MH.GE.90d0 .AND. MH.LE.150d0)THEN
         IF(MA.LT.125d0)THEN
          LOWBOUND=HAAMUS5(I-1,2)+(MH-90d0)/(125d0-90d0)
     .            *(HAAMUS5(I-1,3)-HAAMUS5(I-1,2))
          HIGHBOUND=HAAMUS5(I,2)+(MH-90d0)/(125d0-90d0)
     .             *(HAAMUS5(I,3)-HAAMUS5(I,2))
         ELSE
          LOWBOUND=HAAMUS5(I-1,3)+(MH-125d0)/(150d0-125d0)
     .            *(HAAMUS5(I-1,4)-HAAMUS5(I-1,3))
          HIGHBOUND=HAAMUS5(I,3)+(MH-2d0)/(150d0-125d0)
     .             *(HAAMUS5(I,4)-HAAMUS5(I,3))
         ENDIF
         LIMIT=(LOWBOUND
     .        +(MH-HAAMUS5(I-1,1))/(HAAMUS5(I,1)-HAAMUS5(I-1,1))
     .        *(HIGHBOUND-LOWBOUND))/1d3 ! translated in pb
         PROB(52)=PROB(52)
     .           +DDIM(ggF13(MH)*C2BRHAA*BRMM(7-3*K)**2/LIMIT,1d0)
        ENDIF
       ENDIF

      ENDIF
       ENDDO
        ENDDO

      END


      SUBROUTINE HSM_AA_PHOTONS(PROB)

*     H_125 -> AA -> photons from ATLAS 1509.05051
*     and CMS 2208.01469 + 2209.06197
*     PROB(63) =/= 0: excluded 

      IMPLICIT NONE

      INTEGER I,N,N2,N3
      PARAMETER (N=185,N2=73,N3=12)

      DOUBLE PRECISION PROB(*),S,SMAX
      DOUBLE PRECISION M(N),L(N),M2(N2),L2(N2),M3(N3),L3(N3)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION MHmin,MHmax

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSFIT/MHmin,MHmax

*   ATLAS 1509.05051 fig. 4b
      DATA M/1d1,1.0436d1,1.0744d1,1.0913d1,1.1094d1,1.1260d1,1.1594d1,
     .1.2120d1,1.2521d1,1.2750d1,1.3071d1,1.3438d1,1.3974d1,1.4343d1,
     .1.4674d1,1.5156d1,1.5500d1,1.5746d1,1.5950d1,1.6175d1,1.6394d1,
     .1.6727d1,1.7276d1,1.7734d1,1.8192d1,1.8765d1,1.9224d1,1.9637d1,
     .2.0140d1,2.0713d1,2.1101d1,2.1343d1,2.1479d1,2.1660d1,2.1796d1,
     .2.1946d1,2.2139d1,2.2325d1,2.2606d1,2.3004d1,2.3483d1,2.3861d1,
     .2.4098d1,2.4313d1,2.4513d1,2.4769d1,2.5011d1,2.5410d1,2.5945d1,
     .2.6327d1,2.6785d1,2.7051d1,2.7268d1,2.7520d1,2.7759d1,2.8236d1,
     .2.8561d1,2.8788d1,2.8902d1,2.9013d1,2.9125d1,2.9236d1,2.9355d1,
     .2.9568d1,2.9775d1,3.0112d1,3.0555d1,3.0804d1,3.1182d1,3.1697d1,
     .3.2282d1,3.2715d1,3.3201d1,3.3659d1,3.4173d1,3.4498d1,3.4710d1,
     .3.4989d1,3.5329d1,3.5802d1,3.6004d1,3.6211d1,3.6408d1,3.6568d1,
     .3.6797d1,3.7268d1,3.7594d1,3.7949d1,3.8185d1,3.8414d1,3.8643d1,
     .3.8865d1,3.9083d1,3.9331d1,3.9554d1,3.9789d1,4.0114d1,4.0625d1,
     .4.0935d1,4.1164d1,4.1378d1,4.1529d1,4.1729d1,4.1867d1,4.2046d1,
     .4.2188d1,4.2403d1,4.2651d1,4.3054d1,4.3492d1,4.3797d1,4.3930d1,
     .4.4077d1,4.4237d1,4.4360d1,4.4473d1,4.4584d1,4.4707d1,4.4820d1,
     .4.5059d1,4.5460d1,4.5976d1,4.6205d1,4.6320d1,4.6434d1,4.6546d1,
     .4.6648d1,4.6720d1,4.6870d1,4.6964d1,4.7110d1,4.7215d1,4.7444d1,
     .4.7934d1,4.8249d1,4.8470d1,4.8611d1,4.8823d1,4.9179d1,4.9799d1,
     .4.9974d1,5.0207d1,5.0437d1,5.0663d1,5.0900d1,5.1126d1,5.1455d1,
     .5.1876d1,5.2394d1,5.2604d1,5.2839d1,5.3058d1,5.3289d1,5.3413d1,
     .5.3535d1,5.3658d1,5.3881d1,5.4110d1,5.4485d1,5.5027d1,5.5371d1,
     .5.5600d1,5.5810d1,5.5922d1,5.6031d1,5.6139d1,5.6282d1,5.6418d1,
     .5.6613d1,5.6757d1,5.6959d1,5.7204d1,5.7547d1,5.7948d1,5.8521d1,
     .5.8922d1,5.9180d1,5.9467d1,5.9713d1,6.0055d1,6.0399d1,6.0641d1,
     .6.1099d1,6.1638d1,6.2d1/
      DATA L/3.0495d-4,3.3188d-4,3.5175d-4,3.7958d-4,4.1748d-4,
     .4.5015d-4,4.9288d-4,5.0986d-4,4.7180d-4,4.3749d-4,4.0261d-4,
     .3.8321d-4,3.9233d-4,4.1524d-4,4.5016d-4,4.6345d-4,4.3712d-4,
     .4.0954d-4,3.7648d-4,3.4511d-4,3.2145d-4,3.0275d-4,3.0490d-4,
     .3.1516d-4,3.2702d-4,3.3475d-4,3.2452d-4,3.1097d-4,2.9572d-4,
     .2.9321d-4,3.1935d-4,3.4496d-4,3.7242d-4,4.0961d-4,4.4202d-4,
     .4.7699d-4,5.2462d-4,5.7701d-4,6.1537d-4,6.5265d-4,6.3860d-4,
     .6.0308d-4,5.6361d-4,5.2482d-4,4.8616d-4,4.4202d-4,4.1415d-4,
     .3.8329d-4,3.8954d-4,4.1103d-4,4.0633d-4,3.7242d-4,3.4511d-4,
     .3.1981d-4,2.9690d-4,2.8802d-4,3.1981d-4,3.5175d-4,3.7958d-4,
     .4.0961d-4,4.4202d-4,4.7699d-4,5.1473d-4,5.6613d-4,6.2266d-4,
     .6.5711d-4,6.3463d-4,5.8810d-4,5.4809d-4,5.3743d-4,5.5929d-4,
     .5.8265d-4,5.7396d-4,5.6264d-4,5.8253d-4,6.3463d-4,6.8484d-4,
     .7.2922d-4,7.7615d-4,7.3902d-4,6.7193d-4,6.1092d-4,5.6613d-4,
     .5.2462d-4,4.8002d-4,4.5880d-4,4.9458d-4,5.4330d-4,5.8161d-4,
     .6.2988d-4,6.8101d-4,7.3839d-4,7.9749d-4,8.7313d-4,9.6472d-4,
     .1.0645d-3,1.1396d-3,1.1342d-3,1.0634d-3,9.8790d-4,9.1117d-4,
     .8.4436d-4,7.6770d-4,7.1142d-4,6.4682d-4,5.9940d-4,5.4498d-4,
     .5.0833d-4,4.8501d-4,5.0502d-4,5.5516d-4,5.9940d-4,6.4682d-4,
     .7.1142d-4,7.6770d-4,8.2844d-4,8.9399d-4,9.6472d-4,1.0410d-3,
     .1.1140d-3,1.1812d-3,1.1037d-3,1.0042d-3,9.2868d-4,8.6059d-4,
     .7.9749d-4,7.3902d-4,6.8484d-4,6.2266d-4,5.7701d-4,5.2462d-4,
     .4.8616d-4,4.4428d-4,4.3512d-4,4.7699d-4,5.2462d-4,5.6409d-4,
     .6.2243d-4,6.8064d-4,6.8484d-4,6.3463d-4,5.7701d-4,5.2462d-4,
     .4.7699d-4,4.3556d-4,4.0296d-4,3.7598d-4,3.6342d-4,3.8221d-4,
     .4.0961d-4,4.4202d-4,4.8616d-4,5.3411d-4,5.7701d-4,6.2266d-4,
     .6.7193d-4,7.4853d-4,8.1282d-4,8.8536d-4,8.6768d-4,7.9209d-4,
     .7.2098d-4,6.5926d-4,6.1092d-4,5.6613d-4,5.2462d-4,4.7699d-4,
     .4.4202d-4,4.0188d-4,3.7242d-4,3.4511d-4,3.1805d-4,2.9313d-4,
     .2.7738d-4,2.8251d-4,2.9963d-4,3.1981d-4,3.5175d-4,3.7911d-4,
     .4.1656d-4,4.5876d-4,4.9140d-4,5.0985d-4,4.9083d-4,4.6356d-4/

*   CMS 2208.01469 fig. 6
      DATA M2/15d0,15.5d0,16d0,16.5d0,17d0,17.5d0,18d0,18.5d0,19d0,
     .19.5d0,20d0,20.5d0,21d0,21.5d0,22d0,22.5d0,23d0,23.5d0,24d0,
     .24.5d0,25d0,25.5d0,26d0,26.5d0,27d0,27.5d0,28d0,28.5d0,29d0,
     .29.5d0,30d0,30.5d0,31d0,31.5d0,32d0,32.5d0,33d0,33.5d0,34d0,
     .34.5d0,35d0,35.5d0,36d0,36.5d0,37d0,37.5d0,38d0,38.5d0,39d0,
     .39.5d0,40d0,41d0,42d0,43d0,44d0,45d0,46d0,47d0,48d0,49d0,50d0,
     .51d0,52d0,53d0,54d0,55d0,56d0,57d0,58d0,59d0,60d0,61d0,62d0/
      DATA L2/0.79762d-3,0.81456d-3,0.69173d-3,0.73468d-3,0.67254d-3,
     .0.7172d-3,0.65982d-3,0.68493d-3,0.92037d-3,0.6506d-3,0.63782d-3,
     .0.87733d-3,1.1097d-3,1.0916d-3,0.79827d-3,1.0893d-3,1.2019d-3,
     .1.3994d-3,1.1016d-3,1.1763d-3,1.3079d-3,1.3632d-3,1.2891d-3,
     .1.2803d-3,1.285d-3,0.99275d-3,0.72587d-3,0.70764d-3,0.72875d-3,
     .0.9417d-3,0.93025d-3,0.81304d-3,0.93861d-3,0.96988d-3,0.93809d-3,
     .0.99021d-3,1.065d-3,0.73144d-3,0.71372d-3,0.85651d-3,0.86656d-3,
     .0.83829d-3,0.76955d-3,0.78722d-3,0.78339d-3,0.90652d-3,0.75615d-3,
     .0.772d-3,0.75378d-3,0.71779d-3,0.74616d-3,0.73923d-3,1.0491d-3,
     .0.9434d-3,1.0727d-3,0.84093d-3,0.65384d-3,0.55801d-3,0.63615d-3,
     .0.65291d-3,0.81987d-3,0.59349d-3,0.52921d-3,0.28462d-3,0.27096d-3,
     .0.44662d-3,0.62633d-3,0.49981d-3,0.36537d-3,0.24079d-3,0.3309d-3,
     .0.27697d-3,0.25991d-3/

*   CMS 2209.06197 fig. 2
      DATA M3/0.1d0,0.2d0,0.3d0,0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,
     .1.0d0,1.1d0,1.2d0/
      DATA L3/0.16804d-2,0.20714d-2,0.20585d-2,0.16535d-2,0.12822d-2,
     .0.96398d-3,0.85217d-3,0.11349d-2,0.13944d-2,0.15127d-2,0.12559d-2,
     .0.33493d-2/

      PROB(63)=0d0

      S=0d0
      DO I=1,3
       IF(SMASS(I).GE.MHMIN .AND. SMASS(I).LE.MHMAX)THEN
        S=S+CJ(I)**2*BRHAA(I,1)*BRGG(4)**2
       ENDIF
      ENDDO

      I=1
      DOWHILE(M(I).LE.PMASS(1) .AND. I.LT.N)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. PMASS(1).LE.M(N))THEN
       SMAX=L(I-1)+(PMASS(1)-M(I-1))/(M(I)-M(I-1))*(L(I)-L(I-1))
       PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
      ENDIF

      I=1
      DOWHILE(M2(I).LE.PMASS(1) .AND. I.LT.N2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. PMASS(1).LE.M2(N2))THEN
       SMAX=L2(I-1)+(PMASS(1)-M2(I-1))/(M2(I)-M2(I-1))*(L2(I)-L2(I-1))
       PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
      ENDIF

      I=1
      DOWHILE(M3(I).LE.PMASS(1) .AND. I.LT.N3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. PMASS(1).LE.M3(N3))THEN
       SMAX=L3(I-1)+(PMASS(1)-M3(I-1))/(M3(I)-M3(I-1))*(L3(I)-L3(I-1))
       PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
      ENDIF

      S=0d0
      DO I=2,3
       IF(SMASS(I).GE.MHMIN .AND. SMASS(I).LE.MHMAX)THEN
        S=S+CJ(I)**2*BRHHH(I-1)*BRGG(1)**2
       ENDIF
      ENDDO

      I=1
      DOWHILE(M(I).LE.SMASS(1) .AND. I.LT.N)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. SMASS(1).LE.M(N))THEN
       SMAX=L(I-1)+(SMASS(1)-M(I-1))/(M(I)-M(I-1))*(L(I)-L(I-1))
       PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
      ENDIF

      I=1
      DOWHILE(M2(I).LE.SMASS(1) .AND. I.LT.N2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. SMASS(1).LE.M2(N2))THEN
       SMAX=L2(I-1)+(SMASS(1)-M2(I-1))/(M2(I)-M2(I-1))*(L2(I)-L2(I-1))
       PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
      ENDIF

      I=1
      DOWHILE(M3(I).LE.SMASS(1) .AND. I.LT.N3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. SMASS(1).LE.M3(N3))THEN
       SMAX=L3(I-1)+(SMASS(1)-M3(I-1))/(M3(I)-M3(I-1))*(L3(I)-L3(I-1))
       PROB(63)=PROB(63)+DDIM(S/SMAX,1d0)
      ENDIF

      END


      SUBROUTINE HSM_ZA(PROB)

*     H_125 -> ZA from ATLAS 2004.01678 + 2110.13673
*     PROB(75) =/= 0: excluded 
*     H_125 -> ZA from CMS 2111.01299 + 2311.00130
*     PROB(76) =/= 0: excluded 

      IMPLICIT NONE

      INTEGER I,N1,N2,N3,N4,N5
      PARAMETER (N1=9,N2=16,N3=129,N4=223,N5=30)

      DOUBLE PRECISION PROB(*),S,SMAX,M1(N1),L1(N1),M2(N2),L2(N2)
      DOUBLE PRECISION M3(N3),L3(N3),M4(N4),L4(N4),M5(N5),L5(N5)
      DOUBLE PRECISION gg13,VB13,VH13,TH13
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION MHmin,MHmax

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/HIGGSFIT/MHmin,MHmax

*   ATLAS 2004.01678, tab. 1: A -> had and 0.5 < MA < 4 GeV
      DATA M1/.5d0,.75d0,1d0,1.5d0,2d0,2.5d0,3d0,3.5d0,4d0/
      DATA L1/17d0,20d0,18d0,20d0,24d0,33d0,68d0,120d0,340d0/

*   ATLAS 2110.13673, fig. 17b: Z -> ll + A -> mumu  and 15 < MA < 30 GeV
      DATA M2/15d0,16d0,17d0,18d0,19d0,20d0,21d0,22d0,23d0,
     .24d0,25d0,26d0,27d0,28d0,29d0,30d0/
      DATA L2/0.8843d-3,0.8012d-3,0.7169d-3,0.1239d-2,0.1764d-2,
     .0.1116d-2,0.9808d-3,0.5484d-3,0.3523d-3,0.3404d-3,0.6708d-3,
     .0.1143d-2,0.1116d-2,0.1065d-2,0.9872d-3,0.9059d-3/

*   CMS 2111.01299, fig. 4a: A -> mumu  and 4 < MA < 35 GeV
      DATA M3/0.42045d1,0.42255d1,0.42466d1,0.42678d1,0.42892d1,
     .0.43106d1,0.43322d1,0.43538d1,0.43756d1,0.43975d1,0.44195d1,
     .0.44416d1,0.44638d1,0.44861d1,0.45085d1,0.45311d1,0.45537d1,
     .0.45765d1,0.45994d1,0.46224d1,0.46455d1,0.46687d1,0.46921d1,
     .0.47155d1,0.47391d1,0.47628d1,0.47866d1,0.48105d1,0.48346d1,
     .0.48588d1,0.48831d1,0.49075d1,0.49320d1,0.49567d1,0.49815d1,
     .0.50064d1,0.50314d1,0.50565d1,0.50818d1,0.51072d1,0.51328d1,
     .0.51584d1,0.51842d1,0.52102d1,0.52362d1,0.52624d1,0.52887d1,
     .0.53151d1,0.53417d1,0.53684d1,0.53953d1,0.54222d1,0.54494d1,
     .0.54766d1,0.55040d1,0.55315d1,0.55592d1,0.55870d1,0.56149d1,
     .0.56430d1,0.56712d1,0.56995d1,0.57280d1,0.57567d1,0.57855d1,
     .0.58144d1,0.58435d1,0.58727d1,0.59020d1,0.59315d1,0.59612d1,
     .0.59910d1,0.60210d1,0.60511d1,0.60813d1,0.61117d1,0.61423d1,
     .0.61730d1,0.62039d1,0.62349d1,0.62661d1,0.62974d1,0.63289d1,
     .0.63605d1,0.63923d1,0.64243d1,0.64564d1,0.64887d1,0.65211d1,
     .0.65537d1,0.65865d1,0.66194d1,0.66525d1,0.66858d1,0.67192d1,
     .0.67528d1,0.67866d1,0.68205d1,0.68546d1,0.68889d1,0.69233d1,
     .0.69580d1,0.69928d1,0.70277d1,0.70629d1,0.70982d1,0.71337d1,
     .0.71693d1,0.72052d1,0.72412d1,0.72774d1,0.73138d1,0.73504d1,
     .0.73871d1,0.74241d1,0.74612d1,0.74985d1,0.75360d1,0.75736d1,
     .0.76115d1,0.76496d1,0.76878d1,0.77263d1,0.77649d1,0.78037d1,
     .0.78427d1,0.78819d1,0.79214d1,0.79610d1/
      DATA L3/0.16816d-4,0.18924d-4,0.19186d-4,0.24204d-4,0.23964d-4,
     .0.30845d-4,0.38356d-4,0.37152d-4,0.35454d-4,0.34586d-4,0.36083d-4,
     .0.28130d-4,0.28170d-4,0.22030d-4,0.22775d-4,0.23344d-4,0.22467d-4,
     .0.21972d-4,0.22205d-4,0.21760d-4,0.21928d-4,0.21878d-4,0.18086d-4,
     .0.21770d-4,0.22042d-4,0.21889d-4,0.22010d-4,0.22309d-4,0.28848d-4,
     .0.29377d-4,0.30359d-4,0.22742d-4,0.23007d-4,0.22995d-4,0.23004d-4,
     .0.23191d-4,0.18617d-4,0.22989d-4,0.29250d-4,0.29626d-4,0.29784d-4,
     .0.29815d-4,0.29560d-4,0.30549d-4,0.31224d-4,0.24269d-4,0.24463d-4,
     .0.24486d-4,0.31158d-4,0.31135d-4,0.31814d-4,0.31449d-4,0.40492d-4,
     .0.40665d-4,0.32667d-4,0.40804d-4,0.31182d-4,0.39499d-4,0.50193d-4,
     .0.59469d-4,0.49406d-4,0.59765d-4,0.59751d-4,0.49834d-4,0.50134d-4,
     .0.49918d-4,0.39875d-4,0.31689d-4,0.32459d-4,0.25069d-4,0.31248d-4,
     .0.31200d-4,0.31041d-4,0.24673d-4,0.24651d-4,0.24554d-4,0.32080d-4,
     .0.41238d-4,0.41642d-4,0.34722d-4,0.34571d-4,0.34542d-4,0.34706d-4,
     .0.34198d-4,0.26417d-4,0.20194d-4,0.20107d-4,0.20138d-4,0.20127d-4,
     .0.20149d-4,0.20175d-4,0.20250d-4,0.20324d-4,0.25459d-4,0.32117d-4,
     .0.32174d-4,0.32005d-4,0.50746d-4,0.50518d-4,0.51762d-4,0.52418d-4,
     .0.51531d-4,0.50465d-4,0.50907d-4,0.50965d-4,0.40924d-4,0.42046d-4,
     .0.51977d-4,0.51957d-4,0.42105d-4,0.33123d-4,0.33209d-4,0.41172d-4,
     .0.40610d-4,0.33579d-4,0.26792d-4,0.25648d-4,0.26100d-4,0.26333d-4,
     .0.26613d-4,0.21160d-4,0.20699d-4,0.20928d-4,0.21530d-4,0.21880d-4,
     .0.21918d-4,0.26612d-4,0.26526d-4,0.26613d-4/
      DATA M4/0.11515d2,0.11572d2,0.11630d2,0.11688d2,0.11747d2,
     .0.11806d2,0.11865d2,0.11924d2,0.11983d2,0.12043d2,0.12104d2,
     .0.12164d2,0.12225d2,0.12286d2,0.12348d2,0.12409d2,0.12471d2,
     .0.12534d2,0.12596d2,0.12659d2,0.12723d2,0.12786d2,0.12850d2,
     .0.12914d2,0.12979d2,0.13044d2,0.13109d2,0.13175d2,0.13240d2,
     .0.13307d2,0.13373d2,0.13440d2,0.13507d2,0.13575d2,0.13643d2,
     .0.13711d2,0.13779d2,0.13848d2,0.13918d2,0.13987d2,0.14057d2,
     .0.14127d2,0.14198d2,0.14269d2,0.14340d2,0.14412d2,0.14484d2,
     .0.14557d2,0.14629d2,0.14703d2,0.14776d2,0.14850d2,0.14924d2,
     .0.14999d2,0.15074d2,0.15149d2,0.15225d2,0.15301d2,0.15378d2,
     .0.15454d2,0.15532d2,0.15609d2,0.15687d2,0.15766d2,0.15845d2,
     .0.15924d2,0.16003d2,0.16084d2,0.16164d2,0.16245d2,0.16326d2,
     .0.16408d2,0.16490d2,0.16572d2,0.16655d2,0.16738d2,0.16822d2,
     .0.16906d2,0.16991d2,0.17075d2,0.17161d2,0.17247d2,0.17333d2,
     .0.17420d2,0.17507d2,0.17594d2,0.17682d2,0.17771d2,0.17859d2,
     .0.17949d2,0.18038d2,0.18129d2,0.18219d2,0.18310d2,0.18402d2,
     .0.18494d2,0.18586d2,0.18679d2,0.18773d2,0.18867d2,0.18961d2,
     .0.19056d2,0.19151d2,0.19247d2,0.19343d2,0.19440d2,0.19537d2,
     .0.19635d2,0.19733d2,0.19831d2,0.19931d2,0.20030d2,0.20130d2,
     .0.20231d2,0.20332d2,0.20434d2,0.20536d2,0.20639d2,0.20742d2,
     .0.20846d2,0.20950d2,0.21055d2,0.21160d2,0.21266d2,0.21372d2,
     .0.21479d2,0.21586d2,0.21694d2,0.21803d2,0.21912d2,0.22021d2,
     .0.22131d2,0.22242d2,0.22353d2,0.22465d2,0.22577d2,0.22690d2,
     .0.22804d2,0.22918d2,0.23032d2,0.23147d2,0.23263d2,0.23379d2,
     .0.23496d2,0.23614d2,0.23732d2,0.23851d2,0.23970d2,0.24090d2,
     .0.24210d2,0.24331d2,0.24453d2,0.24575d2,0.24698d2,0.24821d2,
     .0.24946d2,0.25070d2,0.25196d2,0.25322d2,0.25448d2,0.25576d2,
     .0.25703d2,0.25832d2,0.25961d2,0.26091d2,0.26221d2,0.26352d2,
     .0.26484d2,0.26617d2,0.26750d2,0.26883d2,0.27018d2,0.27153d2,
     .0.27289d2,0.27425d2,0.27562d2,0.27700d2,0.27839d2,0.27978d2,
     .0.28118d2,0.28258d2,0.28400d2,0.28542d2,0.28684d2,0.28828d2,
     .0.28972d2,0.29117d2,0.29262d2,0.29409d2,0.29556d2,0.29703d2,
     .0.29852d2,0.30001d2,0.30151d2,0.30302d2,0.30453d2,0.30606d2,
     .0.30759d2,0.30913d2,0.31067d2,0.31222d2,0.31379d2,0.31535d2,
     .0.31693d2,0.31852d2,0.32011d2,0.32171d2,0.32332d2,0.32493d2,
     .0.32656d2,0.32819d2,0.32983d2,0.33148d2,0.33314d2,0.33480d2,
     .0.33648d2,0.33816d2,0.33985d2,0.34155d2,0.34326d2,0.34498d2,
     .0.34670d2,0.34843d2/
      DATA L4/0.83880d-4,0.94605d-4,0.93819d-4,0.93227d-4,0.81170d-4,
     .0.68840d-4,0.48266d-4,0.50217d-4,0.41231d-4,0.33997d-4,0.27690d-4,
     .0.23369d-4,0.23469d-4,0.34484d-4,0.34450d-4,0.34263d-4,0.33610d-4,
     .0.33480d-4,0.33364d-4,0.33275d-4,0.33379d-4,0.23668d-4,0.23691d-4,
     .0.28230d-4,0.28351d-4,0.48646d-4,0.48658d-4,0.48154d-4,0.46850d-4,
     .0.53693d-4,0.53753d-4,0.54571d-4,0.54752d-4,0.32813d-4,0.33043d-4,
     .0.39978d-4,0.58723d-4,0.48397d-4,0.48656d-4,0.47628d-4,0.57537d-4,
     .0.67123d-4,0.67406d-4,0.56695d-4,0.39503d-4,0.55958d-4,0.65635d-4,
     .0.64350d-4,0.45314d-4,0.53731d-4,0.53308d-4,0.53774d-4,0.53566d-4,
     .0.45511d-4,0.52461d-4,0.61298d-4,0.82977d-4,0.92286d-4,0.10318d-3,
     .0.12520d-3,0.13672d-3,0.13651d-3,0.12520d-3,0.11320d-3,0.11275d-3,
     .0.10034d-3,0.79778d-4,0.78615d-4,0.67693d-4,0.57839d-4,0.48444d-4,
     .0.48837d-4,0.34772d-4,0.29875d-4,0.35155d-4,0.30154d-4,0.41169d-4,
     .0.57163d-4,0.67420d-4,0.97809d-4,0.10829d-3,0.10791d-3,0.96529d-4,
     .0.10717d-3,0.86701d-4,0.65595d-4,0.63740d-4,0.38423d-4,0.33341d-4,
     .0.32969d-4,0.32738d-4,0.28374d-4,0.24738d-4,0.27686d-4,0.27843d-4,
     .0.31992d-4,0.37116d-4,0.43340d-4,0.57877d-4,0.57606d-4,0.66159d-4,
     .0.57095d-4,0.56842d-4,0.64834d-4,0.73829d-4,0.64325d-4,0.48277d-4,
     .0.63054d-4,0.70608d-4,0.79519d-4,0.89004d-4,0.98700d-4,0.96459d-4,
     .0.95036d-4,0.95088d-4,0.85697d-4,0.66356d-4,0.66835d-4,0.49918d-4,
     .0.42995d-4,0.29648d-4,0.37809d-4,0.37636d-4,0.49208d-4,0.49251d-4,
     .0.55782d-4,0.64117d-4,0.63798d-4,0.71640d-4,0.71084d-4,0.71333d-4,
     .0.60626d-4,0.67690d-4,0.45206d-4,0.39911d-4,0.45301d-4,0.39276d-4,
     .0.35100d-4,0.45100d-4,0.35351d-4,0.31648d-4,0.35414d-4,0.39215d-4,
     .0.34440d-4,0.47640d-4,0.62094d-4,0.53848d-4,0.60163d-4,0.59802d-4,
     .0.51271d-4,0.58858d-4,0.58773d-4,0.46099d-4,0.35702d-4,0.32144d-4,
     .0.35878d-4,0.44249d-4,0.56363d-4,0.44622d-4,0.56788d-4,0.56553d-4,
     .0.57946d-4,0.57664d-4,0.81611d-4,0.81054d-4,0.79351d-4,0.71357d-4,
     .0.49372d-4,0.49191d-4,0.53819d-4,0.59673d-4,0.47043d-4,0.42957d-4,
     .0.62795d-4,0.77766d-4,0.94592d-4,0.11317d-3,0.11328d-3,0.13244d-3,
     .0.11187d-3,0.12173d-3,0.84237d-4,0.82440d-4,0.66417d-4,0.52862d-4,
     .0.66841d-4,0.60313d-4,0.53360d-4,0.47024d-4,0.59581d-4,0.60566d-4,
     .0.68704d-4,0.12552d-3,0.10492d-3,0.95251d-4,0.11324d-3,0.11354d-3,
     .0.10459d-3,0.11422d-3,0.11401d-3,0.76329d-4,0.77161d-4,0.70855d-4,
     .0.91092d-4,0.83645d-4,0.93542d-4,0.77241d-4,0.76382d-4,0.76732d-4,
     .0.68324d-4,0.77708d-4,0.67232d-4,0.67348d-4,0.60600d-4,0.53897d-4,
     .0.47869d-4,0.47879d-4,0.36148d-4,0.31325d-4,0.22097d-4,0.22701d-4,
     .0.20439d-4,0.21101d-4/

*   CMS 2311.00130, fig. 7: Z -> ll + A -> gamgam  and 1 < MA < 30 GeV
      DATA M5/1d0 ,2d0 ,3d0 ,4d0 ,5d0 ,6d0 ,7d0 ,8d0 ,9d0 ,10d0,11d0,
     .12d0,13d0,14d0,15d0,16d0,17d0,18d0,19d0,20d0,21d0,22d0,23d0,24d0,
     .25d0,26d0,27d0,28d0,29d0,30d0/
      DATA L5/17.89d-3,7.948d-3,7.064d-3,3.880d-3,2.879d-3,1.136d-3,
     .3.997d-3,5.003d-3,2.044d-3,2.212d-3,2.581d-3,2.151d-3,3.748d-3,
     .3.137d-3,2.858d-3,5.702d-3,5.486d-3,4.183d-3,5.569d-3,4.625d-3,
     .4.226d-3,5.445d-3,5.457d-3,3.956d-3,4.144d-3,5.355d-3,4.550d-3,
     .9.954d-3,6.717d-3,6.875d-3/


      PROB(75)=0d0
      PROB(76)=0d0
      gg13=48.61d0
      VB13=3.766d0
      VH13=2.238d0
      TH13=.7029d0

*   ATLAS 2004.01678, tab. 1: A -> had and 0.5 < MA < 4 GeV

      S=0d0
      DO I=1,3
       IF(SMASS(I).GE.MHMIN .AND. SMASS(I).LE.MHMAX)THEN
        S=S+(CJ(I)**2*gg13+CV(I)**2*(VB13+VH13)+CU(I)**2*TH13)
     .    *BRHAZ(I,1)*BRJJ(4)
       ENDIF
      ENDDO

      I=1
      DOWHILE(M1(I).LE.PMASS(1) .AND. I.LT.N1)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. PMASS(1).LE.M1(N1))THEN
       SMAX=L1(I-1)+(PMASS(1)-M1(I-1))/(M1(I)-M1(I-1))*(L1(I)-L1(I-1))
       PROB(75)=PROB(75)+DDIM(S/SMAX,1d0)
      ENDIF

*   ATLAS 2110.13673, fig. 17b: Z -> ll + A -> mumu  and 15 < MA < 30 GeV

      S=0d0
      DO I=1,3
       IF(SMASS(I).GE.MHMIN .AND. SMASS(I).LE.MHMAX)THEN
        S=S+CJ(I)**2*gg13*BRHAZ(I,1)*BRMM(4)*.673d0
       ENDIF
      ENDDO

      I=1
      DOWHILE(M2(I).LE.PMASS(1) .AND. I.LT.N2)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. PMASS(1).LE.M2(N2))THEN
       SMAX=L2(I-1)+(PMASS(1)-M2(I-1))/(M2(I)-M2(I-1))*(L2(I)-L2(I-1))
       PROB(75)=PROB(75)+DDIM(S/SMAX,1d0)
      ENDIF

*   CMS 2111.01299, fig. 4a: A -> mumu  and 4 < MA < 35 GeV

      S=0d0
      DO I=1,3
       IF(SMASS(I).GE.MHMIN .AND. SMASS(I).LE.MHMAX)THEN
        S=S+CJ(I)**2*BRHAZ(I,1)*BRMM(4)
       ENDIF
      ENDDO

      I=1
      DOWHILE(M3(I).LE.PMASS(1) .AND. I.LT.N3)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. PMASS(1).LE.M3(N3))THEN
       SMAX=L3(I-1)+(PMASS(1)-M3(I-1))/(M3(I)-M3(I-1))*(L3(I)-L3(I-1))
       PROB(76)=PROB(76)+DDIM(S/SMAX,1d0)
      ENDIF

      I=1
      DOWHILE(M4(I).LE.PMASS(1) .AND. I.LT.N4)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. PMASS(1).LE.M4(N4))THEN
       SMAX=L4(I-1)+(PMASS(1)-M4(I-1))/(M4(I)-M4(I-1))*(L4(I)-L4(I-1))
       PROB(76)=PROB(76)+DDIM(S/SMAX,1d0)
      ENDIF

*   CMS 2311.00130, fig. 7: Z -> ll + A -> gamgam  and 1 < MA < 30 GeV

      S=0d0
      DO I=1,2
       IF(SMASS(I).GE.MHMIN .AND. SMASS(I).LE.MHMAX)THEN
        S=S+CJ(I)**2*gg13*BRHAZ(I,1)*BRGG(4)*.673d0
       ENDIF
      ENDDO

      I=1
      DOWHILE(M5(I).LE.PMASS(1) .AND. I.LT.N5)
       I=I+1
      ENDDO
      IF(I.GT.1 .AND. PMASS(1).LE.M5(N5))THEN
       SMAX=L5(I-1)+(PMASS(1)-M5(I-1))/(M5(I)-M5(I-1))*(L5(I)-L5(I-1))
       PROB(76)=PROB(76)+DDIM(S/SMAX,1d0)
      ENDIF

      END


      SUBROUTINE HEAVYA_ZHSM(PROB)

*     A -> ZHSM, HSM -> bb from CMS: 1807.02826 + 1903.00941
*     A -> ZHSM, HSM -> tautau from CMS: 1510.01181 + 1910.11634
*     PROB(78) =/= 0: excluded
*     A -> ZHSM, HSM -> bb from ATLAS: 2207.00230
*     A -> ZHSM, HSM -> tautau from ATLAS: 1502.04478
*     PROB(79) =/= 0: excluded

      IMPLICIT NONE

      INTEGER I,J,N1,N2,N4,N5,N7,N8,N9
      PARAMETER (N1=51,N2=53,N4=58,N5=10,N7=14,N8=9,N9=18)

      DOUBLE PRECISION PROB(*),M1(N1),L1(N1),M2(N2),L2(N2)
      DOUBLE PRECISION M4(N4),L4(N4),L4B(N4),M5(N5),L5(N5)
      DOUBLE PRECISION M7(N7),L7(N7),M8(N8),L8(N8)
      DOUBLE PRECISION M9(N9),L9(N9),L9B(N9)
      DOUBLE PRECISION MA,SIGJBB,SIGBBB,SIGJLL,SIGBLL,L,ggF13,bbH13
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION MHmin,MHmax

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/HIGGSFIT/MHmin,MHmax

*   CMS 1903.00941 fig 5a (ggF -> A -> ZHSM, HSM -> bb)
      DATA M1/2.25d2,2.35d2,2.40d2,2.45d2,2.50d2,2.55d2,2.65d2,2.70d2,
     .2.75d2,2.80d2,2.85d2,2.95d2,3.00d2,3.05d2,3.10d2,3.15d2,3.20d2,
     .3.30d2,3.35d2,3.40d2,3.50d2,3.60d2,3.80d2,3.90d2,4.00d2,4.10d2,
     .4.30d2,4.40d2,4.60d2,4.70d2,4.90d2,5.00d2,5.25d2,5.50d2,5.75d2,
     .6.00d2,6.25d2,6.50d2,6.75d2,7.00d2,7.25d2,7.50d2,7.75d2,8.00d2,
     .8.25d2,8.50d2,8.75d2,9.00d2,9.50d2,9.75d2,1.00d3/
      DATA L1/7.8694d-1,3.5916d-1,3.8894d-1,7.0770d-1,8.2091d-1,
     .6.2560d-1,6.7967d-1,6.7240d-1,6.3711d-1,7.1989d-1,6.9598d-1,
     .4.5546d-1,3.3952d-1,2.7001d-1,2.2829d-1,2.1499d-1,2.2683d-1,
     .3.7358d-1,3.9508d-1,3.7865d-1,2.0652d-1,1.5039d-1,2.4399d-1,
     .2.3884d-1,2.0291d-1,1.6002d-1,8.5323d-2,7.0728d-2,6.0887d-2,
     .6.2368d-2,7.9997d-2,7.2467d-2,7.9594d-2,4.1929d-2,3.2668d-2,
     .3.4909d-2,2.9929d-2,2.8912d-2,2.9880d-2,3.9561d-2,4.2100d-2,
     .3.9169d-2,2.5838d-2,1.7032d-2,1.1388d-2,9.8492d-3,1.1023d-2,
     .1.2719d-2,1.7473d-2,1.9602d-2,2.0650d-2/

*   CMS 1903.00941 fig 5b (bbA -> ZHSM, HSM -> bb)
      DATA M2/2.25d2,2.35d2,2.40d2,2.45d2,2.50d2,2.55d2,2.60d2,2.65d2,
     .2.70d2,2.75d2,2.80d2,2.85d2,3.00d2,3.10d2,3.15d2,3.20d2,3.25d2,
     .3.30d2,3.35d2,3.40d2,3.50d2,3.60d2,3.70d2,3.80d2,3.90d2,4.00d2,
     .4.10d2,4.20d2,4.30d2,4.40d2,4.50d2,4.60d2,4.70d2,4.90d2,5.00d2,
     .5.25d2,5.50d2,5.75d2,6.00d2,6.25d2,6.50d2,6.75d2,7.00d2,7.25d2,
     .7.50d2,8.00d2,8.25d2,8.50d2,9.00d2,9.25d2,9.50d2,9.75d2,1.00d3/
      DATA L2/8.8507d-1,4.1274d-1,4.3231d-1,6.4800d-1,8.0138d-1,
     .6.5623d-1,6.8894d-1,7.0975d-1,6.7071d-1,6.3655d-1,7.4587d-1,
     .7.4587d-1,3.7741d-1,2.9134d-1,2.8881d-1,3.1365d-1,3.7972d-1,
     .4.8589d-1,4.9955d-1,4.6142d-1,2.4028d-1,2.2282d-1,2.9598d-1,
     .3.6404d-1,3.4935d-1,3.1142d-1,2.5437d-1,1.7922d-1,1.1937d-1,
     .9.1615d-2,8.0745d-2,7.3751d-2,7.3144d-2,8.9605d-2,8.0537d-2,
     .9.1365d-2,4.8903d-2,3.9779d-2,4.2093d-2,3.6803d-2,3.5535d-2,
     .3.8201d-2,5.0418d-2,5.3468d-2,4.9883d-2,2.3526d-2,1.5182d-2,
     .1.2985d-2,1.7649d-2,2.1813d-2,2.5328d-2,2.7331d-2,2.8566d-2/

*   ATLAS 2207.00230
      DATA M4/220d0,260d0,300d0,340d0,350d0,360d0,370d0,380d0,390d0,
     .400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,490d0,500d0,520d0,
     .540d0,560d0,580d0,600d0,620d0,640d0,660d0,680d0,700d0,720d0,740d0,
     .760d0,780d0,800d0,850d0,900d0,950d0,1000d0,1050d0,1100d0,1150d0,
     .1200d0,1250d0,1300d0,1350d0,1400d0,1450d0,1500d0,1550d0,1600d0,
     .1650d0,1700d0,1750d0,1800d0,1850d0,1900d0,1950d0,2000d0/
*   fig 9a (ggF -> A -> ZHSM) assuming BR(HSM -> bb)=.569
      DATA L4/0.5536d0,0.8408d0,0.5999d0,0.2456d0,0.2435d0,0.2281d0,
     .0.1664d0,0.1598d0,0.1320d0,0.1053d0,0.7299d-1,0.6626d-1,0.5908d-1,
     .0.6331d-1,0.6438d-1,0.7154d-1,0.7734d-1,0.9876d-1,0.1055d0,
     .0.8630d-1,0.5875d-1,0.4374d-1,0.3403d-1,0.3154d-1,0.3173d-1,
     .0.2381d-1,0.2188d-1,0.1895d-1,0.1368d-1,0.1363d-1,0.1488d-1,
     .0.1284d-1,0.1152d-1,0.1191d-1,0.1539d-1,0.1032d-1,0.8630d-2,
     .0.6776d-2,0.7330d-2,0.6611d-2,0.4591d-2,0.3539d-2,0.4030d-2,
     .0.3873d-2,0.3402d-2,0.3085d-2,0.2936d-2,0.2839d-2,0.2584d-2,
     .0.2403d-2,0.2281d-2,0.2263d-2,0.2244d-2,0.2228d-2,0.2272d-2,
     .0.2472d-2,0.2771d-2,0.2934d-2/
*   fig 9b (bbA -> ZHSM) assuming BR(HSM -> bb)=.569
      DATA L4B/0.5154d0,0.4396d0,0.3436d0,0.3184d0,0.3421d0,0.3124d0,
     .0.2367d0,0.2308d0,0.1606d0,0.1323d0,0.1320d0,0.1142d0,0.9966d-1,
     .0.9374d-1,0.8097d-1,0.7279d-1,0.7161d-1,0.6539d-1,0.9450d-1,
     .0.9841d-1,0.7323d-1,0.4778d-1,0.3919d-1,0.3799d-1,0.3653d-1,
     .0.3244d-1,0.2718d-1,0.2097d-1,0.1756d-1,0.1356d-1,0.1208d-1,
     .0.1290d-1,0.1372d-1,0.1208d-1,0.1560d-1,0.9661d-2,0.7789d-2,
     .0.5714d-2,0.5780d-2,0.5282d-2,0.4281d-2,0.3906d-2,0.3427d-2,
     .0.3014d-2,0.3028d-2,0.3162d-2,0.3209d-2,0.3556d-2,0.3283d-2,
     .0.2676d-2,0.2677d-2,0.2483d-2,0.2410d-2,0.2403d-2,0.2381d-2,
     .0.2407d-2,0.2458d-2,0.2421d-2/

*   ATLAS 1502.04478 fig 3a (ggF -> A -> ZHSM, HSM -> tautau)
      DATA M5/220d0,240d0,260d0,300d0,340d0,350d0,400d0,500d0,800d0,
     .1000d0/
      DATA L5/0.9890d-1,0.7876d-1,0.1797d-1,0.6652d-1,0.3920d-1,
     .0.3835d-1,0.4281d-1,0.1773d-1,0.1361d-1,0.1312d-1/

*   CMS 1510.01181 fig 10a (ggF -> A -> ZHSM, Z -> ll and HSM -> tautau)
      DATA M7/220d0,230d0,240d0,250d0,260d0,270d0,280d0,290d0,300d0,
     .310d0,320d0,330d0,340d0,350d0/
      DATA L7/0.1670d-1,0.1446d-1,0.1246d-1,0.9710d-2,0.7174d-2,
     .0.5761d-2,0.5036d-2,0.4638d-2,0.4891d-2,0.5145d-2,0.5942d-2,
     .0.6920d-2,0.6993d-2,0.6775d-2/

*   CMS 1910.11634 fig 5 (pp -> A -> ZHSM, Z -> ee,mm and HSM -> tautau)
      DATA M8/220d0,240d0,260d0,280d0,300d0,320d0,340d0,350d0,400d0/
      DATA L8/26.573d-3,23.461d-3,11.469d-3,5.701d-3,4.0395d-3,
     .4.6992d-3,6.8103d-3,6.262d-3,4.4301d-3/

*   CMS 1807.02826
      DATA M9/800d0,850d0,900d0,950d0,1000d0,1050d0,1100d0,1150d0,
     .1200d0,1250d0,1300d0,1400d0,1500d0,1600d0,1650d0,1700d0,1800d0,
     .2000d0/
*   fig 8 (ggF -> A -> ZHSM, HSM -> bb)
      DATA L9/0.2926d-1,0.3174d-1,0.2599d-1,0.1594d-1,0.1134d-1,
     .0.8068d-2,0.7273d-2,0.8752d-2,0.9082d-2,0.9286d-2,0.1030d-1,
     .0.1061d-1,0.8310d-2,0.6000d-2,0.5739d-2,0.6000d-2,0.7382d-2,
     .0.6855d-2/
*   fig 8 (bbA -> ZHSM, HSM -> bb)
      DATA L9B/0.2638d-1,0.3081d-1,0.2618d-1,0.1667d-1,0.1186d-1,
     .0.8188d-2,0.7382d-2,0.9150d-2,0.9926d-2,0.1007d-1,0.1126d-1,
     .0.1126d-1,0.9708d-2,0.6557d-2,0.6180d-2,0.6509d-2,0.7949d-2,
     .0.7604d-2/

      PROB(78)=0d0
      PROB(79)=0d0

      DO I=1,2
       MA=PMASS(I)
       SIGJBB=0d0
       SIGBBB=0d0
       SIGJLL=0d0
       SIGBLL=0d0

       DO J=1,3
        IF(SMASS(J).GE.MHMIN .AND. SMASS(J).LE.MHMAX)THEN
         SIGJBB=SIGJBB+CJ(I+3)**2*ggF13(MA)*BRAHZ(I,J)*BRBB(J)
         SIGBBB=SIGBBB+CB(I+3)**2*bbH13(MA)*BRAHZ(I,J)*BRBB(J)
         SIGJLL=SIGJLL+CJ(I+3)**2*ggF13(MA)*BRAHZ(I,J)*BRLL(J)
         SIGBLL=SIGBLL+CB(I+3)**2*bbH13(MA)*BRAHZ(I,J)*BRLL(J)
        ENDIF
       ENDDO

       J=1
       DOWHILE(M1(J).LE.MA .AND. J.LT.N1)
	 J=J+1
       ENDDO
       IF(J.GE.2 .AND. MA.LE.M1(N1)) THEN
	L=L1(J-1)+(MA-M1(J-1))/(M1(J)-M1(J-1))*(L1(J)-L1(J-1))
	PROB(78)=PROB(78)+DDIM(SIGJBB/L,1d0)
       ENDIF

       J=1
       DOWHILE(M2(J).LE.MA .AND. J.LT.N2)
	 J=J+1
       ENDDO
       IF(J.GE.2 .AND. MA.LE.M2(N2)) THEN
	L=L2(J-1)+(MA-M2(J-1))/(M2(J)-M2(J-1))*(L2(J)-L2(J-1))
	PROB(78)=PROB(78)+DDIM(SIGBBB/L,1d0)
       ENDIF

       J=1
       DOWHILE(M4(J).LE.MA .AND. J.LT.N4)
	 J=J+1
       ENDDO
       IF(J.GE.2 .AND. MA.LE.M4(N4)) THEN
	L=L4(J-1)+(MA-M4(J-1))/(M4(J)-M4(J-1))*(L4(J)-L4(J-1))
	PROB(79)=PROB(79)+DDIM(SIGJBB/.569d0/L,1d0)
	L=L4B(J-1)+(MA-M4(J-1))/(M4(J)-M4(J-1))*(L4B(J)-L4B(J-1))
	PROB(79)=PROB(79)+DDIM(SIGBBB/.569d0/L,1d0)
       ENDIF

       J=1
       DOWHILE(M5(J).LE.MA .AND. J.LT.N5)
	 J=J+1
       ENDDO
       IF(J.GE.2 .AND. MA.LE.M5(N5)) THEN
	L=L5(J-1)+(MA-M5(J-1))/(M5(J)-M5(J-1))*(L5(J)-L5(J-1))
	PROB(79)=PROB(79)+DDIM(SIGJLL/L,1d0)
       ENDIF

       J=1
       DOWHILE(M7(J).LE.MA .AND. J.LT.N7)
	 J=J+1
       ENDDO
       IF(J.GE.2 .AND. MA.LE.M7(N7)) THEN
	L=L7(J-1)+(MA-M7(J-1))/(M7(J)-M7(J-1))*(L7(J)-L7(J-1))
	PROB(78)=PROB(78)+DDIM(SIGJLL*.101d0/L,1d0)
       ENDIF

       J=1
       DOWHILE(M8(J).LE.MA .AND. J.LT.N8)
	 J=J+1
       ENDDO
       IF(J.GE.2 .AND. MA.LE.M8(N8)) THEN
	L=L8(J-1)+(MA-M8(J-1))/(M8(J)-M8(J-1))*(L8(J)-L8(J-1))
	PROB(78)=PROB(78)+DDIM((SIGJLL+.76d0*SIGBLL)*.673d0/L,1d0)
       ENDIF

       J=1
       DOWHILE(M9(J).LE.MA .AND. J.LT.N9)
	 J=J+1
       ENDDO
       IF(J.GE.2 .AND. MA.LE.M9(N9)) THEN
	L=L9(J-1)+(MA-M9(J-1))/(M9(J)-M9(J-1))*(L9(J)-L9(J-1))
	PROB(78)=PROB(78)+DDIM(SIGJBB/L,1d0)
	L=L9B(J-1)+(MA-M9(J-1))/(M9(J)-M9(J-1))*(L9B(J)-L9B(J-1))
	PROB(78)=PROB(78)+DDIM(SIGBBB/L,1d0)
       ENDIF

      ENDDO

      END


      SUBROUTINE HEAVYHA_ZAH(PROB)

*     H/A -> ZA/H ->bb from CMS arXiv:1911.03781 and
*     PROB(80) =/= 0: excluded
*     H/A -> ZA/H ->bb/WW from ATLAS arXiv:2011.05639
*     PROB(81) =/= 0: excluded

      IMPLICIT NONE

      INTEGER I,J,K,NHZH1,NHZH2,NHZH3,NHZH4,NHZHT1,NHZHT2
      INTEGER NHZHT3,NHZHT4,HZHT1(10000,3),HZHT2(10000,3)
      INTEGER HZHT3(10000,3),HZHT4(10000,3),HHHT(10000,3)

      DOUBLE PRECISION PROB(*),MA,MH,SIG,L,ggF13,bbH13
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5),CB(5),CL(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5),HZH1(10000,3)
      DOUBLE PRECISION HZH2(10000,3),HZH3(10000,3),HZH4(10000,3)

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/CB/CB,CL
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/LHZH/HZH1,HZH2,HZH3,HZH4,HZHT1,HZHT2,HZHT3,HZHT4,
     .      NHZH1,NHZH2,NHZH3,NHZH4,NHZHT1,NHZHT2,NHZHT3,NHZHT4

      PROB(80)=0d0
      PROB(81)=0d0

* CMS

      DO I=1,2
       MA=PMASS(I)
       DO J=1,3
        MH=SMASS(J)
        SIG=(CJ(I+3)**2*ggF13(MA)+CB(I+3)**2*bbH13(MA))
     .      *BRAHZ(I,J)*BRBB(J)*0.0673d0
        IF(SIG.NE.0d0)THEN
	 CALL TABTRI(HZH1,NHZH1,HZHT1,NHZHT1,MA,MH,L)
         IF(L.NE.0d0)PROB(80)=PROB(80)+DDIM(SIG/L,1d0)
        ENDIF
       ENDDO
      ENDDO

      DO I=1,3
       MH=SMASS(I)
       DO J=1,2
        MA=PMASS(J)
        SIG=(CJ(I)**2*ggF13(MH)+CB(I)**2*bbH13(MH))
     .      *BRHAZ(I,J)*BRBB(J+3)*0.0673d0
        IF(SIG.NE.0d0)THEN
	 CALL TABTRI(HZH1,NHZH1,HZHT1,NHZHT1,MH,MA,L)
         IF(L.NE.0d0)PROB(80)=PROB(80)+DDIM(SIG/L,1d0)
        ENDIF
       ENDDO
      ENDDO

* ATLAS ggF/bb

      DO I=1,2
       MA=PMASS(I)
       DO J=1,3
        MH=SMASS(J)
        SIG=CJ(I+3)**2*ggF13(MA)*BRAHZ(I,J)*BRBB(J)
        IF(SIG.NE.0d0)THEN
	 CALL TABTRI(HZH2,NHZH2,HZHT2,NHZHT2,MA,MH,L)
         IF(L.NE.0d0)PROB(81)=PROB(81)+DDIM(SIG/L,1d0)
        ENDIF
       ENDDO
      ENDDO

      DO I=1,3
       MH=SMASS(I)
       DO J=1,2
        MA=PMASS(J)
        SIG=CJ(I)**2*ggF13(MH)*BRHAZ(I,J)*BRBB(J+3)
        IF(SIG.NE.0d0)THEN
	 CALL TABTRI(HZH2,NHZH2,HZHT2,NHZHT2,MH,MA,L)
         IF(L.NE.0d0)PROB(81)=PROB(81)+DDIM(SIG/L,1d0)
        ENDIF
       ENDDO
      ENDDO

* ATLAS bbH/bb

      DO I=1,2
       MA=PMASS(I)
       DO J=1,3
        MH=SMASS(J)
        SIG=CB(I+3)**2*bbH13(MA)*BRAHZ(I,J)*BRBB(J)
        IF(SIG.NE.0d0)THEN
	 CALL TABTRI(HZH3,NHZH3,HZHT3,NHZHT3,MA,MH,L)
         IF(L.NE.0d0)PROB(81)=PROB(81)+DDIM(SIG/L,1d0)
        ENDIF
       ENDDO
      ENDDO

      DO I=1,3
       MH=SMASS(I)
       DO J=1,2
        MA=PMASS(J)
        SIG=CB(I)**2*bbH13(MH)*BRHAZ(I,J)*BRBB(J+3)
        IF(SIG.NE.0d0)THEN
	 CALL TABTRI(HZH3,NHZH3,HZHT3,NHZHT3,MH,MA,L)
         IF(L.NE.0d0)PROB(81)=PROB(81)+DDIM(SIG/L,1d0)
        ENDIF
       ENDDO
      ENDDO

* ATLAS ggF/WW

      DO I=1,2
       MA=PMASS(I)
       DO J=1,3
        MH=SMASS(J)
        SIG=CJ(I+3)**2*ggF13(MA)*BRAHZ(I,J)*BRWW(J)
        IF(SIG.NE.0d0)THEN
	 CALL TABTRI(HZH4,NHZH4,HZHT4,NHZHT4,MA,MH,L)
         IF(L.NE.0d0)PROB(81)=PROB(81)+DDIM(SIG/L,1d0)
        ENDIF
       ENDDO
      ENDDO

      END


      SUBROUTINE HEAVYHA_HSMHA(PROB)

*     H/A -> (HSM -> tautau) + (H/A ->bb) from CMS 2106.10361
*     PROB(82) =/= 0: excluded
*     H/A -> (HSM -> bb) + (H/A ->bb) from CMS 2204.12413
*     PROB(83) =/= 0: excluded
*     H/A -> (HSM -> gamgam) + (H/A ->bb) from CMS 2310.01643
*     PROB(74) =/= 0: excluded

      IMPLICIT NONE

      INTEGER I,J,K,N,N1,N2,N3,N4,NHHH,NHHH2,NHHH3,NHHHT,NHHHT2,NHHHT3
      INTEGER HHHT(10000,3),HHHT2(10000,3),HHHT3(10000,3)

      DOUBLE PRECISION PROB(*),MA,MH,MS,SIG,L,ggF13
      DOUBLE PRECISION HHH(10000,3),HHH2(10000,3),HHH3(10000,3)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION MHmin,MHmax

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/HIGGSFIT/MHmin,MHmax
      COMMON/LHHH/HHH,HHH2,HHH3,HHHT,HHHT2,HHHT3,
     .      NHHH,NHHH2,NHHH3,NHHHT,NHHHT2,NHHHT3

      PROB(82)=0d0
      PROB(83)=0d0
      PROB(74)=0d0

      MA=PMASS(2)
      DO J=1,3
       MH=SMASS(J)
       IF(MH.GE.MHMIN .AND. MH.LE.MHMAX)THEN
        MS=PMASS(1)
        SIG=CJ(5)**2*ggF13(MA)*BRAHA(J)*BRBB(4)
        IF(SIG.NE.0d0)THEN
         CALL TABTRI(HHH,NHHH,HHHT,NHHHT,MA,MS,L)
         IF(L.NE.0d0)PROB(82)=PROB(82)+DDIM(SIG*BRLL(J)/L,1d0)
         CALL TABTRI(HHH2,NHHH2,HHHT2,NHHHT2,MA,MS,L)
         IF(L.NE.0d0)PROB(83)=PROB(83)+DDIM(SIG*BRBB(J)/L,1d0)
         CALL TABTRI(HHH3,NHHH3,HHHT3,NHHHT3,MA,MS,L)
         IF(L.NE.0d0)PROB(74)=PROB(74)+DDIM(SIG*BRGG(J)/L,1d0)
        ENDIF
       ENDIF
      ENDDO

      DO I=2,3
      MA=SMASS(I)
       DO J=1,I-1
        MH=SMASS(J)
        IF(MH.GE.MHMIN .AND. MH.LE.MHMAX)THEN
         DO K=1,I-1
          MS=SMASS(K)
          SIG=CJ(I)**2*ggF13(MA)*BRHHH(I+J+K-3)*BRBB(K)
          IF(SIG.NE.0d0)THEN
           CALL TABTRI(HHH,NHHH,HHHT,NHHHT,MA,MS,L)
           IF(L.NE.0d0)PROB(82)=PROB(82)+DDIM(SIG*BRLL(J)/L,1d0)
           CALL TABTRI(HHH2,NHHH2,HHHT2,NHHHT2,MA,MS,L)
           IF(L.NE.0d0)PROB(83)=PROB(83)+DDIM(SIG*BRBB(J)/L,1d0)
           CALL TABTRI(HHH3,NHHH3,HHHT3,NHHHT3,MA,MS,L)
           IF(L.NE.0d0)PROB(74)=PROB(74)+DDIM(SIG*BRGG(J)/L,1d0)
          ENDIF
         ENDDO
        ENDIF
       ENDDO
      ENDDO

      END


      SUBROUTINE HEAVYH_HSMHSM(PROB)

*     VBF -> H -> HSM+HSM -> 4b from ATLAS 2001.05178
*     PROB(84) =/= 0: excluded 
*     ggF -> H -> HSM+HSM from ATLAS 2311.15956
*     PROB(85) =/= 0: excluded 
*     H -> HSM+HSM from CMS 2403.16926
*     PROB(86) =/= 0: excluded 

      IMPLICIT NONE

      INTEGER I,J,K,N1,N2,N3
      PARAMETER (N1=10,N2=41,N3=32)

      DOUBLE PRECISION PROB(*),MH,MHSM,SIG,L,ggF13,VBF13
      DOUBLE PRECISION M1(N1),L1(N1),M2(N2),L2(N2),M3(N3),L3(N3)
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION CU(5),CD(5),CV(3),CJ(5),CG(5)
      DOUBLE PRECISION BRJJ(5),BREE(5),BRMM(5),BRLL(5)
      DOUBLE PRECISION BRCC(5),BRBB(5),BRTT(5),BRWW(3),BRZZ(3)
      DOUBLE PRECISION BRGG(5),BRZG(5),BRHHH(4),BRHAA(3,3)
      DOUBLE PRECISION BRHCHC(3),BRHAZ(3,2),BRAHA(3),BRAHZ(2,3)
      DOUBLE PRECISION BRHCW(5),BRHIGGS(5),BRNEU(5,5,5),BRCHA(5,3)
      DOUBLE PRECISION BRHSQ(3,10),BRHSL(3,7),BRASQ(2,2),BRASL(2)
      DOUBLE PRECISION BRSUSY(5),WIDTH(5)
      DOUBLE PRECISION MHmin,MHmax

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/REDCOUP/CU,CD,CV,CJ,CG
      COMMON/BRN/BRJJ,BREE,BRMM,BRLL,BRCC,BRBB,BRTT,BRWW,
     .      BRZZ,BRGG,BRZG,BRHHH,BRHAA,BRHCHC,BRHAZ,BRAHA,BRAHZ,
     .      BRHCW,BRHIGGS,BRNEU,BRCHA,BRHSQ,BRHSL,BRASQ,BRASL,
     .      BRSUSY,WIDTH
      COMMON/HIGGSFIT/MHmin,MHmax

*  ATLAS: 2001.05178, fig 5a (VBF -> H -> HSM+HSM -> 4b)
      DATA M1/260d0,280d0,300d0,400d0,500d0,600d0,700d0,800d0,900d0,
     .1000d0/
      DATA L1/1090.5d-3,1327.6d-3,905.3d-3,142.8d-3,51.8d-3,31.2d-3,
     .13.2d-3,13.9d-3,4.3d-3,3.7d-3/

*  ATLAS: 2311.15956, fig 1b (ggF -> H -> HSM+HSM)
      DATA M2/251d0,260d0,270d0,280d0,290d0,300d0,312.5d0,325d0,337.5d0,
     .350d0,375d0,400d0,425d0,450d0,475d0,500d0,550d0,600d0,650d0,700d0,
     .750d0,800d0,850d0,900d0,1000d0,1100d0,1200d0,1300d0,1400d0,1500d0,
     .1600d0,1800d0,2000d0,2250d0,2500d0,2750d0,3000d0,3500d0,4000d0,
     .4500d0,5000d0/
      DATA L2/414.1724d-3,580.88711d-3,595.33678d-3,245.68665d-3,
     .243.61684d-3,284.32297d-3,404.23453d-3,173.30747d-3,262.1178d-3,
     .148.09078d-3,139.20797d-3,57.89657d-3,182.21464d-3,40.92366d-3,
     .142.85358d-3,40.80219d-3,23.08846d-3,12.84065d-3,149.29492d-3,
     .15.06086d-3,44.48321d-3,14.14199d-3,86.74655d-3,16.26264d-3,
     .10.19933d-3,14.74682d-3,10.25159d-3,6.08505d-3,7.95851d-3,
     .7.39896d-3,5.15856d-3,3.26557d-3,2.80147d-3,2.82465d-3,
     .1.97049d-3,0.96079d-3,1.14611d-3,1.08775d-3,2.70823d-3,
     .2.04452d-3,1.44471d-3/

*  CMS: 2403.16926, fig 28a (ggF -> H -> HSM+HSM)
      DATA M3/260d0,270d0,280d0,300d0,320d0,350d0,400d0,450d0,500d0,
     .550d0,600d0,700d0,800d0,900d0,1000d0,1100d0,1200d0,1300d0,1400d0,
     .1500d0,1600d0,1700d0,1800d0,1900d0,2000d0,2200d0,2400d0,2600d0,
     .2800d0,3000d0,3500d0,4000d0/
      DATA L3/288.27d-3,247.81d-3,197.55d-3,153.13d-3,119.57d-3,
     .242.96d-3,87.442d-3,48.115d-3,42.964d-3,44.272d-3,30.214d-3,
     .29.797d-3,14.207d-3,10.512d-3,2.7068d-3,1.4436d-3,1.208d-3,
     .1.1441d-3,1.0884d-3,0.76992d-3,0.56244d-3,0.55468d-3,0.53743d-3,
     .0.52064d-3,0.54115d-3,0.56651d-3,0.49214d-3,0.40295d-3,0.3467d-3,
     .0.26869d-3,0.189d-3,0.21036d-3/

      PROB(84)=0d0
      PROB(85)=0d0
      PROB(86)=0d0

      DO I=2,3
      MH=SMASS(I)
       DO J=1,I-1
        MHSM=SMASS(J)
        IF(MHSM.GE.MHMIN .AND. MHSM.LE.MHMAX)THEN
         SIG=BRHHH(I-1+2*(J-1))
         IF(SIG.NE.0d0)THEN

          K=1
          DO WHILE(M1(K).LE.MH .AND. K.LT.N1)
            K=K+1
          ENDDO
          IF(K.GE.2 .AND. MH.LE.M1(N1)) THEN
           L=L1(K-1)+(MH-M1(K-1))/(M1(K)-M1(K-1))*(L1(K)-L1(K-1))
           PROB(84)=PROB(84)
     .     +DDIM(CV(I)**2*VBF13(MH)*SIG*(BRBB(J)/5.824d-1)**2/L,1d0)
          ENDIF

          K=1
          DO WHILE(M2(K).LE.MH .AND. K.LT.N2)
            K=K+1
          ENDDO
          IF(K.GE.2 .AND. MH.LE.M2(N2)) THEN
           L=L2(K-1)+(MH-M2(K-1))/(M2(K)-M2(K-1))*(L2(K)-L2(K-1))
           PROB(85)=PROB(85)+DDIM(CJ(I)**2*ggF13(MH)*SIG/L,1d0)
          ENDIF

          K=1
          DO WHILE(M3(K).LE.MH .AND. K.LT.N3)
            K=K+1
          ENDDO
          IF(K.GE.2 .AND. MH.LE.M3(N3)) THEN
           L=L3(K-1)+(MH-M3(K-1))/(M3(K)-M3(K-1))*(L3(K)-L3(K-1))
           PROB(86)=PROB(86)+DDIM(CJ(I)**2*ggF13(MH)*SIG/L,1d0)
          ENDIF
         ENDIF
        ENDIF
       ENDDO
      ENDDO

      END


      SUBROUTINE CH_TB(PAR,PROB)

*     gg -> tbH+ [fb]
*     Production cross section from LHCHWGMSSMCharged CERN Twiki [pb]

      IMPLICIT NONE

      INTEGER NGGCH,NGGCHT
      INTEGER GGCHT(10000,3)

      DOUBLE PRECISION PAR(*),PROB(*),GGCH(10000,3),X,Y,Z,CHTBXS
      DOUBLE PRECISION SMASS(3),SCOMP(3,3),PMASS(2),PCOMP(2,2),CMASS
      DOUBLE PRECISION HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC
      DOUBLE PRECISION HCBRBT,HCBRWH(5),HCBRWHT,HCBRNC(5,2)
      DOUBLE PRECISION HCBRSQ(5),HCBRSL(3),HCBRSUSY,HCWIDTH

      COMMON/HIGGSPEC/SMASS,SCOMP,PMASS,PCOMP,CMASS
      COMMON/BRC/HCBRM,HCBRL,HCBRSU,HCBRBU,HCBRSC,HCBRBC,
     .      HCBRBT,HCBRWH,HCBRWHT,HCBRNC,HCBRSQ,HCBRSL,
     .      HCBRSUSY,HCWIDTH
      COMMON/GGCH/GGCH,GGCHT,NGGCH,NGGCHT
      COMMON/CHTBXS/CHTBXS

      X=CMASS
      Y=PAR(3)
      CALL TABTRI(GGCH,NGGCH,GGCHT,NGGCHT,X,Y,Z)
      CHTBXS=Z*HCBRBT*1d3

      END


      DOUBLE PRECISION FUNCTION ggF8(MH)

* SM Higgs ggF production cross section at 8 TeV in pb from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt8TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     . 65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     . 120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     . 190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     . 290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     . 390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     . 490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     . 950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     . 1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     . 1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     . 2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     . 2750d0,2800d0,2850d0,2900d0,2950d0,3000d0/

      DATA X/4.614d3,2.701d3,1.261d3,6.676d2,4.009d2,2.664d2,1.911d2,
     . 1.450d2,1.147d2,9.365d1,7.827d1,6.661d1,5.752d1,5.025d1,4.431d1,
     . 3.940d1,3.526d1,3.175d1,2.873d1,2.611d1,2.383d1,2.184d1,2.008d1,
     . 1.851d1,1.712d1,1.587d1,1.475d1,1.375d1,1.284d1,1.127d1,9.961d0,
     . 8.866d0,7.945d0,7.163d0,6.492d0,5.919d0,5.426d0,5.000d0,4.632d0,
     . 4.314d0,4.040d0,3.807d0,3.608d0,3.444d0,3.314d0,3.222d0,3.175d0,
     . 3.200d0,3.409d0,3.503d0,3.470d0,3.357d0,3.198d0,3.011d0,2.812d0,
     . 2.609d0,2.411d0,2.220d0,2.038d0,1.868d0,1.709d0,1.562d0,1.426d0,
     . 1.301d0,8.236d-1,5.241d-1,3.377d-1,2.207d-1,1.464d-1,9.849d-2,
     . 6.713d-2,4.634d-2,3.237d-2,2.285d-2,1.628d-2,1.172d-2,8.507d-3,
     . 6.222d-3,4.587d-3,3.402d-3,2.540d-3,1.909d-3,1.442d-3,1.095d-3,
     . 8.349d-4,6.400d-4,4.925d-4,3.805d-4,2.951d-4,2.297d-4,1.793d-4,
     . 1.404d-4,1.103d-4,8.685d-5,6.857d-5,5.422d-5,4.311d-5,3.419d-5,
     . 2.722d-5,2.173d-5,1.737d-5,1.390d-5,1.113d-5,8.936d-6,7.190d-6,
     . 5.771d-6,4.665d-6,3.753d-6,3.028d-6,2.455d-6,1.979d-6,1.604d-6,
     . 1.298d-6,1.048d-6/

      ggF8=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        ggF8=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      DOUBLE PRECISION FUNCTION ggF13(MH)

* SM Higgs ggF production cross section at 13 TeV in pb from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt13TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     .65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,120d0,
     .125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,190d0,200d0,
     .210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,290d0,300d0,310d0,
     .320d0,330d0,340d0,350d0,360d0,370d0,380d0,390d0,400d0,410d0,420d0,
     .430d0,440d0,450d0,460d0,470d0,480d0,490d0,500d0,550d0,600d0,650d0,
     .700d0,750d0,800d0,850d0,900d0,950d0,1000d0,1050d0,1100d0,1150d0,
     .1200d0,1250d0,1300d0,1350d0,1400d0,1450d0,1500d0,1550d0,1600d0,
     .1650d0,1700d0,1750d0,1800d0,1850d0,1900d0,1950d0,2000d0,2050d0,
     .2100d0,2150d0,2200d0,2250d0,2300d0,2350d0,2400d0,2450d0,2500d0,
     .2550d0,2600d0,2650d0,2700d0,2750d0,2800d0,2850d0,2900d0,2950d0,
     .3000d0/

      DATA X/6.996d3,4.275d3,2.085d3,1.146d3,7.103d2,4.846d2,3.555d2,
     .2.751d2,2.214d2,1.835d2,1.555d2,1.341d2,1.172d2,1.036d2,9.240d1,
     .8.302d1,7.507d1,6.825d1,6.235d1,5.720d1,5.268d1,4.869d1,4.514d1,
     .4.198d1,3.914d1,3.659d1,3.428d1,3.219d1,3.029d1,2.697d1,2.419d1,
     .2.184d1,1.984d1,1.812d1,1.665d1,1.537d1,1.426d1,1.331d1,1.248d1,
     .1.176d1,1.114d1,1.062d1,1.018d1,9.823d0,9.559d0,9.392d0,9.349d0,
     .9.521d0,1.025d1,1.063d1,1.064d1,1.040d1,1.000d1,9.516d0,8.976d0,
     .8.415d0,7.853d0,7.301d0,6.771d0,6.266d0,5.788d0,5.341d0,4.924d0,
     .4.538d0,3.008d0,2.006d0,1.352d0,9.235d-1,6.398d-1,4.491d-1,
     .3.195d-1,2.301d-1,1.675d-1,1.233d-1,9.159d-2,6.868d-2,5.193d-2,
     .3.958d-2,3.040d-2,2.349d-2,1.828d-2,1.431d-2,1.126d-2,8.913d-3,
     .7.089d-3,5.666d-3,4.546d-3,3.662d-3,2.965d-3,2.408d-3,1.963d-3,
     .1.606d-3,1.317d-3,1.084d-3,8.953d-4,7.415d-4,6.144d-4,5.108d-4,
     .4.267d-4,3.569d-4,2.988d-4,2.509d-4,2.110d-4,1.778d-4,1.499d-4,
     .1.272d-4,1.077d-4,9.121d-5,7.771d-5,6.611d-5,5.636d-5,4.809d-5,
     .4.112d-5,3.502d-5/

      ggF13=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        ggF13=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      DOUBLE PRECISION FUNCTION bbH8(MH)

* SM Higgs bbH production cross section at 8 TeV in pb from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt8TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     . 65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     . 120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     . 190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     . 290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     . 390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     . 490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     . 950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     . 1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     . 1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     . 2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     . 2750d0,2800d0,2850d0,2900d0,2950d0,3000d0/

      DATA X/6.366d1,3.458d1,2.171d1,1.450d1,9.894d0,7.077d0,5.161d0,
     . 3.900d0,2.991d0,2.351d0,1.871d0,1.510d0,1.228d0,1.012d0,8.349d-1,
     . 6.966d-1,5.842d-1,4.955d-1,4.215d-1,3.595d-1,3.094d-1,2.665d-1,
     . 2.315d-1,2.021d-1,1.763d-1,1.550d-1,1.364d-1,1.203d-1,1.066d-1,
     . 8.430d-2,6.734d-2,5.441d-2,4.434d-2,3.636d-2,3.004d-2,2.494d-2,
     . 2.088d-2,1.755d-2,1.484d-2,1.262d-2,1.076d-2,9.229d-3,7.945d-3,
     . 6.846d-3,5.935d-3,5.156d-3,4.494d-3,3.937d-3,3.451d-3,3.033d-3,
     . 2.674d-3,2.361d-3,2.095d-3,1.860d-3,1.653d-3,1.474d-3,1.317d-3,
     . 1.177d-3,1.055d-3,9.463d-4,8.517d-4,7.668d-4,6.913d-4,6.256d-4,
     . 3.832d-4,2.418d-4,1.567d-4,1.037d-4,6.987d-5,4.789d-5,3.326d-5,
     . 2.346d-5,1.672d-5,1.205d-5,8.756d-6,6.423d-6,4.747d-6,3.537d-6,
     . 2.653d-6,2.000d-6,1.517d-6,1.158d-6,8.871d-7,6.828d-7,5.273d-7,
     . 4.096d-7,3.188d-7,2.489d-7,1.949d-7,1.533d-7,1.200d-7,9.539d-8,
     . 7.518d-8,5.992d-8,4.764d-8,3.804d-8,3.028d-8,2.407d-8,1.937d-8,
     . 1.551d-8,1.245d-8,1.001d-8,8.063d-9,6.396d-9,5.217d-9,4.206d-9,
     . 3.421d-9,2.742d-9,2.215d-9,1.798d-9,1.448d-9,1.171d-9,9.427E-10,
     . 7.647E-10/

      bbH8=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        bbH8=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      DOUBLE PRECISION FUNCTION bbH13(MH)

* SM Higgs bbH production cross section at 13 TeV in pb from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt13TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     . 65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     . 120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     . 190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     . 290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     . 390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     . 490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     . 950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     . 1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     . 1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     . 2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     . 2750d0,2800d0,2850d0,2900d0,2950d0,3000d0/

      DATA X/1.138d2,6.297d1,4.011d1,2.729d1,1.931d1,1.388d1,1.029d1,
     . 7.846d0,6.102d0,4.849d0,3.914d0,3.189d0,2.633d0,2.194d0,1.838d0,
     . 1.549d0,1.317d0,1.126d0,9.669d-1,8.379d-1,7.262d-1,6.325d-1,
     . 5.534d-1,4.879d-1,4.304d-1,3.818d-1,3.383d-1,3.018d-1,2.693d-1,
     . 2.175d-1,1.769d-1,1.451d-1,1.200d-1,1.000d-1,8.397d-2,7.092d-2,
     . 6.021d-2,5.133d-2,4.410d-2,3.799d-2,3.287d-2,2.854d-2,2.491d-2,
     . 2.180d-2,1.915d-2,1.689d-2,1.491d-2,1.321d-2,1.174d-2,1.045d-2,
     . 9.324d-3,8.351d-3,7.492d-3,6.731d-3,6.046d-3,5.470d-3,4.941d-3,
     . 4.472d-3,4.057d-3,3.690d-3,3.352d-3,3.055d-3,2.784d-3,2.547d-3,
     . 1.651d-3,1.101d-3,7.518d-4,5.251d-4,3.723d-4,2.692d-4,1.967d-4,
     . 1.457d-4,1.091d-4,8.258d-5,6.299d-5,4.856d-5,3.773d-5,2.950d-5,
     . 2.322d-5,1.863d-5,1.465d-5,1.172d-5,9.329d-6,7.618d-6,6.179d-6,
     . 5.038d-6,4.120d-6,3.380d-6,2.787d-6,2.303d-6,1.910d-6,1.586d-6,
     . 1.332d-6,1.105d-6,9.255d-7,7.769d-7,6.532d-7,5.511d-7,4.742d-7,
     . 3.941d-7,3.343d-7,2.837d-7,2.415d-7,2.059d-7,1.748d-7,1.466d-7,
     . 1.265d-7,1.090d-7,9.444d-8,8.022d-8,6.968d-8,5.972d-8,5.178d-8,
     . 4.467d-8/

      bbH13=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        bbH13=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      DOUBLE PRECISION FUNCTION VBF8(MH)

* SM Higgs VBF production cross section at 8 TeV in pb from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt8TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     .65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     .120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     .190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     .290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     .390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     .490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     .950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     .1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     .1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     .2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     .2750d0,2800d0,2850d0,2900d0,2950d0,3000d0/

      DATA X/5.474d0,5.208d0,4.934d0,4.667d0,4.413d0,4.169d0,3.938d0,
     .3.722d0,3.519d0,3.329d0,3.152d0,2.986d0,2.831d0,2.686d0,2.550d0,
     .2.423d0,2.304d0,2.192d0,2.087d0,1.988d0,1.896d0,1.809d0,1.727d0,
     .1.650d0,1.577d0,1.508d0,1.443d0,1.382d0,1.323d0,1.216d0,1.120d0,
     .1.033d0,9.549d-1,8.839d-1,8.195d-1,7.609d-1,7.074d-1,6.587d-1,
     .6.140d-1,5.731d-1,5.356d-1,5.011d-1,4.692d-1,4.399d-1,4.128d-1,
     .3.877d-1,3.644d-1,3.428d-1,3.228d-1,3.042d-1,2.868d-1,2.707d-1,
     .2.556d-1,2.415d-1,2.284d-1,2.161d-1,2.046d-1,1.938d-1,1.837d-1,
     .1.742d-1,1.652d-1,1.569d-1,1.490d-1,1.415d-1,1.103d-1,8.684d-2,
     .6.892d-2,5.510d-2,4.433d-2,3.588d-2,2.918d-2,2.383d-2,1.954d-2,
     .1.608d-2,1.328d-2,1.099d-2,9.126d-3,7.593d-3,6.331d-3,5.289d-3,
     .4.425d-3,3.709d-3,3.113d-3,2.616d-3,2.200d-3,1.853d-3,1.562d-3,
     .1.318d-3,1.112d-3,9.395d-4,7.941d-4,6.714d-4,5.680d-4,4.807d-4,
     .4.068d-4,3.444d-4,2.917d-4,2.470d-4,2.092d-4,1.772d-4,1.501d-4,
     .1.271d-4,1.076d-4,9.112d-5,7.713d-5,6.527d-5,5.522d-5,4.670d-5,
     .3.948d-5,3.337d-5,2.819d-5,2.380d-5,2.009d-5,1.695d-5/

      VBF8=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        VBF8=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      DOUBLE PRECISION FUNCTION VBF13(MH)

* SM Higgs bbH production cross section at 13 TeV in pb from
* https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt13TeV

      IMPLICIT NONE

      INTEGER I,N
      PARAMETER (N=114)

      DOUBLE PRECISION MH,M(N),X(N)

      DATA M/10d0,15d0,20d0,25d0,30d0,35d0,40d0,45d0,50d0,55d0,60d0,
     .65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0,105d0,110d0,115d0,
     .120d0,125d0,130d0,135d0,140d0,145d0,150d0,160d0,170d0,180d0,
     .190d0,200d0,210d0,220d0,230d0,240d0,250d0,260d0,270d0,280d0,
     .290d0,300d0,310d0,320d0,330d0,340d0,350d0,360d0,370d0,380d0,
     .390d0,400d0,410d0,420d0,430d0,440d0,450d0,460d0,470d0,480d0,
     .490d0,500d0,550d0,600d0,650d0,700d0,750d0,800d0,850d0,900d0,
     .950d0,1000d0,1050d0,1100d0,1150d0,1200d0,1250d0,1300d0,1350d0,
     .1400d0,1450d0,1500d0,1550d0,1600d0,1650d0,1700d0,1750d0,1800d0,
     .1850d0,1900d0,1950d0,2000d0,2050d0,2100d0,2150d0,2200d0,2250d0,
     .2300d0,2350d0,2400d0,2450d0,2500d0,2550d0,2600d0,2650d0,2700d0,
     .2750d0,2800d0,2850d0,2900d0,2950d0,3000d0/

      DATA X/1.121d1,1.074d1,1.025d1,9.769d0,9.299d0,8.847d0,8.419d0,
     .8.011d0,7.627d0,7.264d0,6.924d0,6.603d0,6.301d0,6.016d0,5.748d0,
     .5.496d0,5.258d0,5.034d0,4.822d0,4.623d0,4.434d0,4.255d0,4.086d0,
     .3.925d0,3.773d0,3.629d0,3.492d0,3.362d0,3.239d0,3.010d0,2.802d0,
     .2.612d0,2.440d0,2.282d0,2.138d0,2.006d0,1.884d0,1.772d0,1.669d0,
     .1.573d0,1.485d0,1.403d0,1.326d0,1.256d0,1.190d0,1.128d0,1.071d0,
     .1.017d0,9.666d-1,9.194d-1,8.752d-1,8.337d-1,7.947d-1,7.580d-1,
     .7.235d-1,6.909d-1,6.602d-1,6.312d-1,6.038d-1,5.778d-1,5.533d-1,
     .5.301d-1,5.081d-1,4.872d-1,3.975d-1,3.274d-1,2.719d-1,2.275d-1,
     .1.915d-1,1.622d-1,1.380d-1,1.180d-1,1.013d-1,8.732d-2,7.551d-2,
     .6.550d-2,5.698d-2,4.970d-2,4.345d-2,3.807d-2,3.343d-2,2.941d-2,
     .2.592d-2,2.288d-2,2.023d-2,1.791d-2,1.588d-2,1.410d-2,1.253d-2,
     .1.115d-2,9.926d-3,8.849d-3,7.896d-3,7.052d-3,6.303d-3,5.638d-3,
     .5.046d-3,4.520d-3,4.051d-3,3.633d-3,3.259d-3,2.925d-3,2.627d-3,
     .2.360d-3,2.122d-3,1.908d-3,1.717d-3,1.545d-3,1.391d-3,1.252d-3,
     .1.128d-3,1.016d-3,9.156d-4,8.253d-4/

      VBF13=0d0

       I=1
       DOWHILE(M(I).LE.MH .AND. I.LT.N)
        I=I+1
       ENDDO
       IF(I.GT.1 .AND. MH.LE.M(N))THEN
        VBF13=X(I-1)+(MH-M(I-1))/(M(I)-M(I-1))*(X(I)-X(I-1))
       ENDIF

      END


      SUBROUTINE TABTRI(TAB,NTAB,TRI,NTRI,X,Y,Z)

      IMPLICIT NONE
      INTEGER I,NTAB,NTRI,TRI(10000,3)
      DOUBLE PRECISION TAB(10000,3),X,Y,Z,D,S1,S2,S3
      DOUBLE PRECISION X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3

      Z=0d0
      I=0
      DO WHILE(Z.EQ.0d0 .AND. I.LT.NTRI)
       I=I+1
       X1=TAB(TRI(I,1),1)
       Y1=TAB(TRI(I,1),2)
       Z1=TAB(TRI(I,1),3)
       X2=TAB(TRI(I,2),1)
       Y2=TAB(TRI(I,2),2)
       Z2=TAB(TRI(I,2),3)
       X3=TAB(TRI(I,3),1)
       Y3=TAB(TRI(I,3),2)
       Z3=TAB(TRI(I,3),3)
       S1=(X1-X)*(Y2-Y)-(Y1-Y)*(X2-X)
       IF(S1.NE.0d0)S1=DSIGN(1d0,S1)
       S2=(X2-X)*(Y3-Y)-(Y2-Y)*(X3-X)
       IF(S2.NE.0d0)S2=DSIGN(1d0,S2)
       S3=(X3-X)*(Y1-Y)-(Y3-Y)*(X1-X)
       IF(S3.NE.0d0)S3=DSIGN(1d0,S3)
       IF((S1.EQ.S2.AND.S1.EQ.S3).OR.(S1.EQ.0d0.AND.S2.EQ.S3).OR.
     .    (S2.EQ.0d0.AND.S1.EQ.S3).OR.(S3.EQ.0d0.AND.S1.EQ.S2).OR.
     .    (S1.EQ.0d0.AND.S2.EQ.S1).OR.(S2.EQ.0d0.AND.S3.EQ.S2).OR.
     .    (S3.EQ.0d0.AND.S1.EQ.S3))THEN
	D=(Y1-Y2)*(X1-X3)-(Y1-Y3)*(X1-X2)
	Z=Z1+((Z1-Z2)*((Y-Y1)*(X1-X3)+(X-X1)*(Y3-Y1))
     .     +(Z1-Z3)*((Y-Y1)*(X2-X1)+(X-X1)*(Y1-Y2)))/D
       ENDIF
      END DO

      END
