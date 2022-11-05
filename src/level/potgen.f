!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!
c***********************************************************************
c Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
c***********************************************************************
c  Version of 06/07/10 with QPAR for EMO & Rref for EMO & MLR and
c   Aubert-Frecon with retardation and  2*B(r) for MLR(Li2(A))
c***********************************************************************
      SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
     1                                                        NCN,CNN)
c** Generate analytic potential  VV(i)  as specified by the choice
c  of parameter IPOTL (see comments in PREPOT (& in main program))
c** All potentials generated in units cm-1 with absolute asymptote at
c  (input) energy VLIM for distance array  X0(i) Angstroms.
c** Return with NCN equal to power of asymptotically dominant inverse
c  power term in long range part of potential
c** Born-Oppenheimer correction functions in IPOTL=3 option may have up
c  to NBOB+1 terms.  ||    ****** last updated  06 June 2010***********
c-----------------------------------------------------------------------
      INTEGER NBOB
      PARAMETER (NBOB=20)
      INTEGER  I,J,IBOB,IAN1,IAN2,IMN1,IMN2,MN1R,MN2R,IORD,IORDD,IPOTL,
     1  PAD,QAD,PNA,NU1,NU2,NT1,NT2,NCMAX,PPAR,QPAR,NCN,NSR,NLR,NVARB,
     2  NPP,LNPT,GNS,GEL, NCMM,IVSR,LVSR,IDSTT,KDER,MM1, MMLR(9)
      CHARACTER*2 NAME1,NAME2
      REAL*8  A0,A1,A2,A3,ALFA,Asw,Rsw,BETA,BINF,B1,B2,BT,CSAV,U1INF,
     1 U2INF,T1INF,T2INF,YPAD,YQAD,YQADSM,YPNA,YPNASM,ABUND,CNN,
     2 DSCM,DX,DX1,FCT,FC1,FC2,FG1,FG2,MASS1,MASS2,RMASS1,RMASS2,REQ,
     3 Rref,Rinn,Rout,SC1,SC2,SG1,SG2,VLIM,VMIN,XDF,X1,XS,XL,XP1,ZZ,ZP,
     4 ZQ,ZME,Aad1,Aad2,Ana1,Ana2,Rad1,Rad2,Rna1,Rna2,fad1e,fad2e,FSW,
     5 ULR,ULRe,rhoAB,REQP,DM(9),DMP(9),DMPP(9),CMM(9),T0,C6adj,C9adj,
     6 RM3,RET,RETsig,RETpi,RETp,RETm,BFCT,PPOW,PVSR,
     7 U1(0:NBOB),U2(0:NBOB),T1(0:NBOB),T2(0:NBOB),PARM(50),
     8 XO(NPP),VV(NPP),RM2(NPP), bTT(-1:2),cDS(-2:0),bDS(-2:0)
      SAVE IBOB,IPOTL,IORD,IORDD,PPAR,QPAR,PAD,QAD,PNA,NSR,
     1 NLR,MMLR,NVARB,NCMM
      SAVE DSCM,REQ,Rref,PARM,U1,U2,T1,T2,CSAV,BINF,ALFA,Rsw,ZME,
     2 Aad1,Aad2,Ana1,Ana2,Rad1,Rad2,Rna1,Rna2,Rinn,Rout,ULR,ULRe,CMM
c** Damping function parameters for printout .....
      DATA bTT/2.44d0,2.78d0,3.126d0,3.471d0/
      DATA bDS/3.3d0,3.69d0,3.95d0/
      DATA cDS/0.423d0,0.40d0,0.39d0/
      SAVE bTT, bDS, cDS
c** Electron mass, as per 2006 physical constants
      DATA ZME/5.4857990943d-4/
c
      IF(LNPT.GT.0) THEN
c** Parameter definitions listed preceeding CALL in subroutine PREPOT
c-----------------------------------------------------------------------
!         READ(5,*) IPOTL, PPAR, QPAR, NSR, NLR, IBOB
!         READ(5,*) DSCM, REQ, Rref
          IF(IPOTL.GE.4) THEN
c** For MLR, DELR or Tiemann-polynomial potentials .....
!             READ(5,*) NCMM, IVSR, IDSTT, rhoAB
!             READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
              ENDIF
c-----------------------------------------------------------------------
          IF(IPOTL.EQ.1) NVARB= 0
          IF(IPOTL.EQ.2) THEN
              NVARB= NLR+2
              IORD= NLR
              IF(PPAR.EQ.0) NVARB= NLR
              ENDIF
          IF(IPOTL.EQ.3) THEN
              IORD= MAX(NSR,NLR)
              NVARB= IORD+1
              IF(PPAR.LE.0) NVARB=2
              ENDIF
          IF(IPOTL.EQ.4) THEN
              IF(NSR.LE.0) NSR= NLR
              IORD= MAX(NSR,NLR)
              NVARB= IORD+ 1
              ENDIF
          IF(IPOTL.EQ.5) THEN
              IORD= MAX(NSR,NLR)
              NVARB= IORD+ 3
              ENDIF
          IF(IPOTL.EQ.6) NVARB= 5
          IF(IPOTL.EQ.7) NVARB= NLR+ 4
c-----------------------------------------------------------------------
!         IF(NVARB.GT.0) READ(5,*) (PARM(I), I=1,NVARB)
c-----------------------------------------------------------------------
          IF(IBOB.GT.0) THEN
c-----------------------------------------------------------------------
!             READ(5,*) MN1R, MN2R, PAD, QAD, NU1, NU2, PNA, NT1, NT2
c-----------------------------------------------------------------------
              NCMAX= MAX0(NU1,NU2,NT1,NT2)
              IF(NCMAX.LT.0) THEN
                  IBOB= 0
                ELSE
c** If appropriate, read parameters & prepare to add mass-dep. BOB corrn
                  CALL MASSES(IAN1,IMN1,NAME1,GEL,GNS,MASS1,ABUND)
                  CALL MASSES(IAN1,MN1R,NAME1,GEL,GNS,RMASS1,ABUND)
                  CALL MASSES(IAN2,IMN2,NAME2,GEL,GNS,MASS2,ABUND)
                  CALL MASSES(IAN2,MN2R,NAME2,GEL,GNS,RMASS2,ABUND)
c  For simplicity, first zero out all correction function coefficients
                  DO  I=0,NCMAX
                      U1(I)= 0.d0
                      U2(I)= 0.d0
                      T1(I)= 0.d0
                      T2(I)= 0.d0
                      ENDDO
                  FC1= 0.d0
                  FC2= 0.d0
                  FG1= 0.d0
                  FG2= 0.d0
                  U1INF= 0.d0
                  U2INF= 0.d0
                  T1INF= 0.d0
                  T2INF= 0.d0
                  Aad1= 0.d0
                  Aad2= 0.d0
                  Ana1= 0.d0
                  Ana1= 0.d0
                  Rad1= 0.d0
                  Rad2= 0.d0
                  Rna1= 0.d0
                  Rna2= 0.d0
c=======================================================================
c** Read actual BOB polynomial expansion coefficients
c=======================================================================
                  IF(NU1.GE.0) THEN
c... use Huang/Le Roy form for atom-1 adiabatic potential BOB radial fx.
c-----------------------------------------------------------------------
!                     READ(5,*) U1INF,(U1(I), I=0,NU1)
c-----------------------------------------------------------------------
                      WRITE(6,630) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1        1,U1INF,PAD,PAD,PAD,PAD,PAD,PAD,NU1,QAD,QAD,QAD,QAD,QAD,
     2                                         NU1+1,(U1(I),I= 0,NU1)
                      FC1= 1.d0 - RMASS1/MASS1
                      ENDIF
c
                  IF(NU2.GE.0) THEN
c... use Huang/Le Roy form for atom-2 adiabatic potential BOB radial fx.
c-----------------------------------------------------------------------
!                     READ(5,*) U2INF,(U2(I), I=0,NU2)
c-----------------------------------------------------------------------
                      WRITE(6,630) 2,MASS2,MN2R,NAME2,IMN2,NAME2,
     1        1,U2INF,PAD,PAD,PAD,PAD,PAD,PAD,NU2,QAD,QAD,QAD,QAD,QAD,
     2                                         NU2+1,(U2(I),I= 0,NU2)
                      FC2= 1.d0 - RMASS2/MASS2
                      ENDIF
c
                  IF(NT1.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-1 ...
c-----------------------------------------------------------------------
!                     READ(5,*) T1INF,(T1(I), I=0,NT1)
c-----------------------------------------------------------------------
                      WRITE(6,634) 1,MASS1,MN1R,NAME1,IMN1,NAME1,
     1 1,T1INF,PNA,PNA,PNA,PNA,PNA,PNA,NT1,PNA,NT1+1,(T1(I),I= 0,NT1)
                      FG1= RMASS1/MASS1
                      ENDIF
c
                  IF(NT2.GE.0) THEN
c... use Huang/Le Roy centrifugal BOB radial function for atom-2 ...
c-----------------------------------------------------------------------
!                     READ(5,*) T2INF,(T2(I), I=0,NT2)
c-----------------------------------------------------------------------
                      WRITE(6,634) 2,MASS2,MN2R,NAME2,IMN2,NAME2,
     1 2,T2INF,PNA,PNA,PNA,PNA,PNA,PNA,NT2,PNA,NT2+1,(T2(I),I= 0,NT2)
                      FG2= RMASS2/MASS2
                      ENDIF
                  U1INF= U1INF*FC1
                  U2INF= U2INF*FC2
                  T1INF= T1INF*FG1
                  T2INF= T2INF*FG2
c... Now generates scaled expansion parameters all at the same time!
                  DO  I=0,NCMAX
                      U1(I)= U1(I)*FC1
                      U2(I)= U2(I)*FC2
                      T1(I)= T1(I)*FG1
                      T2(I)= T2(I)*FG2
                      ENDDO
                ENDIF
              ENDIF
          ENDIF
c
c=======================================================================
c** Generate a  Lennard-Jones(NSR,NLR)  potential here.
c=======================================================================
      IF(IPOTL.EQ.1) THEN
          XS= NSR
          XL= NLR
          XDF= DSCM/(XS-XL)
          IF(LNPT.GT.0) WRITE(6,600) NSR,NLR,DSCM,REQ
          CNN= XS*XDF*REQ**NLR
          NCN= NLR
          DO  I= 1,NPP
              VV(I)= (XL*(REQ/XO(I))**NSR - XS*(REQ/XO(I))**NLR)*XDF
     1                  +VLIM
              ENDDO
          ENDIF
c
      IF(IPOTL.EQ.2) THEN
c=======================================================================
c** Generate Seto-modified form of Surkus' GPEF function which includes
c  Dunham, SPF and OT forms as special cases.
c=======================================================================
          VMIN= VLIM
          A0= DSCM
          X1= 1.d0
          FCT= PARM(NLR+1)
          IF((PPAR.NE.0).AND.(DABS(FCT).GT.0.d0)) THEN
              FCT= 1.d0/PARM(NLR+1)
              DO  J=1,IORD
                  X1= X1+ PARM(J)*FCT**J
                  ENDDO
              DSCM= DSCM*X1*FCT**2 + VMIN
              ENDIF
          IF(PPAR.EQ.1) THEN
c  Cases with power =1 (including Dunham, SPF & O-T expansions).
              IF(DABS(PARM(NLR+1)).LE.0.d0) THEN
c ... print for Dunham expansion ...
                  WRITE(6,612) PARM(NLR+2),REQ,VMIN,A0,NLR,
     1                                              (PARM(I),I= 1,NLR)
                  NCN= -99
                  CNN= 0.d0
                  ENDIF
              IF(DABS(PARM(NLR+2)).LE.0.d0) THEN
c ... print for Simons-Parr-Finlan expansion ...
                  WRITE(6,614) PARM(NLR+1),REQ,DSCM,A0,NLR,
     1                                              (PARM(I),I= 1,NLR)
                  NCN= 1
                  ENDIF
              IF(DABS(PARM(NLR+2)-PARM(NLR+1)).LE.0.d0) THEN
c ... print for Ogilvie-Tipping expansion ...
                  WRITE(6,616) PARM(NLR+2),REQ,DSCM,A0,NLR,
     1                                              (PARM(I),I= 1,NLR)
                  NCN= 1
                  ENDIF
              ENDIF
          IF((PPAR.NE.0).AND.((PPAR.NE.1).OR.
     1                   ((DABS(PARM(NLR+2)-PARM(NLR+1)).GT.0.d0).AND.
     2                 (DABS(PARM(NLR+2)*PARM(NLR+1)).GT.0.d0)))) THEN
c ... print for general GPEF expansion variable ...
              IF(PPAR.LT.0) THEN
c ... for negative PPAR, convert to equivalent positive PPAR case
                  PPAR= -PPAR
                  A1= PARM(NLR+2)
                  PARM(NLR+2)= -PARM(NLR+1)
                  PARM(NLR+1)= -A1
                  ENDIF
              WRITE(6,618) PPAR,PPAR,PARM(NLR+1),PPAR,PARM(NLR+2),
     1                         PPAR,REQ,DSCM,A0,NLR,(PARM(I),I= 1,NLR)
              NCN= PPAR
              ENDIF
          IF(PPAR.EQ.0) THEN
c** For case of simple power series in  R  itself
              NCN= -1
              WRITE(6,620) NLR,VMIN,(PARM(I),I= 1,NLR)
              DO  I= 1, NPP
                  ZP= 1.d0
                  A1= VMIN
                  DO  J= 1,NLR
                      ZP= ZP*XO(I)
                      A1= A1+ PARM(J)*ZP
                      ENDDO
                  VV(I)= A1
                  ENDDO
              VLIM= VV(NPP)
              RETURN
              ENDIF
c ... otherwise - generate potential as a GPEF-type expansion
          DO  I= 1, NPP
              ZZ= (XO(I)**PPAR - REQ**PPAR)/(PARM(NLR+1)*XO(I)**PPAR
     1                                        + PARM(NLR+2)*REQ**PPAR)
              A1= 1.d0
              ZP= 1.d0
              DO  J=1, NLR
                  ZP= ZP*ZZ
                  A1= A1+ PARM(J)*ZP
                  ENDDO
              VV(I)= A0*ZZ*ZZ*A1 + VMIN
              ENDDO
          IF(DABS(PARM(NLR+1)).GT.0) THEN
c ...Reset asymptote to avoid spurious  E > VLIM  warnings (e.g. for HO)
              VLIM= VV(NPP)
            ELSE
              VLIM= VV(1)
            ENDIF
          ENDIF
c
c=======================================================================
c** Generate a simple Morse, or Extended (EMOp) Morse potential, or as
c  special cases, Coxon's GMO or Wei Hua's generalized Morse
c=======================================================================
      IF(IPOTL.EQ.3) THEN
          IF(Rref.LE.0.d0) Rref= REQ
          BETA= PARM(1)
          NCN= 99
          IF(LNPT.GE.0) THEN
              IF(PPAR.GT.0) THEN
c... Normal case is Morse or EMO
                  IF(IORD.EQ.0) THEN
                      WRITE(6,606) DSCM,REQ,BETA
                    ELSE
                      WRITE(6,608) PPAR,DSCM,REQ,Rref,IORD,PPAR,PPAR,
     1                            PPAR,PPAR,NVARB,(PARM(i),i= 1,NVARB)
                      IF(NSR.LT.NLR) WRITE(6,609) NSR
                      IF(NLR.LT.NSR) WRITE(6,611) NLR
                    ENDIF
                ELSE
c... Option to generate Wei Hua's extended 4-parameter Morse-type potl.
                  CSAV= PARM(2)
                  WRITE(6,605) DSCM,REQ,CSAV,BETA
                ENDIF
              ENDIF
c  Loop over distance array XO(I)
          DO  I= 1,NPP
c... for Wei Hua's extended Morse function ...
              IF(PPAR.LE.0) THEN
                  VV(I)= DSCM*((1.d0 - DEXP(-BETA*(XO(I)-REQ)))/(1.d0
     1                - CSAV*DEXP(-BETA*(XO(I)-REQ))))**2 - DSCM+ VLIM
                ELSE
c... for proper Morse or EMO function ...
                  IF(IORD.GE.1) THEN
                      ZZ= (XO(I)- Rref)/(XO(I)+ Rref)
c... for proper LeRoy-Huang yp(r) expansion ...
                      IF(PPAR.GT.1) ZZ= (XO(i)**PPAR - Rref**PPAR)/
     1                                  (XO(i)**PPAR + Rref**PPAR)
                      BETA= 0.d0
                      IF(ZZ.GT.0) IORDD= NLR
                      IF(ZZ.LE.0) IORDD= NSR
                      DO  J= IORDD,0,-1
                          BETA= BETA*ZZ+ PARM(J+1)
                          ENDDO
                      ENDIF
                  VV(I)=  DSCM*(1.d0 - DEXP(-BETA*(XO(I)-REQ)))**2
     1                                                    - DSCM+ VLIM
                ENDIF
              ENDDO
          ENDIF
c
c=======================================================================
c** Generate an MLR potential [as per J.Chem.Phys. 131, 204309 (2009)]
c=======================================================================
      IF(IPOTL.EQ.4) THEN
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= CMM(1)
              ULRe= 0.d0
              IF((NCMM.GE.3).AND.(MMLR(2).LE.0)) THEN
c** For special Aubert-Frecon Li2(^2P + ^2S) {3,0,6,6} type cases
                  RM3= 1.d0/REQ**3
c** Include QED retardation function in C3 terms for Li2(A) !!
c  NOTE ... the numerical factor here is  2\pi/\lambda  for this case
                  RET= 9.36423830d-4*REQ
                  RETSig= DCOS(RET) + (RET)*DSIN(RET)
                  RETPi= RETSig - RET**2 *DCOS(RET)
                  RETp= RETSig + 0.5d0*RETPi
                  RETm= RETSig - 0.5d0*RETPi
                  C6adj= CMM(3) + CMM(1)**2/(4.d0*DSCM)
                  C9adj = CMM(1)*C6adj/(2.d0*DSCM)
                  IF(MMLR(2).EQ.0) THEN
                      T0= CMM(1)*RETm*RM3
                      ULRe= 0.5d0*(-CMM(2) + CMM(1)*RETp*RM3)
                      IF(NCMM.GE.3) THEN
                          T0= T0 + C6adj*RM3**2
                          ULRe= ULRe + 0.5d0*C6adj*RM3**2
                          IF(NCMM.GE.4) THEN
                              T0= T0+ CMM(4)*(RM3/REQ)**2
                              ULRe= ULRe + 0.5d0*CMM(4)*(RM3/REQ)**2
     1                                                  + C9adj*RM3**3
                              ENDIF
                          ENDIF
                      T0= T0/3.d0
                      ULRe= ULRe+0.5d0*DSQRT((T0-CMM(2))**2+ 8.d0*T0**2)
                      ENDIF
                  IF(MMLR(2).EQ.-1) THEN
                      CALL AF3X3LEV(REQ,CMM(2),CMM(1),C6adj,
     1                                               CMM(4),DSCM,ULRe)
                      ULRe= ULRe + C9adj*RM3**3
                      ENDIF
                ELSE
                  IF(rhoAB.GT.0.d0) THEN
                      KDER= 0
                      CALL dampF(REQ,rhoAB,NCMM,MMLR,IVSR,IDSTT,KDER,
     1                                                    DM,DMP,DMPP)
                      ENDIF
                  DO  J= 1,NCMM
                      IF(rhoAB.LE.0.d0) THEN
                          ULRe= ULRe + CMM(J)/REQ**MMLR(J)
                        ELSE
                          ULRe= ULRe + DM(J)*CMM(J)/REQ**MMLR(J)
                        ENDIF
                      ENDDO
                ENDIF
              BINF= DLOG(2.d0*DSCM/ULRe)
              WRITE(6,602) NCN,PPAR,QPAR,DSCM,REQ
c... use THEOCHEM/Huang form:  \beta(yp)= Binf*yp + [1-yp]*{power series in yq}
              WRITE(6,607) PPAR,PPAR,QPAR,NSR,NLR,IORD+1,
     1                                       (PARM(J),J= 1,IORD+1)
              IF(Rref.GT.0) THEN
                  WRITE(6,613) Rref
                ELSE
                  WRITE(6,615) REQ
                  Rref= REQ
                ENDIF
              IF(rhoAB.GT.0.d0) THEN
                  PVSR= 0.5d0*IVSR
                  IF(IDSTT.GT.0) THEN
                      PVSR= 0.5d0*IVSR
                      WRITE(6,664) rhoAB,PVSR,bDS(IVSR),cDS(IVSR),PVSR
                    ELSE
                      LVSR= IVSR/2
                      WRITE(6,666) rhoAB,LVSR,bTT(LVSR)
                    ENDIF
                ELSE
                  WRITE(6,668)
                ENDIF
              WRITE(6,617) BINF,MMLR(1),CMM(1),MMLR(1)
              IF(NCMM.GT.1) THEN
                  MM1= 2
                  IF(MMLR(2).EQ.0) THEN
                      MM1= 3
                      WRITE(6,623) MMLR(2),CMM(2),MMLR(2)
                      ENDIF
                  DO  I= MM1,NCMM
                      IF(MMLR(I).LE.9) WRITE(6,619) MMLR(I),CMM(I)
     1                                                    ,MMLR(I)
                      IF(MMLR(I).GT.9) WRITE(6,621) MMLR(I),CMM(I)
     1                                                    ,MMLR(I)
                      ENDDO
                  ENDIF
              ENDIF
c  Loop over distance array XO(I)
          DO  I= 1,NPP
              ZZ= (XO(i)**PPAR- REQ**PPAR)/(XO(i)**PPAR+ REQ**PPAR)
              ZP= (XO(i)**PPAR-Rref**PPAR)/(XO(i)**PPAR+Rref**PPAR)
              ZQ= (XO(i)**QPAR-Rref**QPAR)/(XO(i)**QPAR+Rref**QPAR)
              BETA= 0.d0
              IF(ZZ.GT.0) IORDD= NLR
              IF(ZZ.LE.0) IORDD= NSR
              DO  J= IORDD,0,-1
                  BETA= BETA*ZQ+ PARM(J+1)
                  ENDDO
c  Calculate MLR exponent coefficient
              BETA= BINF*ZP + (1.d0- ZP)*BETA
              ULR= 0.d0
c** Calculate local value of uLR(r)
              IF((NCMM.GE.3).AND.(MMLR(2).LE.0)) THEN
c... For special Aubert-Frecon Li2(^2P + ^2S) {3,0,6,6} type case
                  RM3= 1.d0/XO(I)**3
c... Include QED retardation function in C3 terms for Li2(A) !!
c  NOTE ... the numerical factor here is  2\pi/\lambda  for Li2(A)
                  RET= 9.36423830d-4*XO(I)
                  RETSig= DCOS(RET) + (RET)*DSIN(RET)
                  RETPi= RETSig - RET**2 *DCOS(RET)
                  RETp= RETSig + 0.5d0*RETPi
                  RETm= RETSig - 0.5d0*RETPi
                  IF(MMLR(2).EQ.0) THEN
c... Aubert-Frecon 2x2 case - for A-state Li2
                      T0= CMM(1)*RETm*RM3
                      ULR= 0.5d0*(-CMM(2) + CMM(1)*RETp*RM3)
                      IF(NCMM.GE.3) THEN
                          T0= T0 + C6adj*RM3**2
                          ULR= ULR + 0.5d0*C6adj*RM3**2
                          IF(NCMM.GE.4) THEN
                              T0= T0+ CMM(4)*(RM3/XO(I))**2
                              ULR= ULR + 0.5d0*CMM(4)*(RM3/XO(I))**2
     1                                                  + C9adj*RM3**3
                              ENDIF
                          ENDIF
                      T0= T0/3.d0
                      ULR= ULR+ 0.5d0*DSQRT((T0- CMM(2))**2+ 8.d0*T0**2)
                      ENDIF
                  IF(MMLR(2).EQ.-1) THEN
c... for Aubert-Frecon 3x3 case yielding lowest (c state) energy
                      CALL AF3X3LEV(XO(I),CMM(2),CMM(1),C6adj,
     1                                                CMM(4),DSCM,ULR)
                      ULR= ULR + C9adj*RM3**3
                      ENDIF
                ELSE
c** For the 'regular' simple inverse-power sum case.
                  IF(rhoAB.GT.0.d0) CALL dampF(XO(I),rhoAB,NCMM,MMLR,
     1                                    IVSR,IDSTT,KDER,DM,DMP,DMPP)
                  DO  J= 1,NCMM
                      IF(rhoAB.LE.0.d0) THEN
                          ULR= ULR + CMM(J)/XO(I)**MMLR(J)
                        ELSE
                          ULR= ULR + DM(J)*CMM(J)/XO(I)**MMLR(J)
                        ENDIF
                      ENDDO
                ENDIF
              BETA= (ULR/ULRe)*DEXP(-BETA*ZZ)
              VV(I)= DSCM*(1.d0 - BETA)**2 - DSCM + VLIM
              ENDDO
          ENDIF
c
c=======================================================================
c** Generate a DELR potential [as per JCP 119, 7398 (2003)]
c=======================================================================
      IF(IPOTL.EQ.5) THEN
          IF(LNPT.GT.0) THEN
              rhoAB= PARM(NVARB-1)
ccc           IVSR= 1                               % Read IVSR explicitly
ccc           IF(PARM(NVARB).LE.0.d0) IVSR= 2       % Read IVSR explicitly
              PPOW= PPAR
              REQP= REQ**PPOW
              A1= 0.0d0
              B1= 0.0d0
c... first, get  AA & BB and their derivatives!
              KDER= 1
              CALL dampF(REQ,rhoAB,NCMM,MMLR,IVSR,IDSTT,KDER,DM,DMP,
     1                                                           DMPP)
              DO  J= 1,NCMM
                  A1= A1+ CMM(J)*DM(J)/REQ**MMLR(J)*(1.d0+ DMP(J)/
     1                        (PARM(1)*DM(J)) - MMLR(J)/(PARM(1)*REQ))
                  B1= B1+ CMM(J)*DM(J)/REQ**MMLR(J)*(2.d0+ DMP(J)/
     1                        (PARM(1)*DM(J)) - MMLR(J)/(PARM(1)*REQ))
                  ENDDO
              A1 = A1 + DSCM
              B1 = B1 + 2.0d0*DSCM
              WRITE(6,650) PPAR,DSCM,REQ,NSR,NLR,(PARM(I),I= 1,IORD+1)
              WRITE(6,652) PPAR,PPAR,PPAR,PPAR,PPAR
              IF(Rref.GT.0.d0) WRITE(6,654) Rref
              IF(Rref.LE.0.d0) WRITE(6,656) REQ
              WRITE(6,658) A1,B1,NCMM,(MMLR(J),CMM(J),J= 1,NCMM)
              IF(IDSTT.GT.0) THEN
                  PVSR= 0.5d0*IVSR
                  WRITE(6,664) rhoAB,PVSR,bDS(IVSR),cDS(IVSR),PVSR
                ELSE
                  LVSR= IVSR/2
                  WRITE(6,666) rhoAB,LVSR,bTT(LVSR)
                ENDIF
              ENDIF
c** Now ... generate potential function array for DELR form
          DO  I= 1, NPP
              ZZ= (XO(I)**PPOW -REQP)/(XO(I)**PPOW +REQP)
              BETA= 0.d0
              IF(ZZ.GT.0) IORDD= NLR
              IF(ZZ.LE.0) IORDD= NSR
c ... calculate the exponent
              DO  J= IORDD,0,-1
                  BETA= BETA*ZZ+ PARM(J+1)
                  ENDDO
              BETA= DEXP(-BETA*(XO(I)-REQ))
c ... calculate the (damped) long-range tail
              A3= 0.0d0
              KDER= 0
              CALL dampF(XO(I),rhoAB,NCMM,MMLR,IVSR,IDSTT,KDER,DM,DMP,
     1                                                           DMPP)
              DO  J= 1, NCMM
                  A3= A3+ DM(J)*CMM(J)/XO(I)**MMLR(J)
                  ENDDO
              VV(I)=  (A1*BETA - B1)*BETA + A3 + VLIM
              ENDDO
          ENDIF
c
      IF(IPOTL.EQ.6) THEN
c=======================================================================
c** For generalized  HFD(m= MMLR(j), j=1,NCMM) potential with reduced
c  form   VBAR = ALFA*x**PARM(5) * exp[-BETR*x - PARM(4)*x**2] - D(x)*
c      [CMM(1)/x**MMLR(1) + CMM(2)/x**sMMLR(2) + CMM(3)/x**MMLR(3) + ...
c  x=r/R_e ,  VBAR= V/De   and    D(x) = 1 for x > PARM(2)   and
c      D(x)= exp[-PARM(1)*(PARM(2)/x - 1)**PARM(3)] for  x < PARM(2)
c=======================================================================
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= CMM(1)*DSCM*REQ**MMLR(1)
              A1= PARM(1)
              A2= PARM(2)
              A3= PARM(3)
              B2= PARM(4)
              DX= 1.d0
              DX1= 0.d0
              IF(A2.GT.1.d0) THEN
                  DX= DEXP(-A1*(A2- 1.d0)**A3)
                  DX1= A1*A2*A3*DX*(A2- 1.d0)**(A3- 1.d0)
                  ENDIF
              ALFA= 0.d0
              DO  J= 1, NCMM
                  ALFA= ALFA+ CMM(J)
                  ENDDO
              ALFA= ALFA*DX -1.D0
              IF(ALFA.LE.0.d0) THEN
                  WRITE(6,622) ALFA,(MMLR(J),CMM(J),J= 1, NCMM)
c                 STOP
                  ENDIF
              B1= 0.d0
              DO  J= 1, NCMM
                  B1= B1+ (MMLR(J)*DX - DX1)*CMM(J)
                  ENDDO
              B1= B1/ALFA + PARM(5) - 2.d0*B2
              ALFA= ALFA*DEXP(B1+B2)
              WRITE(6,624) PARM(5),B1,B2,ALFA*DSCM,
     1                                     (MMLR(J),CMM(J),J= 1, NCMM)
              WRITE(6,626) DSCM,REQ,A1,A2,A3
              ENDIF
          DO  I= 1,NPP
              X1= XO(I)/REQ
              XP1= 0.0D0
              IF((B1*X1+ B2*X1**2).LT.170.D0) XP1= DEXP(-X1*(B1+ B2*X1))
              XP1= XP1*X1**PARM(5)
              FC1= 0.d0
              DO  J= 1, NCMM
                  FC1= FC1 + CMM(J)/X1**MMLR(J)
                  ENDDO
              IF(X1.LT.A2) FC1= FC1*DEXP(-A1*(A2/X1- 1.d0)**A3)
              VV(I)= DSCM*(ALFA*XP1- FC1) + VLIM
              ENDDO
          ENDIF
c
      IF(IPOTL.EQ.7) THEN
c=======================================================================
c** Generate Tiemann-type polynomial potential attached to inverse-power
c  tail and 1/R^{12} (or exponential) inner wall [PRA 63, 012710 (2000)].
c  Polynomial expansion variable is  z= [R - Rm]/[R + b*Rm] where
c  expansion has constant and linear terms.  The read-in DSCM= De (well
c  depth), but  Rm (read in as REQ) is not precisely Re (for a1 .neq. 0).
c  NCMM= number of inverse-power long-range terms;
c  NVARB= (polynomial order) + 4.  [PPAR and NSR are dummy parameters]
c** Read-in parameters PARM(i) are in order: the  (NLR+1)  polynomial
c  coefficients  a(0) - a(NLR), the expansion variable denominator
c  factor b=PARM(NLR+2), and the the inner and outer bounds on the
c  polynomial domain, Tiemann's Rinn= PARM(NLR+3) & Rout= PARM(NLR+4),
c  respectively.  The powers and coefficients (-ve if attractive) of the
c  NCMM inverse-power long-range terms are MMCM(j) and CMM(j).
c=======================================================================
          IF(LNPT.GT.0) THEN
              NCN= MMLR(1)
              CNN= -CMM(1)
              A0= VLIM- DSCM
              BT= PARM(NLR+2)
              Rinn= PARM(NLR+3)
              Rout= PARM(NLR+4)
c** Determine analytic function attaching smoothly to inner wall of
c  polynomial expansion at  R= Rinn < Rm
              ZZ= (Rinn - REQ)/(Rinn+ BT*REQ)
              ZP= 1.d0
              A1= PARM(1)
              A2= 0.d0
              DO  J= 1,NLR
                  A2= A2+ J*ZP*PARM(J+1)
                  ZP= ZP*ZZ
                  A1= A1+ ZP*PARM(J+1)
                  ENDDO
              A2= A2*(REQ+ BT*REQ)/(Rinn + BT*REQ)**2
c* If inward extrapolation is exponential:   A1*exp(-A2*(R-Rinn))
c             A2= -A2/A1
c* If inward extrapolation is inverse-power:   A1 + A2/R**12
              A2= -A2*Rinn**13/12.d0
              A1= A1 - A2/Rinn**12 + VLIM - DSCM
c** With long-range tail an NCMM-term inverse-power sum, add 1 additional
c   higher-power term to ensure continuity (not smoothness) at  Rout
c** NOTE attractive long-range terms have negative (-) coefficients!
              ZZ= (Rout - REQ)/(Rout+ BT*REQ)
              ZP= 1.d0
              B1= PARM(1)
              DO  J= 1,NLR
                  ZP= ZP*ZZ
                  B1= B1+ ZP*PARM(J+1)
                  ENDDO
              A3= DSCM
              DO  J= 1,NCMM
                  A3= A3+ CMM(J)/Rout**MMLR(J)
                  ENDDO
              PPAR= NCMM+ 1
              MMLR(PPAR)= MMLR(NCMM)+ 2
              CMM(PPAR)= (B1-A3)*Rout**MMLR(PPAR)
c*** Print for Tiemann-type potential
              IF(LNPT.GE.0) THEN
                  WRITE(6,640) DSCM,REQ,PARM(NLR+2),NLR,NLR+1,
     1                                            (PARM(J),J= 1,NLR+1)
ccc               IF(XO(1).LT.Rinn) WRITE(6,642) PARM(NLR+3),A1,A2,A0
                  IF(XO(1).LT.Rinn) WRITE(6,642) PARM(NLR+3),A1,A2
                  IF(XO(NPP).GT.Rout) WRITE(6,644) PARM(NLR+4),
     1                                     (CMM(J),MMLR(J),J= 1, PPAR)
                  ENDIF
              ENDIF
c ... now generate potential as a Tiemann-type expansion
          DO  I= 1, NPP
              IF(XO(I).LE.Rinn) THEN
c ... for exponential inward extrapolation ...
c                 VV(I)= A1*DEXP(-A2*(XO(I)- Rinn)) + A0
c ... for   A + B/R**12  inward extrapolation ...
                  VV(I)= A1 + A2/XO(I)**12
                ELSEIF(XO(I).LE.Rout) THEN
                  ZZ= (XO(I) - REQ)/(XO(I) + BT*REQ)
                  A3= A0 + PARM(1)
                  ZP= 1.d0
                  DO  J= 1,NLR
                      ZP= ZP*ZZ
                      A3= A3+ PARM(J+1)*ZP
                      ENDDO
                  VV(I)= A3
                ELSEIF(XO(I).GT.Rout) THEN
                  A3= VLIM
                  DO  J= 1, PPAR
                      A3= A3+ CMM(J)/XO(I)**MMLR(J)
                      ENDDO
                  VV(I)= A3
                ENDIF
              ENDDO
          ENDIF
c
      IF(IBOB.GT.0) THEN
c=======================================================================
c** If appropriate, generate Born-Oppenheimer breakdown correction
c      functions to rotationless and/or centrifugal potential(s) using
c      LeRoy/Huang radial functions ...
c=======================================================================
          DO  I=1,NPP
              YPAD= (XO(I)**PAD- REQ**PAD)/(XO(I)**PAD+ REQ**PAD)
              YQAD= (XO(I)**QAD- REQ**QAD)/(XO(I)**QAD+ REQ**QAD)
              YPNA= (XO(I)**PNA- REQ**PNA)/(XO(I)**PNA+ REQ**PNA)
              SC1= U1INF*YPAD
              SC2= U2INF*YPAD
              SG1= T1INF*YPNA
              SG2= T2INF*YPNA
              YQADSM= (1.d0- YPAD)
              YPNASM= (1.d0- YPNA)
c ... finally, accumulate overall BOB terms ... all at the same time!
              DO  J= 0,NCMAX
                  SC1= SC1+ YQADSM*U1(J)
                  SC2= SC2+ YQADSM*U2(J)
                  SG1= SG1+ YPNASM*T1(J)
                  SG2= SG2+ YPNASM*T2(J)
                  YQADSM= YQADSM*YQAD
                  YPNASM= YPNASM*YPNA
                  ENDDO
              RM2(I)= (1.d0+ SG1+ SG2)/XO(i)**2
              VV(I)= VV(I) + SC1 + SC2
              ENDDO
          IF((IAN1.EQ.3).AND.(IAN2.EQ.3).AND.(IMN1.NE.IMN2).AND.
     1                                            (MMLR(2).EQ.0)) THEN
c!! For mixed isotopopogue {6,7}Li_2(A) state, shift asymptote!
              DO  I= 1,NPP
                  RM3= (2.d0/3.d0)*CMM(1)/XO(I)**3
                  VV(I)= VV(I) + RM3 - DSQRT(RM3**2 + 3.085959756d-02)
                  ENDDO
              ENDIF
          ENDIF
      RETURN
  600 FORMAT(/' Lennard-Jones(',I2,',',I2,') potential with   De=',
     1  F10.3,'(cm-1)   Re =',F10.6,'(A)')
  601 FORMAT(' *** Input error:  NSR=',i3,'  NL=',i3,'   and  NVARB=',
     1 i3,'  inconsistent so  STOP !!!!')
  602 FORMAT(/' MLR(n=',i1,'; p=',I1,', q=',I1,') Potential with:   De='
     1 ,F10.3,'[cm-1]    Re=',F12.8,'[A]')
  605 FORMAT(/' Potential is a Hua-Wei 4-parameter Morse type function w
     1ith   De =',F11.4/11x,'Re =',F12.9,'   C=',f7.4,'   &   beta=',
     1  F13.10,' [1/Angstroms]')
  606 FORMAT(/' Potential is a simple Morse function with   De =',F11.4,
     1  '    Re =',F12.9/39x,'and   beta =',F13.10,' [1/Angstroms]')
  607 FORMAT('   with exponent coefficient   beta(r)= beta{INF}*y',I1,
     1  ' + [1-y',i1,']*Sum{beta_i*y',i1,'^i}'/6x,'exponent coefft. powe
     2r series orders',I4,' for  R < Re  and',I4,' for  R > Re'/6x,
     3  'and',i3,' coefficients:',1PD16.8,2D16.8:/(10x,4D16.8:))
  608 FORMAT(/' EMO_',i1,' Potential with   De=',F11.4,'    Re=',F11.8,
     1 '   Rref=',F11.8/3x,'Exponent coeft: order-',i2,
     2 ' power series in  y=(r**',i1,' - Rref**',i1,')/(r**',i1,
     3 ' + Rref**',i1,')'/'   with',I3,' coefficients:',1x,1PD18.9,
     4 2D18.9:/(7X,4D18.9:))
  609 FORMAT('   where for  r < Re  polynomial order is truncated to ord
     1er   NSR=',i2)
  610 FORMAT(/' Potential is Generalized Morse Oscillator with   De=',
     1 F10.3,'   Re=',F11.8/4x,'Exponent factor "beta" is',i3,' order po
     2wer series in (r-Re) with coefficients:'/4x,1PD18.9,3D18.9:/
     3 (4X,4D18.9:))
  611 FORMAT('   where for  r > Re  polynomial order is truncated to ord
     1er   NLR=',i2)
  612 FORMAT(/' Potential is a Dunham expansion in  (r-Re)/(',f5.2,
     1  ' * Re)  with   Re=',f12.9/'  V(Re)=',f12.4,'    a0=',1PD16.9,
     2  '   and',i3,'  a_i coefficients:'/(5D16.8))
  613 FORMAT(6x,'with radial variables  y_p & y_q  defined w.r.t.',
     1 '  Rref=',F10.7)
  615 FORMAT(6x,'radial variables  y_p & y_q  defined w.r.t.',
     1 '  Rref= Re=' F10.7)
  617 FORMAT('   while  betaINF=',f12.8,'  & uLR defined by  C',i1,' =',
     1  1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  619 FORMAT(50x,'C',I1,' =',1PD13.6,'[cm-1 Ang','^',0P,I1,']')
  621 FORMAT(50x,'C',I2,'=',1PD13.6,'[cm-1 Ang','^',0P,I2,']')
  623 FORMAT(3x,'Use Aubert-Frecon model for uLR(r) with',
     1  8x,'C',I1,' =',1PD13.6,'[cm-1 Ang^',i1,']'/8x,'including retarda
     2tion & B(r) in u^{tot} plus Delta(V)_{gu}^{(6,7)})')
  614 FORMAT(/' Potential is an SPF expansion in  (r-Re)/(',F5.2,
     1  '* r)  with   Re=',f12.9/5x,'De=',g18.10,'   b0=',
     2  1PD16.9,'   and',i3,'  b_i  coefficients:'/(5D16.8))
  616 FORMAT(/' Potential is an O-T expansion in  (r-Re)/[',f5.2,
     1  '*(r+Re)]  with   Re=',f12.9/5x,'De=',G18.10,
     2  '   c0=',1PD16.9,'   and',i3,'  c_i coefficients:'/(5D16.8))
  618 FORMAT(/' Potential is a general GPEF expansion in  (r**',i1,
     1  ' - Re**',i1,')/(',SP,F5.2,'*r**',SS,i1,SP,F6.2,'*Re**',SS,i1,
     2  ')'/5x,'with   Re=',f12.9,'   De=',g18.10,'   g0=',1PD16.9/
     3  5x,'and',i3,'  g_i coefficients:  ',3D16.8/(5D16.8:))
  620 FORMAT(/' Potential is a power series in  r  of  order',i3,
     1 ' with   V(r=0)=',f11.4/3x,'& coefficients (from linear term):',
     2 1P2d16.8:/(5x,4D16.8:))
  622 FORMAT(/' *** ERROR in generating HFD potential *** generate   ALF
     1A=',G15.7,'  from reduced  Cm  coefficients:'/
     2   (3x,3('   C',i2,' =',1PD15.7:)) )
  624 FORMAT(/' Potential is Generalized HFD with exponent factors   gam
     1ma=',f9.6/'   beta1=',f12.8,'   beta2=',f9.6,'   A=',1PD16.9,
     2 "   & reduced Cm's:"/(3x,3('   C',i2,' =',D15.7:)) )
  626 FORMAT('   De=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]   and'/
     1  '     Damping function  D(r)= exp[ -',0Pf6.4,'*(',f7.4,
     1  '/X -1.0)**',f5.2,']')
  630 FORMAT(/' BOB adiabatic potential correction for atom-',I1,
     1 '  of mass ',f15.11/'   consists of mass factor  [1- MASS(',I3,
     2 A2,')/MASS(',I3,A2,')]  multiplying all of:'/5x,'u',i1,'INF=',
     3 f11.6,'  times  y',i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1,
     4 ' + Re**',i1,')]'/5x,'plus  [1 - y',i1,']  times an order',I3,
     5 ' polynomial in'/7x,'y',i1,'=[(r**',i1,' - Re**',i1,')/(r**',i1,
     6 ' + Re**',i1,')]  with the ',i3,' coefficients:'/(3x,4G17.9:))
  632 FORMAT(/' BOB adiabatic potential correction for atom-',I1,
     1  '  of mass ',f15.11/'   consists of mass factor  m{electron}*[1/
     2MASS(',I3,A2,') - 1/MASS(',I3,A2,')]'/5x,'multiplying   u',i1,
     3 'INF=',1PD17.9,'  times [1 - fsw(r)/fsw(Re)]'/'   plus  fsw(r)  t
     4imes an order',0P,i3,' polynomial in z{O-T} with coefficients:' /
     5  (3x,1P4D17.9:))
  634 FORMAT(/' BOB centrifugal correction for atom-',I1,'  of mass ',
     1 f15.11/3x,'consists of mass factor  [MASS(',I3,A2,')/MASS(',I3,
     2 A2,')]  multiplying all of:'/5x,'q',i1,'INF=',F11.6,' times  y',
     3 i1,'= [(r**',i1,' - Re**',i1,')/(r**',i1,' + Re**',i1,')]'/
     4 5x,'plus [1 - y',i1,'] times an order',I3,' polynomial in y',i1,
     5 ' with the',i3,' coefficients:'/(3x,4G17.9:))
  636 FORMAT(3x,'where   fsw(r) = 1/[1 - exp{',f7.4,'*(r -',f7.4,')}]')
  638 FORMAT(/' BOB centrifugal correction for atom-',I1,'  of mass ',
     1 f15.11/3x,'consists of mass factor   [mass{electron}/MASS(',I3,
     2 A2,')]'/'   multiplying   q',i1,'INF=',1PD17.9,'  times [1 - fsw(
     3r)/fsw(Re)]'/ '   plus  fsw(r)  times an order',0P,i3,' polynomial
     4 in z{O-T} with coefficients:'/ (3x,4G17.9:))
  640 FORMAT(/' Tiemann-type potential with   De=',F11.4,'   Rm=',f9.6,
     1 '   is a power series'/10x,'in  (r - Re)/(r ',SP,F9.5,
     2 '*Re) of order',SS,I3,'  with the',I3,' coefficients:'/(5D16.8))
c 642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',1PD13.6,'*exp[-',
c    1  0P,F9.6,'*(r - Rinn)] ',SP,F10.3)
  642 FORMAT(' where for  r < Rinn=',F7.4,'   V=',SP,F12.4,1x,1PD13.6,
     1  '/R**12' )
  644 FORMAT('  and  for  r > Rout=',F7.3,'   V= VLIM ',
     1 (SP,1PD14.6,'/r**',SS,I2):/(39x,SP,1PD14.6,'/r**',SS,I2))
  650 FORMAT(/' DELR(p=',i2,') Potential with   De=', F11.4,'[cm-1]   Re
     1=',F11.8,'[A]   where'/3x,'exponent coefft. has power series order
     2s',I4,' for  R < Re  and',I3,' for  R > Re'/6x,'with polynomial co
     3efficients',8x,1PD17.8,D17.8/(8x,4D17.8))
  652 FORMAT(6x,'where the radial variable   y_',I1,'= (r**',I1,' - Rref
     4**',i1,')/(r**',I1,' + Rref**',i1, ')')
  654 FORMAT(10x,'is defined w.r.t.   Rref=',F11.8)
  656 FORMAT(10x,'is defined w.r.t.   Rref= Re= ',F11.8)
  658 FORMAT(3x,'Generate   A(DELR)=',1Pd17.9,'   B(DELR)=',D17.9/
     1 6x,'from uLR defined by',I2," inverse-power terms with coeffts (+
     2've repulsive):"/(5x,3(5x,'C',0P,i2,' =',1Pd14.6:)))
  664 FORMAT(4x,'uLR inverse-power terms incorporate DS-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',SS,f4.1,'}'/8x,'Dm(r)= [1 - exp(-',f5.2,
     3 '(rhoAB*r)/m -',f6.3,'(rhoAB*r)^2/sqrt{m})]^{m',SP,F4.1,'}')
  666 FORMAT(4x,'uLR inverse-power terms incorporate TT-type damping wit
     1h   rhoAB=',f10.7/8x,'defined to give very short-range  Dm(r)*Cm/r
     2^m  behaviour   r^{',I2,'}'/8x,'Dm(r)= [1 - exp(-bTT*r)*SUM{(bTT*r
     3)^k/k!}]   where   bTT=',f6.3,'*rhoAB')
  668 FORMAT(4x,'uLR inverse-power terms incorporate NO damping function
     1s')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE dampF(r,rhoAB,NCMM,MMLR,IDF,IDSTT,KDER,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen version of the incomplete
c    gamma function damping function, for  m= 1 to MMAX.
c---------------------- RJL Version of 06 July 2010 --------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* r - the radial distance in Angsroms (!)
c* RHOab  'universal' scaling coefficient used for systems other than H_2
c       RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c              where I_p^A is the ionization potential of atom A
c              and I_p^H is the ionization potential of atomic hydrogen
c* NCMM  the number of inverse-power terms to be considered
c* MMLR  are the powers of the NCMM inverse-power terms
c* IDF requires damping to be defined s.th.  Dm(r)/r^m --> r^{IDF/2}
c* IDSTT specifies damping function type:  > 0  use Douketis et al. form
c                               if  IDSTT .LE. 0  use Tang-Toennies form
c* KDER:  if KDER.GT.0  the first derivative is also calculated
c*        if KDER.GT.1  the second derivative is also calculated
c-----------------------------------------------------------------------
c                 Upon Output
c  DM(m) - The value of the damping function for the long range term
c          C_MMLR(m)/r^MMLR(m)    {m= 1, NCMM}
c  DMP(m) - The first derivative of the damping function  DM(m)
c  DMPP(m) - The second derivative of the damping function  DM(m)
c-----------------------------------------------------------------------
      INTEGER NCMM,NCMMax,MMLR(NCMM),IDF,IDSTT,KDER,IDFF,FIRST,
     1  Lsr,m,MM,MMAX
      REAL*8 r,rhoAB,bTT(-2:2),cDS(-4:0),bDS(-4:0),aTT,br,XP,YP,
     1  TK, DM(NCMM),DMP(NCMM),DMPP(NCMM),SM(-3:25),
     2  bpm(20,-2:0), cpm(20,-2:0),ZK
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J.Chem.Phys. 50, 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA bTT/2.10d0,2.44d0,2.78d0,3.126d0,3.471d0/
       DATA bDS/2.50d0,2.90d0,3.3d0,3.69d0,3.95d0/
       DATA cDS/0.468d0,0.446d0,0.423d0,0.40d0,0.39d0/
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpm
c------------------------------------------------------------------------
      IF(RHOab.LE.0) THEN
          WRITE(6,602) RHOab
c         STOP
          ENDIF
      IF(IDSTT.LE.0) THEN
c===========================================
c** For Tang-Toennies type damping functions
c===========================================
          IF((IDF.LT.-4).OR.(IDF.GT.4)) THEN
                WRITE(6,600) IDSTT,IDF
c               STOP
                ENDIF
          Lsr= IDF/2
          MMAX= MMLR(NCMM) + Lsr - 1
          aTT= RHOab*bTT(Lsr)
          br= aTT*r
          XP= DEXP(-br)
          SM(-3)= 0.d0
          SM(-2)= 0.d0
          SM(-1)= 0.d0
          SM(0)=  1.d0
          TK= 1.d0
          IF(br.GT.0.5d0) THEN
              DO  m= 1,MMAX
                  TK= TK*br/DFLOAT(m)
                  SM(m)= SM(m-1)+ TK
                  ENDDO
              DO m= 1, NCMM
                  MM= MMLR(m) - 1 + Lsr
                  DM(m)= 1.d0 - XP*SM(MM)
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*(SM(MM) - SM(MM-1))
                      IF(KDER.GT.1) DMPP(m)= -aTT*aTT*XP*(SM(MM)
     1                                     - 2.d0*SM(MM-1) + SM(MM-2))
                      ENDIF
                  ENDDO
c-----------------------------------------------------------------------
c  The above section handles the calculation of the value of the damping
c  function for most values of r.  However, at very small r that algorithm
c  becomes unstable due to numerical noise.  To avoid this, if the
c  argument is very small it is re-evaluated as a finite sum ...
c-----------------------------------------------------------------------
            ELSE
              MMAX= MMAX+5
              DO  m= 1, MMAX
c... NOTE that here SM(m) is the m'th term  (b*r)^m/m!  [not a sum]
                  SM(m)= SM(m-1)*br/DFLOAT(m)
                  ENDDO
              DO  m= 1, NCMM
                  MM= MMLR(m) + Lsr
                  DM(m)= XP*(SM(MM)+ SM(MM+1)+ SM(MM+2)+ SM(MM+3)
     1                                                     + SM(MM+4))
                  IF(KDER.GT.0) THEN
                      DMP(m)= aTT*XP*SM(m-1)
                      IF(KDER.GT.1)DMPP(m)= aTT*aTT*XP*(SM(m-2)-SM(m-1))
                      ENDIF
                  ENDDO
            ENDIF
          ENDIF
c
      IF(IDSTT.GT.0) THEN
c=======================================================================
c** For Douketis-Scoles-Marchetti-Zen-Thakkar type damping function ...
c=======================================================================
          IF((IDF.LT.-4).OR.(IDF.GT.0)) THEN
              WRITE(6,600) IDSTT,IDF
c             STOP
              ENDIF
          IF(FIRST.EQ.1) THEN
              DO m= 1, 20
                  DO  IDFF= -2,0
                      bpm(m,IDFF)= bDS(IDFF)/DFLOAT(m)
                      cpm(m,IDFF)= cDS(IDFF)/DSQRT(DFLOAT(m))
                      ENDDO
                  ENDDO
              FIRST= 0
              ENDIF
          br= rhoAB*r
          DO m= 1, NCMM
              MM= MMLR(m)
              XP= DEXP(-(bpm(MM,IDF) + cpm(MM,IDF)*br)*br)
              YP= 1.d0 - XP
              ZK= MM-1.d0
              DM(m)= YP**(MM-1)
c... Actually ...  DM(m)= YP**(MM + IDF/2)  :  set it up this way to
c   avoid taking exponential of a logarithm for fractional powers (slow)
              IF(IDF.EQ.-4) THEN
                  ZK= ZK- 1.d0
                  DM(m)= DM(m)/YP
                  ENDIF
              IF(IDF.EQ.-3) THEN
                  ZK= ZK- 0.5d0
                  DM(m)= DM(m)/DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.-1) THEN
                  ZK= ZK+ 0.5d0
                  DM(m)= DM(m)*DSQRT(YP)
                  ENDIF
              IF(IDF.EQ.0) THEN
                  ZK= MM
                  DM(m)= DM(m)*YP
                  ENDIF
              IF(KDER.GT.0) THEN
                  TK= bpm(MM,IDF) + 2.d0*cpm(MM,IDF)*br
                  DMP(m) = ZK*XP*rhoAB*TK*DM(m)/YP
                  IF(KDER.GT.1) THEN
c ... if desired ... calculate second derivative [for DELR case] {check this!}
                      DMPP(m)= (ZK-1.d0)*XP*TK*DMP(m)/YP
     1               - DMP(m)*TK + DMP(m)*2.d0*cpm(MM,IDF)*rhoAB**2/TK
                      ENDIF
                  ENDIF
              ENDDO
          ENDIF
      RETURN
  600 FORMAT(/,' *** ERROR ***  For  IDSTT=',i3,'   IDF=',i3,'  no dampi
     1ng function is defined')
  602 FORMAT( /,' ***ERROR ***  rhoAB=', F7.4,'  yields an invalid Dampi
     1ng Function definition')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
