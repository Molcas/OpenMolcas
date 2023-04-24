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
!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
!** Subroutine to calculates quasibound level tunneling lifetime/width
!** For relevant theory see Le Roy & Liu [J.Chem.Phys.69,3622-31(1978)]
!  and Connor & Smith [Mol.Phys. 43, 397 (1981)] and Huang & Le Roy
!  [J.Chem.Phys. 119, 7398 (2003); Erratum, ibid, 127, xxxx (2007)]
!** Final level width calculation from Eq.(4.5) of Connor & Smith.
!------------------ Corrected: 12 March 2007 --------------------------
      SUBROUTINE WIDTHas(KV,JROT,E,EO,DSOC,V,S,RVB,SDRDY,VMX,RMIN,H,    &
     &  BFCT,IWR,ITP1,ITP2,ITP3,INNER,NPP,GAMA)
!++ "WIDTH" calls subroutine "LEVQAD" ++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,IMM,INNER,IRM,ITP1,ITP1P,ITP1P1,ITP2,ITP2M,ITP2M2,     &
     &  ITP2P1,ITP2P2,ITP3,IWR,JROT,KV,KVI,KVO,M,M2,NPP,NN,NST
      REAL*8  ANS1,ANS2,ARG,BFCT,COR,D1,D2,D3,DFI,DSGB,DSGN,DSOC,DWEB,  &
     &  OMEGJC,E,EO,EMSC,EMV,G1,G2,G3,GAMA,GAMALG,H,H2,HBW,HBWB,PI,     &
     &  PMX,RMIN,RMINN,RMX,RT,R1,R2,R3,SM,TAU,TAULG,TI,TUN0,U1,U2,      &
     &  VMAX,VMX,XJ,XX,V(NPP),S(NPP),RVB(NPP),SDRDY(NPP)
      CHARACTER*5 LWELL(2)
      DATA PI/3.141592653589793D0/
      DATA LWELL/'INNER','OUTER'/
      IMM = 0
      PMX = 0
      RMINN= RMIN- H
      H2= H*H
!** First - locate innermost turning point ...
      DO  I= 1,ITP2
          ITP1= I
          IF(V(I).LT.E) GOTO 40
          ENDDO
      GAMA= 0.d0
      GO TO 250
!** ITP1 is first mesh point to right of innermost turning point.
   40 ITP1P= ITP1+ 1
      ITP1P1= ITP1P+ 1
      IRM= ITP1- 1
!** Calculate JWKB tunneling probability from quadrature over barrier
!** (ITP2 is first point inside barrier - as determined in QBOUND)
      ITP2P1= ITP2+ 1
      ITP2P2= ITP2+ 2
!** ITP2M is the last mesh point before the 2-nd turning point.
      ITP2M= ITP2- 1
      ITP2M2= ITP2- 2
      G1= V(ITP2M)- E
      G2= V(ITP2)- E
      G3= V(ITP2P1)- E
      R1= RVB(ITP2M)
      R2= RVB(ITP2)
      R3= RVB(ITP2P1)
      WRITE(6,*) R1,R2,R3 ! make them "referenced"
!** Quadrature over barrier starts here.
      CALL LEVQAD(G1,G2,G3,H,RT,ANS1,ANS2)
      SM= ANS2*SDRDY(ITP2)**2/H
!c    SM= ANS2/H
      IF(G3.LT.0.d0) GO TO 218
      SM= SM+ 0.5d0*DSQRT(G3)*SDRDY(ITP2)**2
      PMX= VMX
      M2= ITP2P2
  204 DO  I=M2,ITP3
          M= I
          G3= V(I)- E
          IF(V(I).GT.PMX) PMX=V(I)
          IF(G3.LT.0.d0) GO TO 210
          SM= SM+ DSQRT(G3)*SDRDY(I)**2
          ENDDO
      IF(V(M).GT.V(M-1)) THEN
          IF(IWR.NE.0) WRITE(6,602) KV,JROT
          GO TO 250
          ENDIF
      RMX= RMINN+ M*H
      U1= DSQRT(G3/(V(M)- DSOC))
      U2= DSQRT((E- DSOC)/(V(M)- DSOC))
      SM= SM- 0.5d0*DSQRT(G3)+ (DLOG((1.d0+U1)/U2)-U1)*RMX*             &
     &                                             DSQRT(V(M)- DSOC)/H
      XJ= (DSQRT(1.d0+ 4.d0*(V(M)-DSOC)*(RMX/H)**2)- 1.d0)*0.5d0
      IF(IWR.NE.0) WRITE(6,603) JROT,EO,XJ,RMX
      GO TO 218
  210 IF(M.LT.ITP3) THEN
!** If encounter a double-humped barrier, take care here.
          IF(IWR.NE.0) WRITE(6,609) KV,JROT,EO,M
          KVO= 0
          DSGN= DSIGN(1.d0,S(M-1))
!** Find the effective quantum number for the outer well
          DO  I= M,ITP3
              DSGB= DSGN
              DSGN= DSIGN(1.d0,S(I))
              IF((DSGN*DSGB).LT.0.d0) KVO=KVO+1
              ENDDO
          KVI= KV- KVO
          IF(INNER.EQ.0) THEN
!** For levels of outer well, get correct width by changing ITP1
              ITP1= M
              IF(IWR.GT.0) WRITE(6,610) KVO,LWELL(2)
              GO TO 40
              ENDIF
          IF(IWR.GT.0) WRITE(6,610) KVI,LWELL(1)
!** For "inner-well" levels, locate outer barrier
          DO  I= M,ITP3
              M2= I
              G3= V(I)- E
              IF(G3.GE.0.d0) GO TO 204
              ENDDO
          GO TO 218
          ENDIF
      G1= V(M) - E
      G2= V(M-1)- E
      G3= V(M-2)- E
      R1= RVB(M)
      R2= RVB(M-1)
      R3= RVB(M-2)
      CALL LEVQAD(G1,G2,G3,H,RT,ANS1,ANS2)
      SM= SM- 0.5d0*DSQRT(G3)*SDRDY(M-2)**2 -DSQRT(G2)*SDRDY(M-1)       &
     &   + ANS2*SDRDY(M)**2/H                                           &
!c   1   + ANS2/H
     &   + ANS2/H
  218 EMSC= -SM/PI
      IF(INNER.GT.0) VMX= PMX
      VMAX= VMX/BFCT
!** Tunneling factors calculated here ** TUN0 is simple WKB result
!  as in Child's eqs.(57c) & (59).
! .....  EPSRJ= -2.* PI* EMSC
      TUN0= 0.5d0*DEXP(2.d0*PI*EMSC)
! ... for permeability calculate Connor-Smith's Eq.(3.7) \omega=OMEGJC
      OMEGJC= DSQRT(1.d0+ 2.d0*TUN0) - 1.d0
! ... alternate calculation to give better precision for small TUN0
      IF(TUN0.LT.1.d-5) OMEGJC= TUN0*(1.d0-0.5d0*TUN0*(1.d0-TUN0))
      OMEGJC= 4.d0*OMEGJC/(OMEGJC + 2.d0)
!** Quadrature for JWKB calculation of vibrational spacing in well HBW
      D1= E- V(IRM)
      D2= E- V(ITP1)
      D3= E- V(ITP1P)
      R1= RVB(IRM)
      R2= RVB(ITP1)
      R3= RVB(ITP1P)
      CALL LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
      SM= ANS1*SDRDY(ITP1)**2/H
!c    SM= ANS1/H
      IF(D3.LT.0.d0) GO TO 228
      SM= SM+ 0.5d0/DSQRT(D3)
      DO  I= ITP1P1,ITP2M2
          IMM= I
          EMV= E- V(I)
          IF(EMV.LT.0.d0) GO TO 222
          SM= SM+ SDRDY(I)**2/DSQRT(EMV)
          ENDDO
      D3= E- V(ITP2M2)
      D2= E- V(ITP2M)
      D1= E- V(ITP2)
      R1= RVB(ITP2)
      R2= RVB(ITP2M)
      R3= RVB(ITP2M2)
      GO TO 226
!** If encounter a double-minimum well, take care here.
  222 D1= EMV
      D2= E- V(IMM-1)
      D3= E- V(IMM-2)
      R1= RVB(IMM)
      R2= RVB(IMM-1)
      R3= RVB(IMM-2)
      IF(IWR.NE.0) WRITE(6,605) KV,JROT,EO
!c226 CALL LEVQAD(D1,D2,D3,R1,R2,R3,H,RT,ANS1,ANS2)
  226 CALL LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
      SM= SM-0.5d0*SDRDY(IMM-2)/DSQRT(D3) + ANS1*SDRDY(IMM-1)**2/H
!c    SM= SM-0.5d0*SDRDY(IMM-2)/DSQRT(D3) + ANS1/H
!** Get HBW in same energy units (1/cm) associated with BFCT
  228 HBW=2.d0*PI/(BFCT*SM)
!** HBW fix up suggested by Child uses his eqs.(48)&(62) for HBW
!** Derivative of complex gamma function argument calculated as
!  per eq.(6.1.27) in Abramowitz and Stegun.
      NST= INT(DABS(EMSC)*1.D2)
      NST= MAX0(NST,4)
      ARG= -1.963510026021423d0
      DO  I= 0,NST
          NN= I
          XX= I + 0.5d0
          TI= 1.d0/(XX*((XX/EMSC)**2 + 1.d0))
          ARG= ARG+TI
          IF(DABS(TI).LT.1.D-10) GO TO 233
          ENDDO
! ... and use continuum approximation for tail of summation (???)
  233 COR= 0.5d0*(EMSC/(NN+1.d0))**2
      ARG= ARG+ COR- COR**2
!** Now use WKL's Weber fx. approx for (?) derivative of barrier integral ..
      DWEB= (EO-VMAX)*BFCT/(H2*EMSC)
      DFI= (DLOG(DABS(EMSC)) - ARG)*BFCT/(H2*DWEB)
      HBWB= 1.d0/(1.d0/HBW + DFI/(2.d0*PI))
!** Width from formula (4.5) of  Connor & Smith, Mol.Phys.43,397(1981)
! [neglect time delay integral past barrier in their Eq.(4.16)].
      IF(EMSC.GT.-25.D0) THEN
          GAMA= (HBWB/(2.d0*PI))* OMEGJC
          TAU= 0.D0
          IF(GAMA.GT.1.D-60) TAU= 5.308837457D-12/GAMA
!** GAM0 = TUN0*HBW/PI  is the simple WKB width GAMMA(0) discussed by
!  Le Roy & Liu in J.C.P.69,3622(1978).
          IF(IWR.GT.0) WRITE(6,601) TAU,GAMA,HBWB,VMAX
        ELSE
          GAMALG= DLOG10(HBWB/(2.d0*PI))+2.d0*PI*EMSC/2.302585093D0
          TAULG= DLOG10(5.308837457D-12)-GAMALG
          IF(IWR.GT.0) WRITE(6,611) TAULG,GAMALG,HBWB,VMAX
        ENDIF
  250 RETURN
  601 FORMAT('    Lifetime=',1PD10.3,'(s)   Width=',D10.3,'   dG/dv=',  &
     & 0PF7.2,'   V(max)=',F9.2)
  602 FORMAT(' *** WARNING ***  For   v =',I3,'   J =',I3,'   cannot cal&
     &culate width since barrier maximum beyond range')
  603 FORMAT(' *** For  J=',I3,'  E=',F9.2,'  R(3-rd) beyond range so ap&
     &prox. tunneling calc. uses'/8X,'pure centrifugal potential with  J&
     &(app)=',F7.2,'  for  R > R(max)=',F7.2)
  605 FORMAT(' **** CAUTION *** Width estimate only qualitative, as have&
     & a double-minimum well for   E(v=',I3,', J=',I3,')=',F15.7/15X,   &
     & 'a more stable result may be obtained by searching for the quasib&
     &ound levels using option: INNER > 0 .')
  609 FORMAT(' *** CAUTION - Permeability estimate not exact as have a d&
     &ouble-humped barrier:  E(v=',I3,', J=',I3,') =',G15.8,I6)
  610 FORMAT(16X,'(NOTE: this has the node count of a   v=',I3,2X,A5,   &
     & '-well level')
  611 FORMAT(12X,'Log10(lifetime/sec)=',F10.5,' ;   Log10(width/cm-1)=',&
     & F10.5,'   Spacing=',G12.5,'   V(max)=',G14.7,'(cm-1)')
      END
