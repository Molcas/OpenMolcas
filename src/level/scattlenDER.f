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
c***** Subroutine SCATTLEN, last updated 30 April 2011 ****
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** SCATTLEN solves the radial Schrodinger equation in dimensionless
c  form  d^2WF/dy^2 = - [(VLIM-V(R))*(r')^2 - F(y)]*WF(y) ,  where WF(I)
c  is the wave function,  y  the reduced radial vble. y_p(r), and  VLIM
c  the energy at the potential asymptote, to determine the scattering
c  length  SL.
c** Integrate by Numerov method over NPP mesh points with increment
c  H=YH across range from YMIN  to  YMAX= 1. Then uses interpoaltion
c  over last 5 \phi(y) values to the \phi'(y=1) in order to generate
c  SL from log-derivative equation.  After completing the calculation
c  using mesh YH, repeat it with mesh 2*YH and thn use Richardson
c  extrapolation (RE) to improve the final result.  Scheme used requires
c  PRV= 1.d0.
c** Input potential asymptote VLIM has have units (cm-1).
c** On entry, the input potential V(I) must include the centrifugal
c  term and the factor:  'BFCT'=2*mu*(2*pi*YH/hbar)**2  (1/cm-1)
c   = ZMU[u]*YH[Angst]**2/16.85762920 (1/cm-1)  which is also
c  (internally) incorporated into VLIM.  Note that inclusion of the
c  squared integration increment YH**2 saves arithmetic in the
c  innermost loop of the algorithm.
c-----------------------------------------------------------------------
      SUBROUTINE SCATTLEN(JROT,SL,VLIM,V,WF,BFCT,YMIN,YH,NPP,CNN,NCN,
     1                            IWR,LPRWF)
!     USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
      USE LEVEL_COMMON
c-----------------------------------------------------------------------
c** Output scattering length SL [Angst] normalized wave function WF(I)
c  and range, NBEG .le. I .le. NEND  over which WF(I) is defined. Define
c  WF(I)=0  outside the range (NBEG,NEND), which is defined by requiring
c  abs(WF(I)) < RATST=1.D-9  outside.
c** If(LPRWF.gt.0) print [WRITE(6,xxx)] wavefx WF(I) every LPRWF-th point.
c* If(LPRWF.lt.0) every |LPRWF|-th point of the wave function to Channel
c      10 starting at R(NBEG)
c** If(IWR.ne.0) print error & warning descriptions
c  If (IWR.ge.2) also show end-of-range wave function amplitudes
c  If (IWR.ge.3) print also intermediate trial eigenvalues, etc.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!!
      INTEGER NDIMR
!     PARAMETER (NDIMR= 200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive.
!     PARAMETER (NDIMR= 131074)
!     REAL*8 PRV,ARV,RVB(NDIMR),YVB(NDIMR),DRDY2(NDIMR),FAS(NDIMR),
!    1                                         SDRDY(NDIMR),VBZ(NDIMR)
      REAL*8 PRV,ARV
!     REAL*8, ALLOCATABLE :: RVB(:),YVB(:),DRDY2(:),FAS(:),
!    1                                         SDRDY(:),VBZ(:)
      COMMON /BLKAS/PRV,ARV!,RVB,YVB,DRDY2,SDRDY,FAS,VBZ
c!!
      INTEGER  I,ITP1,ITP1P,IWR,J,JPSIQ,JROT,LPRWF,
     1  LNPT0,NCN,NPP,NBEG,NBEG2,NPR,NP2,NODE,NNH
      REAL*8  BFCT,DSOC,GI,GN,HT,RATIN,RATST,SB,SI,SL,SL2,SLcor,
     x  sumSL,C4BAR,
     1  YH,RINC,YMIN,YMINN,RSTT,WF(NPP),V(NPP),VLIM,Y1,Y2,Y3,
     2  GB,
     3  CNN,
     4  PHIp1,PHIp2,PHIp3,PHIp4,Z4,WF0,WF1,WF2,WF3,WF4,
     5  sumVV

      DATA RATST/1.D-9/,NP2/2/,LNPT0/0/
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SAVE NP2,LNPT0
      NDIMR= 131074
!     CALL MMA_ALLOCATE(RVB,NDIMR,LABEL='RVB')
!     CALL MMA_ALLOCATE(YVB,NDIMR,LABEL='YVB')
!     CALL MMA_ALLOCATE(DRDY2,NDIMR,LABEL='DRDY2')
!     CALL MMA_ALLOCATE(FAS,NDIMR,LABEL='FAS')
!     CALL MMA_ALLOCATE(SDRDY,NDIMR,LABEL='SDRDY')
!     CALL MMA_ALLOCATE(VBZ,NDIMR,LABEL='VBZ')
      WRITE(6,*) LNPT0,NP2,NDIMR !Make sure they are used and referenced
      WF4=0
      IF(DABS(PRV-1.d0).GT.0.d0) THEN
c** Scattering length calculation assumes  PRV=1  s.th.  FAS= 0.0
          WRITE(6,620) PRV
          SL= 0.d0
          RETURN
          ENDIF
      Z4= 4.d0
      YMINN= YMIN-YH
      HT= 1.d0/12.D+0
      DSOC= VLIM*BFCT
      RATIN= 0.d0
      NBEG= 1
      C4BAR= 0.d0
      IF(NCN.EQ.4) C4BAR= BFCT*CNN/(2.d0*ARV)**2
c** Begin by checking that Numerov is stable at innermost end of range ...
   10 GN= V(NBEG) - DSOC*DRDY2(NBEG)
      IF(GN.GT.10.D0) THEN
c** If potential has [V(i)-E] so high that H is (locally) too ;arge,
c  then shift inner starting point outward.
          NBEG= NBEG+1
          IF(NBEG.LT.NPP) GO TO 10
          IF(IWR.NE.0) WRITE(6,600)
          GO TO 999
          ENDIF
      IF(IWR.NE.0) THEN
          IF(NBEG.GT.1) WRITE(6,602) JROT,NBEG,YVB(NBEG)
          IF(GN.LE.0.d0) WRITE(6,604) JROT,NBEG,V(NBEG)/BFCT
          ENDIF
      NNH= (NPP-NBEG)/2
      IF((NPP-NBEG).GT.(2*NNH)) THEN
c** If necessary, adjust NBEG by 1 to ensure interval has an even
c   no. mesh points in order to simplify RE correction step
          NBEG= NBEG+1
          GN= V(NBEG) - DSOC*DRDY2(NBEG)
          ENDIF
c** Initialize outward wave function with a node:  WF(NBEG) = 0.
      SB= 0.d0
      SI= 1.d0
      GI= V(NBEG+1) - DSOC*DRDY2(NBEG+1)
      Y1= SB*(1.d0- HT*GN)
      Y2= SI*(1.d0- HT*GI)
      WF(NBEG)= SB
      WF(NBEG+1)= SI
      NODE= 0
c     sumSL= SI*(GI/SDRDY(NBEG+1))
c    1                      *(1.D0 + YVB(NBEG+1))/(1.D0 - YVB(NBEG+1))
      sumSL= SI*GI*(1.D0 + YVB(NBEG+1))
c** Actual outward integration loops start here
      DO  I= NBEG+2,NPP
          Y3= Y2+Y2-Y1+GI*SI
          GI= V(I) - DSOC*DRDY2(I)
          SI= Y3/(1.d0- HT*GI)
          WF(I)= SI
cc        sumSL= sumSL+ (GI*SI/SDRDY(I)) *(1.d0+YVB(I))/(1.d0 - YVB(I))
          sumSL= sumSL+ GI*SI*(1.d0+YVB(I))
          IF(DABS(SI).GE.1.D+36) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically forbidden
c  region where  V(I) .gt. E
              SI= 1.d0/SI
              sumSL= sumSL*SI
              DO  J= NBEG,I
                  WF(J)= WF(J)*SI
                  ENDDO
              Y2= Y2*SI
              Y3= Y3*SI
              SI= 1.d0
              ENDIF
          ITP1= I
c** Exit from this loop at onset of classically allowed region
          IF(GI.LE.0.d0) GO TO 20
          Y1= Y2
          Y2= Y3
          ENDDO
      IF(IWR.NE.0) WRITE(6,606) JROT,NPP
      GO TO 999
   20 ITP1P= ITP1+1
      DO  I= ITP1P, NPP-1
c** Now - integrate automatically to second-last mesh point ...
          Y1= Y2
          Y2= Y3
          Y3= Y2+Y2-Y1+GI*SI
          GB= GI
          GI= GB ! Make sure GB is "referened".
          GI= V(I) - DSOC*DRDY2(I)
          SB= SI
          SI= Y3/(1.d0- HT*GI)
          sumSL= sumSL+ GI*SI*(1.d0+YVB(I))
c** perform node count ...
          IF(SI*SB.LE.0.d0) THEN
              IF(DABS(SI).GT.0.d0) NODE= NODE+1
              ENDIF
          WF(I)= SI
          ENDDO
c** Finally ... complete integration to very last mesh point at  y= 1,
      Y1= Y2
      Y2= Y3
      Y3= Y2+Y2-Y1+GI*SI
      GB= GI
      IF(NCN.GT.4) GI= 0.d0
      IF(NCN.EQ.4) GI= -C4BAR
      SB= SI
      SI= Y3/(1.d0- HT*GI)
      WF(NPP)= SI
c** Now generate a value for  \phi'(y=1) from the WF values using
c  "Newton's formula for forward interpolation" as described in
c  Sect. 1.4 of K. Smith "Calculation of Atomic Collision Processes"`
      PHIp1= (WF(NPP)- WF(NPP-1))/YH
      PHIp2= PHIp1 + (WF(NPP) -2.d0*WF(NPP-1) + WF(NPP-2))/(2.d0*YH)
      PHIp3= PHIp2 + (WF(NPP) - 3.d0*WF(NPP-1) + 3.d0*WF(NPP-2)
     1                                      - WF(NPP-3))/(3.d0*YH)
      PHIp4= PHIp3 + (WF(NPP) - 4.d0*WF(NPP-1) + 6.d0*WF(NPP-2)
     1                     - 4.d0*WF(NPP-3) + WF(NPP-4))/(4.d0*YH)
      SL= ARV*(2.d0*PHIp4/WF(NPP) - 1.d0)
      WRITE(6,608) SL,PHIp4/SI,PHIp1,PHIp2,PHIp3,PHIp4
cc=====================================================================

c** If desired, calculate partial derivatives of scattering length
c  w.r.t. parameters.
c** DF*H  is the integral of  (WF(I))**2 dR
c!!   IF(NPARM.GT.0) THEN
c!!       DO  J= 1, NPARM
c!!           DADPARM(J)= 0.d0
c!!           ENDDO
c!!       DO  I= NBEG,NPP
c!!           DF= DRDY2(I)*WF(I)**2
c!!           DO  J= 1,NPARM
c!!               DADPARM(J)= DADPARM(J) + DF*DVDP(I,J)
c!!               ENDDO
c!!           ENDDO
c!!       DO  J= 1, NPARM
c!!           DADPARM(J)= DADPARM(J)*YH
c!!           ENDDO
c!!
      IF((DABS(RATIN).GT.RATST).AND.(YMIN.GT.0.d0))
     1                                         WRITE(6,614) JROT,RATIN
      IF(LPRWF.LT.0) THEN
c** If desired, write every |LPRWF|-th point of the wave function
c  to a file on channel-10, starting at the NBEG-th mesh point.
          JPSIQ= -LPRWF
          NPR= 1+(NPP-NBEG)/JPSIQ
          RINC= YH*JPSIQ
          RSTT= YMINN+NBEG*YH
c** Write every JPSIQ-th point of the wave function, beginning at mesh
c  point NBEG & distance RSTT where
c  the NPR values written separated by mesh step RINC=JPSIQ*YH
          WRITE(10,701) JROT,NPR,RSTT,RINC,NBEG,JPSIQ
          WRITE(10,702) (YVB(I),WF(I),I=NBEG,NPP,JPSIQ)
          ENDIF
c
c** Now ... re-do SL calculation with twice the step size to allow
c  Richardson Extraoplation correction extimation ...
c** Initialize outward wave function with a node:  WF(NBEG) = 0.
      SB= 0.d0
      SI= 1.d0
      GN= Z4*(V(NBEG) - DSOC*DRDY2(NBEG))
      GI= Z4*(V(NBEG+2) - DSOC*DRDY2(NBEG+2))
      Y1= SB*(1.d0- HT*GN)
      Y2= SI*(1.d0- HT*GI)
      WF1= SB
      WF0= SI
      NBEG2= NBEG+2
      NBEG=NBEG2-2 ! Make sure NBEG2 is referenced
c     sumSL= SI*(GI/SDRDY(NBEG+1))
c    1                      *(1.D0 + YVB(NBEG+1))/(1.D0 - YVB(NBEG+1))
      sumSL= SI*GI*(1.D0 + YVB(NBEG+2))
c** Actual outward integration loops start here
      DO  I= NBEG+4,NPP,2
          WF1= WF0
          Y3= Y2+Y2-Y1+GI*SI
          GI= (V(I) - DSOC*DRDY2(I))*Z4
          SI= Y3/(1.d0- HT*GI)
          WF0= SI
cc        sumSL= sumSL+ (GI*SI/SDRDY(I)) *(1.d0+YVB(I))/(1.d0 - YVB(I))
          sumSL= sumSL+ GI*SI*(1.d0+YVB(I))
          IF(DABS(SI).GE.1.D+36) THEN
c** Renormalize to prevent overflow of  WF(I)  in classically forbidden
c  region where  V(I) .gt. E
              SI= 1.d0/SI
              WF1= WF1*SI
              sumSL= sumSL*SI
              Y2= Y2*SI
              Y3= Y3*SI
              SI= 1.d0
              WF0= SI
              ENDIF
          ITP1= I
c** Exit from this loop at onset of classically allowed region
          IF(GI.LE.0.d0) GO TO 40
          Y1= Y2
          Y2= Y3
          ENDDO
      IF(IWR.NE.0) WRITE(6,606) JROT,NPP
      GO TO 999
   40 ITP1P= ITP1+2
      WF2= WF1
      WF3= WF2
      DO  I= ITP1P, NPP, 2
c** Now - integrate automatically to second-last mesh point ...
          Y1= Y2
          Y2= Y3
          Y3= Y2+Y2-Y1+GI*SI
          GB= GI
          GI= (V(I) - DSOC*DRDY2(I))*Z4
          IF(I.EQ.NPP) THEN
              IF(NCN.GT.4) GI= 0.d0
              IF(NCN.EQ.4) GI= -C4BAR*Z4
              ENDIF
ccc   IF(NCN.GT.4) GI= 0.d0
ccc... HEY ... RJ should go & figure out how to treat the C4 case!
ccc  C4bar = BFCT*C4/(4*ARV**2)   ???
          WF4= WF3
          WF3= WF2
          WF2= WF1
          WF1= WF0
          SI= Y3/(1.d0- HT*GI)
          WF0= SI
          sumSL= sumSL+ GI*SI*(1.d0+YVB(I))
          ENDDO
c** Now generate a value for  \phi'(y=1) from the WF values using
c  "Newton's formula for forward interpolation" as described in
c  Sect. 1.4 of K. Smith "Calculation of Atomic Collision Processes"`
      PHIp1= (WF0- WF1)/(2.d0*YH)
      PHIp2= PHIp1 + (WF0 - 2.d0*WF1 + WF2)/(4.d0*YH)
      PHIp3= PHIp2 + (WF0 - 3.d0*WF1 + 3.d0*WF2 - WF3)/(6.d0*YH)
      PHIp4= PHIp3 + (WF0 - 4.d0*WF1 + 6.d0*WF2 - 4.d0*WF3 + WF4)/
     1                                                   (8.d0*YH)
c...  SL2  is scattering length associated with mesh of  2*YH
      SL2= ARV*(2.d0*PHIp4/WF0 - 1.d0)
      WRITE(6,608) SL2,PHIp4/SI,PHIp1,PHIp2,PHIp3,PHIp4
c** Finally - user Ricardson expraolation of results for mesh  YH  and
c   2*YH  to obtain final optimum  SL estimate!
      SLcor= SL + (SL-SL2)/15.d0
      WRITE(6,610)  YH, SLCOR, SL2,SL
      WRITE(8,610)  YH, SLCOR, SL2,SL
cc    WRITE(6,612) NODE-1
      SL= SLcor
c** Now ... use second-last mesh point to normalize wavefunction to
c  correspond to asymptotic normalization  \psi(r) \sim r .
      SI= (RVB(NPP-1)-SL)/(WF(NPP-1)*SDRDY(NPP-1))
      SUMVV= 0.d0
      DO  I= NBEG, NPP-1
          WF(I)= WF(I)*SI
c ... and calculate expectation values of  V(r)  in cm-1
          SUMVV= SUMVV+ DRDY2(I)*V(I)*WF(I)**2
          ENDDO
      SUMVV= SUMVV/YH
      WRITE(6,616)  SUMVV, BFCT
  616 FORMAT(' Expectation value of  V(r) is:', 1PD17.8,'   BFCT=',
     1  D17.8)

      WRITE(6,612) NODE-1
!     CALL MMA_DEALLOCATE(RVB)
!     CALL MMA_DEALLOCATE(YVB)
!     CALL MMA_DEALLOCATE(DRDY2)
!     CALL MMA_DEALLOCATE(FAS)
!     CALL MMA_DEALLOCATE(SDRDY)
!     CALL MMA_DEALLOCATE(VBZ)
      RETURN
c** ERROR condition if  E.gt.V(R)  at outer end of integration range.
c** Return in error mode
  999 JROT= -1
      RETURN
  600 FORMAT(/' *** ERROR in potential array ... V(I) everywhere',
     1 ' too big to integrate with given  increment')
  602 FORMAT(' *** For  J=',I3,"  integration can't start till past"/
     1  23x,'mesh point',I5,' (yp=',0pf6.2,'),  so YMIN smaller than nee
     2ded')
  604 FORMAT('   NOTE:  for  J=',I3,'   V(',i3,')=',F12.4,' .LE. 0.0')
  606 FORMAT(/' *** ERROR *** for   J =',I3,'  Innermost turning point n
     1ot found by   M = MSAVE =',I5)
  608 FORMAT(/' Calculate  SL=',1PD21.13,'   log-derivative(y=1)=',
     1 D20.12/'     with slope convergence:',D21.13/(28x,D21.13))
  610 FORMAT(/' YH=',f10.7,'  gives  SL(RE)=',1PD21.13,':  SL2=',
     1  D21.13/55x,'SL=',D21.13)
  612 FORMAT(/' Last bound level of this potential is   v=',i3////)
  614 FORMAT(' *** CAUTION *** For  J=',I3,'   WF(first)/WF(Max)=',D9.2,
     1  '  suggests  YMIN  may be too large')
  620 FORMAT(/' *** ERROR in scattlen ***  Input  PRV=',F7.3,
     1   '  .NE. 1')
  701 FORMAT(/2x,'For   J=',I3,',  wave function at',I6,' points.'/
     1  7x,'R(1-st)=',F12.8,'   mesh=',F12.8,'   NBEG=',I4,
     2  '   |LPRWF|=',I3)
  702 FORMAT((1X,4(0Pf9.5,1PD13.5)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
