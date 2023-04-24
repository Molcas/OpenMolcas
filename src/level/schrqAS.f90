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
!*****  subroutine SCHRQ, last updated  21 August 2007 *****
!   !!!!   Form modified for handling Stolyarov radial variable  !!!!!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** SCHRQ solves radial Schrodinger equation in dimensionless form
!  d^2WF/dy^2 = - [(E-V(R))*BFCT*(r')^2 - F(y)]*WF(R) ,  where WF(I) is
!  the wave function and  y  the reduced radial vble.  y_p(r).
!** Integrate by Numerov method over NPP mesh points with increment
!  H=YH across range beginning at YMIN .
!** Input trial energy EO, eigenvalue convergence criterion EEPS
!  potential asymptote VLIM, and all returned energies (EO, GAMA & VMAX)
!  have units (cm-1).
!** On entry, the input potential V(I) must include the centrifugal
!  term, the factor:  'BFCT'=2*mu*(YH/hbar)**2 (1/cm-1) as well as the
!  Stolyarov conversion factors (r')^2 and F(y).
!  BFCT is also internally incorporated into EO, VLIM & EEPS.
!* Note that these reduced quantities (& the internal eigenvalue E)
!  contain a factor of the squared integration increment  YH**2 .
!  This saves arithmetic work in the innermost loop of the algorithm.
!** For energy in (cm-1), BFCT=ZMU(u)*H(Angst)**2/16.85762920 (1/cm-1)
!** INNODE > 0  specifies that wavefx. initiates at YMIN with a node
!     (normal default case);  INNODE.le.0  specifies  zero slope  at
!     YMIN (for finding symmetric eigenfunctions of symmetric potential
!     with potential mid-point @ YMIN).
!** INNER specifies wave function matching condition: INNER = 0  makes
!     matching of inward & outward solutions occur at outermost turning
!     point;  INNER > 0 makes matching occur at innermost turning point.
! * Normally use  INNER=0 ,  but to find inner-well levels of double
!     minimum potential, set  INNER > 0 .
!-----------------------------------------------------------------------
      SUBROUTINE SCHRQas(KV,JROT,EO,GAMA,VMAX,VLIM,V,WF,BFCT,EEPS,YMIN, &
     &                        YH,NPP,NBEG,NEND,INNODE,INNER,IWR,LPRWF)
!     USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
      USE LEVEL_COMMON
!-----------------------------------------------------------------------
!** Output vibrational quantum number KV, eigenvalue EO, normalized
!  wave function WF(I), and range, NBEG .le. I .le. NEND  over
!  which WF(I) is defined. *** Have set  WF(I)=0  outside this range.
!* (NBEG,NEND), defined by requiring  abs(WF(I)) < RATST=1.D-9  outside.
!** If(LPRWF.NE.0) write every LPRWF-th value of wavefunction WF(I) to
!   a file on channel-10 (i.e., WRITE(10,XXX)), starting at YVB(NBEG)
!   with step size  |LPRWF|*YH.
!** For energies above the potential asymptote VLIM, locate quasibound
!  levels using Airy function boundary condition and return the level
!  width GAMA and barrier height VMAX, as well as EO.
!** ERROR condition on return is  KV < 0 ; usually KV=-1, but return
!  KV=-2 if error appears to arise from too low trial energy.
!** If(IWR.ne.0) print error & warning descriptions
!  If (IWR.gt.0) also print final eigenvalues & node count.
!  If (IWR.ge.2) also show end-of-range wave function amplitudes
!  If (IWR.ge.3) print also intermediate trial eigenvalues, etc.
!** If input KV.ge.998 , tries to find highest bound level, and
!  trial energy should be only slightly less than VLIM.
!-----------------------------------------------------------------------
!++ "SCHRQ" calls subroutines "QBOUND" and "WIDTH", and the latter
!++ calls "LEVQAD" .
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!
      IMPLICIT NONE
!     INTEGER NDIMR
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
!!!
      INTEGER  I,IBEGIN,ICOR,INNODE,INNER,IT,ITER,ITP1,ITP1P,           &
     &  ITP2,ITP3,IWR,J,J1,J2,JPSIQ,JQTST,JROT,KKV,KV,KVIN,LPRWF,M,     &
     &  MS,MSAVE,NPP,NBEG,NDN,NEND,NPR
!
      REAL*8  BFCT,DE,DEP,DEPRN,DF,DOLD,DSOC,E,EEPS,EO,EPS,F,GAMA,      &
     &  GI,GB,H,HT,PROD,PPROD,RATIN,RATOUT,RATST,REND,YH,YMIN,          &
     &  YMINN,RR,SB,SI,SN,SRTGI,SRTGB,SM,VLIM,VMAX,VMX,VPR,WKBTST,XEND, &
     &  XPR,XPW,DXPW,Y1,Y2,Y3,YIN,YM,YOUT,WF(NPP),V(NPP)
      DATA RATST/1.D-9/,XPW/23.03d0/
      DATA NDN/10/
!     NDIMR = 131074
!     CALL MMA_ALLOCATE(RVB,NDIMR,LABEL='RVB')
!     CALL MMA_ALLOCATE(YVB,NDIMR,LABEL='YVB')
!     CALL MMA_ALLOCATE(DRDY2,NDIMR,LABEL='DRDY2')
!     CALL MMA_ALLOCATE(FAS,NDIMR,LABEL='FAS')
!     CALL MMA_ALLOCATE(SDRDY,NDIMR,LABEL='SDRDY')
!     CALL MMA_ALLOCATE(VBZ,NDIMR,LABEL='VBZ')
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING
!     WRITE(6,*) 'After entering schrq.f we have:'
!     WRITE(6,*) 'KV=',KV
!     WRITE(6,*) 'JROT=',JROT
!     WRITE(6,*) 'EO=',EO
!     WRITE(6,*) 'GAMA=',GAMA
!     WRITE(6,*) 'VMAX=',VMAX
!     WRITE(6,*) 'VLIM=',VLIM
!     DO I=1,3
!      WRITE(6,*) 'V=',V(I)
!      WRITE(6,*) 'WF=',WF(I)
!     ENDDO
!     WRITE(6,*) 'BFCT=',BFCT
!     WRITE(6,*) 'EEPS=',EEPS
!     WRITE(6,*) 'YMIN=',YMIN
!     WRITE(6,*) 'YH=',YH
!     WRITE(6,*) 'NPP=',NPP
!     WRITE(6,*) 'NBEG=',NBEG
!     WRITE(6,*) 'NEND=',NEND
!     WRITE(6,*) 'INNODE=',INNODE
!     WRITE(6,*) 'INNER=',INNER
!     WRITE(6,*) 'IWR=',IWR
!     WRITE(6,*) 'LPRWF=',LPRWF
      YOUT = 0
      YM = 0
      YIN = 0
      ITP1P = 0
      DXPW= XPW/NDN
      ICOR= 0
      KVIN= KV
      KV= -1
      YMINN= YMIN-YH
      GAMA= 0.d0
      VMAX= VLIM
      VMX= VMAX*BFCT
      H= YH
      HT= 1.d0/12.D+0
      E= EO*BFCT
      EPS= EEPS*BFCT
      DSOC= VLIM*BFCT
      DE= 0.d0
      RATIN= 0.d0
      RATOUT= 0.d0
      IF(IWR.GT.2) THEN
          IF(KVIN.GE.998) THEN
              WRITE(6,610) EO
          ELSE
              WRITE(6,601) KVIN,JROT,EO,INNER
          ENDIF
          WRITE(6,602)
      ENDIF
      NEND= NPP
      GI= 0.d0
      GB= 0.d0
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
!     WRITE(6,*) 'NEND=',NEND
!     WRITE(6,*) 'V(1)=',V(1)
!     WRITE(6,*) 'V(NEND)=',V(NEND)
!     WRITE(6,*) 'E=',E
!     WRITE(6,*) 'DSOC=',DSOC
      JQTST = 0
      WRITE(6,*) 'JQTST=',JQTST ! Make sure it is "referenced"
!** Start iterative loop; try to converge for up to 15 iterations.
! Actually allow only 10 because garble "finds" v=10 with 12 iterations.
      DO 90 IT= 1,10
          ITER= IT
! OPTIONALLY write when debugging:
!         WRITE(6,*) 'INNER=',INNER,'If >0, GO TO 50'
          IF(INNER.GT.0) GO TO 50
   10     IF(E.GT.DSOC) THEN
!** For quasibound l,vels, initialize wave function in "QBOUND"
          ENDIF
          IF(ITER.LE.2) THEN
!** For  E < DSOC  begin inward integration by using JWKB to estimate
!  optimum (minimum) inward starting point which will still give
!  RATOUT < RATST = exp(-XPW) (ca. 1.d-9) [not needed after 1st 2 ITER]
              NEND= NPP - 1
              GB= VBZ(NEND) - E
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
!             WRITE(6,*) 'VBZ(NEND)=',VBZ(NEND)
!             WRITE(6,*) 'GB=',GB
! ... first do rough inward search for outermost turning point
              DO  M= NEND-NDN,1,-NDN
                  ITP2= M
                  GI= VBZ(M) - E
!                 WRITE(6,*) 'VBZ(M)=',VBZ(M)
!                 WRITE(6,*) 'E=',E
!                 WRITE(6,*) 'GI=',GI,'If <= 0, GO TO 12'
                  IF(GI.LE.0.D0.AND.E.LE.0.d0) GO TO 12
!                 IF(GI.LE.0.D0) GO TO 12
                  GB= GI
              ENDDO
              IF(IWR.NE.0) WRITE(6,611) JROT,EO
              GO TO 999
   12         SM= GB/(GI-GB)
              SM= 0.5d0*(1.d0+ SM)*DSQRT(GB)
              ITP2= ITP2+ 2*NDN
! OPTIONALLY write when debugging:
!         WRITE(6,*) 'ITP2=',ITP2,'If >= NEND, GO TO 20'
              IF(ITP2.GE.NEND) GO TO 20
! ... now integrate exponent till JWKB wave fx. would be negligible
              DO  M= ITP2,NPP-1,NDN
                  NEND= M
                  SM= SM + DSQRT(VBZ(M) - E)*SDRDY(M)**2
! OPTIONALLY write when debugging:
!         WRITE(6,*) 'SM=',SM,'If SM > DXPW, GO TO 18'
                  IF(SM.GT.DXPW) GO TO 18
              ENDDO
   18         CONTINUE
          ENDIF
!** Now, checking that {[V-E](r')**2 + FAS} small enuf that Numerov,
!  stable, and if necessary, step inward till  {[V-E](r')**2 - F} < 10
   20     GB= V(NEND) - E*DRDY2(NEND)
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING:
!     WRITE(6,*) 'NEND=',NEND
!     WRITE(6,*) 'V(NEND)=',V(NEND)
!     WRITE(6,*) 'DRDY2(NEND)=',DRDY2(NEND)
!     WRITE(6,*) 'GB=',GB
          IF(GB.GT.10.D0) THEN
!** If potential has [V-E] so high that H is (locally) much too large,
!  then shift outer starting point inward & use WKB starting condition.
!  [extremely unlikely condition w. WKB initialization]
              NEND= NEND-1
! OPTIONALLY write when debugging:
!         WRITE(6,*) 'NEND=',NEND,'If >1, GO TO 20'
              IF(NEND.GT.1) GO TO 20
              IF(IWR.NE.0) WRITE(6,613)
              GO TO 999
              ENDIF
          IF((ITER.LE.1).AND.(IWR.GE.2).AND.(NEND.LT.NPP-1))            &
     &                            WRITE(6,6609) JROT,EO,NEND,YVB(NEND)
          IF(NEND.EQ.NPP-1) THEN
!!! Initialize with node if at end of range  (YMAX= 1)
              NEND= NPP
              SB= 0.d0
              Y1= 0.d0
              SI= 1.d0
              GI= GB
              GO TO 40
          ENDIF
!** For truly bound state initialize wave function as 1-st order WKB
!   solution increasing inward
          GB= V(NEND) - E*DRDY2(NEND)
          GI= V(NEND-1) - E*DRDY2(NEND-1)
          MS= NEND-1
          IF(GI.LT.0.d0) GO TO 998
! Below is an even stronger condition to go to 998. Basically print an error if the level above 0cm=-1
! Comment the three lines below (IF statement) if you want to allow levels above dissociation:
!         WRITE(6,*) 'EO=',EO
!         WRITE(6,*) 'GI=',GI
!         IF(EO.GT.0.d0)  THEN
!             WRITE(6,*) 'Level is not bound!'
!             GO TO 998
!         ENDIF
          SRTGB= DSQRT(VBZ(NEND) - E)
          SRTGI= DSQRT(VBZ(NEND-1) - E)
          SB= 1.d0
          SI= SB*DSQRT(SRTGB/SRTGI)*                                    &
     &        DEXP((SRTGB+SRTGI)*0.5d0*(RVB(NEND)-RVB(NEND-1))/YH)
          IF(SB.GT.SI) THEN
! WOOPS - JWKB gives inward DEcreasing solution, so initialize with node
              IF(IWR.NE.0) WRITE(6,618) JROT,EO,SB/SI
              SI= 1.d0
              SB= 0.d0
              GI= V(NEND-1) - E*DRDY2(NEND-1)
              ENDIF
   40     M= NEND-1
          Y1= (1.d0-HT*GB)*SB
          Y2= (1.d0-HT*GI)*SI
          WF(NEND)= SB
          WF(NEND-1)= SI
          MS= NEND
          IBEGIN= 3
          IF(INNER.GT.0) IBEGIN= ITP1+2
!** Actual inward integration loop starts here
          DO  I= IBEGIN,NEND
              M= M-1
              Y3= Y2+Y2-Y1+GI*SI
              GI= V(M) - E*DRDY2(M)
              SB= SI
              SI= Y3/(1.d0-HT*GI)
              WF(M)= SI
              IF(DABS(SI).GE.1.D+17) THEN
!** Renormalize to prevent overflow of  WF(I)  in classically
!  forbidden region where  (V(I) .gt. E)
                  SI= 1.d0/SI
                  DO  J= M,MS
                      WF(J)= WF(J)*SI
                      ENDDO
!cc               MS= M
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SB= SB*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
!** Test for outermost maximum of wave function.
! ... old matching condition - turning point works OK & is simpler.
!c            IF((INNER.EQ.0).AND.(DABS(SI).LE.DABS(SB))) GO TO 44
!** Test for outer well turning point
!             WRITE(6,*) 'GI=',GI
              IF((INNER.EQ.0).AND.(GI.LT.0.d0)) GO TO 44
              ENDDO
          IF(INNER.EQ.0) THEN
!** Error mode ... inward propagation finds no turning point
              KV= -2
              IF(IWR.NE.0) WRITE(6,616) KV,JROT,EO
              GO TO 999
              ENDIF
!** Scale outer part of wave function before proceding
   44     SI= 1.d0/SI
          YIN= Y1*SI
          MSAVE= M
          RR= YMINN+MSAVE*H
          RATOUT= WF(NEND-1)*SI
          DO  J= MSAVE,NEND
              WF(J)= WF(J)*SI
              ENDDO
          IF(INNER.GT.0) GO TO 70
!-------------------------------------------------------------------
!** Set up to prepare for outward integration **********************
   50     NBEG= 2
          IF(INNODE.LE.0) THEN
!** Option to initialize with zero slope at beginning of the range
              SB= 1.d0
              GB= V(1) - E*DRDY2(1)
              Y1= SB*(1.d0-HT*GB)
              Y2= Y1+GB*SB*0.5d0
              GI= V(2) - E*DRDY2(2)
              SI= Y2/(1.d0-HT*GI)
            ELSE
!** Initialize outward integration with a node at beginning of range
   60         GB= V(NBEG) - E*DRDY2(NBEG)
              IF(GB.GT.10.D0) THEN
!** If potential has [V(i)-E] so high that H is (locally) much too
!  large, then shift inner starting point outward.
                  NBEG= NBEG+1
                  IF(NBEG.LT.NPP) GO TO 60
                  IF(IWR.NE.0) WRITE(6,613)
                  GO TO 999
                  ENDIF
              IF(NBEG.EQ.2) NBEG= 1
              IF((ITER.LE.1).AND.(IWR.NE.0)) THEN
                  IF(NBEG.GT.1) WRITE(6,609) JROT,EO,NBEG,YVB(NBEG)
                  IF(GB.LE.0.d0) WRITE(6,604) JROT,EO,NBEG,V(NBEG)/BFCT
                  ENDIF
!** Initialize outward wave function with a node:  WF(NBEG) = 0.
              SB= 0.d0
              SI= 1.d0
              GI= V(NBEG+1) - E*DRDY2(NBEG+1)
              Y1= SB*(1.d0 - HT*GB)
              Y2= SI*(1.d0 - HT*GI)
            ENDIF
!
          WF(NBEG)= SB
          WF(NBEG+1)= SI
          IF(INNER.GT.0) MSAVE= NPP
!** Actual outward integration loops start here
          DO  I= NBEG+2, MSAVE
              Y3= Y2 + Y2 - Y1 + GI*SI
              GI= V(I) - E*DRDY2(I)
              SI= Y3/(1.d0- HT*GI)
              WF(I)= SI
              IF(DABS(SI).GE.1.D+17) THEN
!** Renormalize to prevent overflow of  WF(I)  in classically forbidden
!  region where  V(I) .gt. E
                  SI= 1.d0/SI
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
              ITP1P= I
!** Exit from this loop at onset of classically allowed region
              IF(GI.LE.0.d0) GO TO 62
              ENDDO
          MS= MSAVE
          IF((INNER.EQ.0).AND.(GB.LE.0.d0)) GO TO 66
          IF(IWR.NE.0) WRITE(6,612) KVIN,JROT,EO,MSAVE
          GO TO 999
!** ITP1 is last point of AS-forbidden region & ITP1P 1'st point in allowed
   62     ITP1= ITP1P - 1
          MS= ITP1
          IF(INNER.GT.0) GO TO 66
          DO  I= ITP1P+1,MSAVE
              Y3= Y2 + Y2 - Y1 + GI*SI
              GI= V(I) - E*DRDY2(I)
              SI= Y3/(1.d0- HT*GI)
              WF(I)= SI
              IF(DABS(SI).GT.1.D+17) THEN
!** Renormalize to prevent overflow of  WF(I) , as needed.
                  SI= 1.d0/SI
                  DO  J= NBEG,I
                      WF(J)= WF(J)*SI
                      ENDDO
                  Y2= Y2*SI
                  Y3= Y3*SI
                  SI= 1.d0
                  ENDIF
              Y1= Y2
              Y2= Y3
              ENDDO
          MS= MSAVE
!** Finished outward integration.  Normalize w.r.t. WF(MSAVE)
   66     SI= 1.d0/SI
          YOUT= Y1*SI
          YM= Y2*SI
          RATIN= WF(NBEG+1)*SI
          DO  I= NBEG,MS
              WF(I)= WF(I)*SI
              ENDDO
          IF(INNER.GT.0) GO TO 10
!----- Finished numerical integration ... now correct trial energy
!** DF*H  is the integral of  (WF(I))**2 dR
   70     DF= 0.d0
          DO  J= NBEG,NEND
              DF= DF+ DRDY2(J)*WF(J)**2
              ENDDO
!** Add edge correction to DF assuming wave function dies off as simple
!  exponential past R(NEND);  matters only if WF(NEND) unusually large.
!!!       IF((E.LE.DSOC).AND.(WF(NEND).NE.0)) THEN
!!!
!!! huh ... how do I fix this for AS ??? - or is it no longer necessary ??
!!!
!!!           IF((KVIN.GE.-10).AND.(WF(NEND-1)/WF(NEND).GT.1.d0))
!!!  1              DF= DF+ WF(NEND)**2/(2.d0*DLOG(WF(NEND-1)/WF(NEND)))
!!!           ENDIF
!!!
!!!. note that by construction, at this point  WF(MSAVE)= 1.0
          F=  - YOUT - YIN + 2.d0*YM + GI
          DOLD= DE
          IF(DABS(F).LE.1.D+30) THEN
              DE= F/DF
            ELSE
              F= 9.9D+30
              DF= F
              DE= DABS(0.01D+0 *(DSOC-E))
            ENDIF
          IF(IWR.GT.2) THEN
              DEPRN = DE/BFCT
              XEND= YMINN+NEND*H
              WRITE(6,*) 'XEND=',XEND,YMINN ! Make them "referenced"
!** RATIN & RATOUT  are wave fx. amplitude at inner/outer ends of range
!  relative to its value at outermost extremum.
              WRITE(6,603) IT,EO,F,DF,DEPRN,MSAVE,RR,RATIN,RATOUT,      &
     &                                                       NBEG,ITP1
              ENDIF
!** Test trial eigenvalue for convergence
          IF(DABS(DE).LE.DABS(EPS)) GO TO 100
          E= E+DE
!** KV.ge.998  Option ... Search for highest bound level.  Adjust new
!  trial energy downward if it would have been above dissociation.
          IF((KVIN.GE.998).AND.(E.GT.VMX)) E= VMX- 2.d0*(VMX-E+DE)
          EO= E/BFCT
          IF((IT.GT.4).AND.(DABS(DE).GE.DABS(DOLD)).AND.                &
     &                                       ((DOLD*DE).LE.0.d0)) THEN
!** Adjust energy increment if having convergence difficulties.  Not
!  usually needed except for some quasibounds extremely near  VMAX .
              ICOR= ICOR+1
              DEP= DE/BFCT
              IF(IWR.NE.0) WRITE(6,617) IT,DEP
              DE= 0.5d0*DE
              E= E-DE
              EO= E/BFCT
              ENDIF
   90     CONTINUE
!** End of iterative loop which searches for eigenvalue ************
!-------------------------------------------------------------------*
!** Convergence fails, so return in error condition
      E= E-DE
      EO= E/BFCT
      DEPRN= DE/BFCT
      IF(IWR.NE.0) WRITE(6,620) KVIN,JROT,ITER,DEPRN
      GO TO 999
  100 IF(IWR.NE.0) THEN
          IF(IWR.GE.3) WRITE(6,619)
          IF((DABS(RATIN).GT.RATST).AND.(INNODE.GT.0)                   &
     &                 .AND.(YMIN.GT.0.d0)) WRITE(6,614) JROT,EO,RATIN
          IF((E.LT.DSOC).AND.(DABS(RATOUT).GT.RATST)) THEN
              WKBTST=0.5d0*DABS(V(NEND)-V(NEND-1))/DSQRT((V(NEND)-E)**3)
              IF(WKBTST.GT.1.d-3)WRITE(6,615)JROT,EO,RATOUT,RATST,WKBTST
              ENDIF
          ENDIF
      KKV = 0
!** Perform node count on converged solution
      PROD= WF(ITP1)*WF(ITP1-1)
      J1= ITP1+1
      J2= NEND-1
      DO  J= J1, J2
          PPROD= PROD
          PROD= WF(J)*WF(J-1)
          IF((PPROD.LE.0.d0).AND.(PROD.GT.0.d0)) KKV= KKV+1
          ENDDO
      KV = KKV

!     write(12,699) kv,jrot,EO,nend
! 699 format('   v=',i3,'    J='i3,'   E=',f10.3,'   NEND=',i6)

!** Normalize & find interval (NBEG,NEND) where WF(I) is non-negligible
      SN= 1.d0/DSQRT(H*DF)
      DO  I= NBEG,NEND
          WF(I)= WF(I)*SN
          ENDDO
      IF(ITP1.LE.1) GO TO 120
      J= ITP1P
      DO  I= 1,ITP1
          J= J-1
          IF(DABS(WF(J)).LT.RATST) GO TO 110
          ENDDO
  110 NBEG= J
      IF(NBEG.LE.1) GO TO 120
      J= J-1
      DO  I= 1,J
          WF(I)= 0.d0
          ENDDO
!** Move NEND inward to where wavefunction "non-negligible"
  120 J= NEND-1
      DO  I= NBEG,NEND
          IF(DABS(WF(J)).GT.RATST) GO TO 130
          J= J-1
          ENDDO
  130 NEND= J+1
      IF(NEND.LT.NPP) THEN
!** Zero out wavefunction array at distances past NEND
          DO  I= NEND+1,NPP
              WF(I)= 0.d0
              ENDDO
          ENDIF
      IF(LPRWF.LT.0) THEN
!** If desired, write every |LPRWF|-th point of wave function to a file
!  on channel-10, starting at mesh point # NBEG for radial distance
!  YVB(NBEG), with the NPR values separated by mesh step  JPSIQ*YH
          JPSIQ= -LPRWF
          NPR= 1+(NEND-NBEG)/JPSIQ
!** Write every JPSIQ-th point of the wave function for level  v=KV
!  J=JROT , beginning at mesh point NBEG & distance RSTT where
          WRITE(10,701) KV,JROT,EO,NPR,YVB(NBEG),YH*JPSIQ,NBEG,JPSIQ
          WRITE(10,702) (YVB(I),WF(I),I=NBEG,NEND,JPSIQ)
          ENDIF
      IF(IWR.EQ.1) WRITE(6,607) KV,JROT,EO
      IF(IWR.GE.2) THEN
          REND= ARV*((1.d0+YVB(NEND-1))/(1.d0-YVB(NEND-1)))**(1.d0/PRV)
          RATIN= RATIN*SDRDY(NBEG+1)/SDRDY(MSAVE)
          RATOUT= RATOUT*SDRDY(NEND-1)/SDRDY(MSAVE)
          WRITE(6,607) KV,JROT,EO,ITER,RR,RATIN,NBEG,REND,RATOUT,NEND-1
          ENDIF
!** For quasibound levels, calculate width in subroutine "WIDTH"
      IF((E.GT.DSOC).AND.(KVIN.GT.-10)) CALL WIDTHas(KV,JROT,E,EO,DSOC, &
     &  VBZ,WF,RVB,SDRDY,VMX,YMIN,H,BFCT,IWR,ITP1,ITP2,ITP3,INNER,NPP,  &
     &  GAMA)
      RETURN
!** ERROR condition if  E.gt.V(R)  at outer end of integration range.
  998 XPR= YMINN+MS*H
      VPR= V(MS)/BFCT
      IF(IWR.NE.0) WRITE(6,608) EO,MS,VPR,XPR,IT
!** Return in error mode
  999 KV= -1
!     CALL MMA_DEALLOCATE(RVB)
!     CALL MMA_DEALLOCATE(YVB)
!     CALL MMA_DEALLOCATE(DRDY2)
!     CALL MMA_DEALLOCATE(FAS)
!     CALL MMA_DEALLOCATE(SDRDY)
!     CALL MMA_DEALLOCATE(VBZ)
      RETURN
  601 FORMAT(/' Solve for  v=',I3,'   J=',I3,'   ETRIAL=',1PD15.7,      &
     &   '  INNER=',i2,'   WF(1st)  WF(NEND)' )
  602 FORMAT('ITER    ETRIAL',8X,'F(E)      DF(E)     D(E)',            &
     & 5X,'M    yp(M)   /WF(M)    /WF(M)  NBEG  ITP1'/                  &
     &  1X,96('-'))
  603 FORMAT(I3,1PD15.7,3D10.2,0P,I6,F7.3,1P2D9.1,0P,I5,I6)
  604 FORMAT('   NOTE:  for  J=',I3,'   EO=',F12.4,' .ge. V(',i3,')=',  &
     &  F12.4)
  607 FORMAT('E(v=',I3,',J=',I3,')=',F15.8,I3,' iterations',            &
     & '  yp(M)=',F6.3,'  WF(NBEG)/WF(M)=',1PD8.1,0P,'   NBEG=',i5/40x, &
     & 'R(NEND)=',f9.2,'   WF(NEND)/WF(M)=',1PD8.1,0P,'   NEND=',i5)
  608 FORMAT(' *** SCHRQ Error:  E=',F9.2,' > V(',I5,')=',F9.2,         &
     &  '  at  Rmax=',F6.2,'  for  IT=',I2)
 6609 FORMAT(' *** For  J=',I3,'   E=',1PD15.7,"  integration can't",   &
     & ' start till inside'/21x,'mesh point',I6,' (yp=',0pf8.4,         &
     &  '),  so YMAX larger than needed')
  609 FORMAT(' *** For  J=',I3,'   E=',1PD15.7,"  integration can't",   &
     & ' start till past'/23x,'mesh point',I5,' (yp=',0pf6.2,           &
     &  '),  so YMIN smaller than needed')
  610 FORMAT(/' Seek highest bound level:   ETRIAL =',1PD9.2,17x,       &
     &  'WF(1st)  WF(NEND)')
  611 FORMAT(' *** SCHRQ inward search at   J=',i3,'   E=',f11.2,       &
     &  ' finds no classical region')
  612 FORMAT(/' *** ERROR *** for   v =',I3,'   J =',I3,'   E =',       &
     &  F12.4,'  Innermost turning point not found by   M = MSAVE =',I5)
  613 FORMAT(/' *** ERROR in potential array ... V(I) everywhere',      &
     & ' too big to integrate with given  increment')
  614 FORMAT(' *** CAUTION *** For  J=',I3,'  E=',G15.8/16x,            &
     & 'WF(first)/WF(Max)=',D9.2,'  suggests  YMIN  may be too large')
  615 FORMAT(' ** CAUTION ** For  J=',I3,'  E=',1PD13.6,                &
     & '  WF(NEND)/WF(Max)=',D8.1,' >',D8.1/4X,'& initialization ',     &
     & 'quality test ',1PD8.1,' > 1.D-3   so RMAX may be too small')
  616 FORMAT(' ** WARNING *** For  v=',I2,', J=',I3,' at  E=',G14.7,    &
     &  ':  inward propagation finds no turning point ... v=-2 means:   &
     & Trial energy is too low (!), or potential is too weak' )
  617 FORMAT(' *** SCHRQ has a convergence problem, so for  IT=',I2,    &
     & '  cut  DE=',1PD10.2,'  in HALF' )
  618 FORMAT(' *** For  J=',I3,'  E=',F9.2,'  JWKB start gives  SB/SI=',&
     &  1PD10.3,'  so use a node.')
  619 FORMAT(1X,96('-'))
  620 FORMAT(' *** CAUTION for  v=',I3,'  J=',I3,"  SCHRQ doesn't conver&
     &ge by  ITER=",I2,'  DE=',1PD9.2)
  701 FORMAT(/2x,'Level  v=',I3,'   J=',I3,'   E=',F12.4,' ,  wave funct&
     &ion at',I6,' points.'/7x,'R(1-st)=',F12.8,'   mesh=',F12.8,       &
     &  '   NBEG=',I4,'   |LPRWF|=',I3)
  702 FORMAT((4(f10.6,f11.7)))
      END
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

!***********************************************************************
      SUBROUTINE QBOUNDas(KV,JROT,E,EO,VMX,DSOC,VBZ,SDRDY,RVB,YMIN,YH,  &
     &  GB,GI,SB,SI,NPP,ITP2,ITP3,IWR,IQTST,BFCT,IT)
!***********************************************************************
!** Subroutine to initialize quasibound level wave function as Airy
!  function at third turning point (if possible). For the theory see
!  J.Chem.Phys. 54, 5114 (1971),  J.Chem.Phys. 69, 3622-31 (1978)
!----------------------------------------------------------------------
!** IQTST  is error flag. *** If (IQTST.lt.0) initialization fails
!  so eigenvalue calculation aborts *** (IQTST.gt.0) for successful
!  Airy function initialization. *** (IQTST=0) if Airy function
!  initialization prevented because 3-rd turning point beyond
!  range, so that WKB initialization is used.
!----------------------------------------------------------------------
      INTEGER I,II,IQTST,IT,ITP2,ITP3,IWR,J,JROT,KV,NPP
      REAL*8  A1,A2,A13,A23,BFCT,C1A,C2A,DSOC,E,EO,FBA,FIA,FJ,GB,GI,    &
     &  GBA,GIA,YH,RH,YMIN,YMINN,SB,SI,SL,VBZ(NPP),SDRDY(NPP),RVB(NPP), &
     &  VMX,VMXPR
      DATA C1A/0.355028053887817D0/,C2A/0.258819403792807D0/
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IQTST=1
      YMINN= YMIN- YH
      WRITE(6,*) 'YMINN=',YMINN ! Make sure YMINN is referenced.
!** Start by searching for third turning point.
      J=NPP-1
      IF(VBZ(J).GT.E) GO TO 30
      DO  I=NPP-2,1,-1
          J=J-1
          IF(VBZ(J).GT.E) GO TO 10
          ENDDO
      IQTST= -9
      WRITE(6,602) JROT,EO
      RETURN
!** ITP3 is the first mesh point outside classically forbidden region
   10 ITP3= J+1
!** Check that there is a classically allowed region inside this point
!  and determine height of barrier maximum.
      II=J
      VMX=DSOC
      DO  I=2,J
          II=II-1
          IF(VBZ(II).LE.E) GO TO 20
          IF(VBZ(II).GT.VMX) VMX= VBZ(II)
          ENDDO
!** Energy too high (or too low): find only one turning point.
      VMXPR= VMX/BFCT
      IF(IWR.NE.0) WRITE(6,604) JROT,EO,VMXPR/BFCT,RVB(J)
      IQTST=-1
      RETURN
!** ITP2 is first mesh point inside forbidden region on left of barrier
   20 ITP2= II+1
!** Now ... continue to set up r3(E) boundary condition ...
      RH= RVB(ITP3)- RVB(ITP3-1)
      GB= (VBZ(ITP3) - E)*(RH/YH)**2
      GI= (VBZ(ITP3-1) - E)*(RH/YH)**2
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING
!     WRITE(6,*) 'ITP3=',ITP3
!     WRITE(6,*) 'VBZ(ITP3-1)=',VBZ(ITP3-1)
!     WRITE(6,*) 'E=',E
!     WRITE(6,*) 'RH=',RH
!     WRITE(6,*) 'YH=',YH
      FJ= GI/(GI-GB)
      WRITE(6,*) FJ ! make sure it's "referenced"
!** Treat quasibound levels as bound using outer boundary condition
!  of Airy function at third turning point ... as discussed by
!  R.J.Le Roy and R.B.Bernstein  in  J.Chem.Phys. 54,5114(1971).
!  Uses series expansions of Abramowitz & Stegun Eq.(10.4.3)
      SL= (GI-GB)**(1.d0/3.d0)/RH
      A1= GI/(SL*RH)**2
      A2= GB/(SL*RH)**2
      A13= A1*A1*A1
      A23= A2*A2*A2
      FIA= 1.d0+ A13*(A13*(A13+72.D0)+2160.D0)/12960.D0
      GIA= A1+A1*A13*(A13*(A13+90.D0)+3780.D0)/45360.D0
      FBA= 1.d0+ A23*(A23*(A23+72.D0)+2160.D0)/12960.D0
      GBA= A2+A2*A23*(A23*(A23+90.D0)+3780.D0)/45360.D0
!** Airy function  Bi(X)  at points straddling 3-rd turning point
      SI= (C1A*FIA+C2A*GIA)/SDRDY(ITP3-1)
      SB= (C1A*FBA+C2A*GBA)/SDRDY(ITP3)
      GI= VBZ(ITP3-1) - E
      GB= VBZ(ITP3) - E
      IF(SB.GE.SI) THEN
!** In case of big error - switch to node at ITP3
          SB= 0.d0
          SI= 1.d0
          IF(IWR.NE.0) WRITE(6,606) KV,JROT,EO,IT
          ENDIF
      RETURN
!
!** If 3-rd turning point beyond range start with WKB wave function
!  at end of range.
   30 IF(IWR.NE.0) WRITE(6,608) JROT,EO
      ITP3= NPP-1
      IQTST= 0
      VMX= VBZ(ITP3)
      II= ITP3
!... and determine barrier maximum ....
      DO  I= 2,ITP3
          II= II-1
          VMXPR= VBZ(II)
          IF(VMXPR.LT.VMX) GO TO 40
          VMX= VMXPR
          ENDDO
      IF(IWR.NE.0) WRITE(6,610)
      IQTST= -9
   40 RETURN
  602 FORMAT(' *** QBOUND fails for   E(J=',i3,')=',f9.3,'  Find no turn&
     &ing point')
  604 FORMAT(' For J=',I3,'  ETRY=',F11.4,' > VMAX=',F11.4,             &
     &  '  find onee turn point:  R=',F6.2)
  606 FORMAT(' *** CAUTION ***  v=',I3,'   J=',I3,'   E=',1PD13.6,      &
     & '   IT=',I2/5x,'Airy initialization unstable so place node just p&
     &ast  R(3-rd)' )
  608 FORMAT(' *** For  J=',I3,'  E=',F9.2,                             &
     &  '  R(3-rd) > RMAX  & E < V(N)  so try WKB B.C. @ RMAX')
  610 FORMAT(" **** QBOUND doesn't work ... no classically allowed regio&
     &n accessible at this energy.")
      END
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
