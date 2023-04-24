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
      SUBROUTINE CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,V,WF0,RM2,RCNST)
      USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Subroutine solving the linear inhomogeneous differential equations
!  formulated by J.M. Hutson [J.Phys.B14, 851 (1982)] for treating
!  centrifugal distortion as a perturbation, to determine centrifugal
!  distortion constants of a diatomic molecule.  Uses the algorithm of
!  J. Tellinghuisen [J.Mol.Spectrosc. 122, 455 (1987)].  The current
!  version calculates Bv, Dv, Hv, Lv, Mv, Nv and Ov and writes them out,
!  but does not return values to the calling program.
!
!** On entry:   EO    is the eigenvalue (in units [cm-1])
!               NBEG & NEND  the mesh point range over which the input
!               NDIMR  is dimension of arrays  V(i), WF0(i), ... etc.
! wavefunction  WF0  (in units 1/sqrt(Ang))  has non-negligible values
!               BvWn  is the numerical factor (hbar^2/2mu) [cm-1 Ang^2]
!               YH    is the integration stepsize
!               WARN  is an integer flag: > 0 print internal warnings,
!               V(i)  is the effective potential (including centrifugal
!                     term if calculation performed at  J > 0) in
!                     'internal' units, including the factor  YH**2/BvWN
!               RM2(i) is the array  (r')^2/(distance**2)
!** On exit:    RCNST(i)  is the set of 7 rotational constants: Bv, -Dv,
!                       Hv, Lv, Mv, Nv & Ov
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Dimension:  potential arrays  and  vib. level arrays.
!!!
      INTEGER NDIMR
!     PARAMETER (NDIMR= 200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive.
!     PARAMETER (NDIMR= 131074)
!     REAL*8 PRV,ARV,RVB(NDIMR),YVB(NDIMR),DRDY2(NDIMR),FAS(NDIMR),
!    1                                         SDRDY(NDIMR),VBZ(NDIMR)
      REAL*8 PRV,ARV
      REAL*8, ALLOCATABLE :: RVB(:),YVB(:),DRDY2(:),FAS(:),             &
     &                                         SDRDY(:),VBZ(:)
      COMMON /BLKAS/PRV,ARV!,RVB,YVB,DRDY2,SDRDY,FAS,VBZ
!!!
      INTEGER I,M,IPASS,M1,M2,NBEG,NEND,WARN
!      REAL*8 V(NEND),WF0(NEND),RM2(NEND),P(NDIMR),WF1(NDIMR),
!     1                                            WF2(NDIMR),RCNST(7)
      REAL*8 V(NEND),WF0(NEND),RM2(NEND),RCNST(7)
      REAL*8, ALLOCATABLE :: P(:),WF1(:),WF2(:)
      REAL*8 BvWN,DV,DVV,HVV,HV2,LVV,LV2,MVV,MV2,NVV,OVV,EO,E,YH,YH2,   &
     &  ZTW,AR,R2IN,G2,G3,P0,P1,P2,P3,PI,PIF,PRS,PRT,V1,V2,V3,Y1,Y2,Y3, &
     &  TSTHv,TSTLv,TSTMv,AMB,AMB1,AMB2,                                &
     &  OV,OV01,OV02,OV03,OV11,OV12,OV13,OV22,OV23,OV33,                &
     &  PER01,PER02,PER03,PER11,PER12,PER13,PER22,PER23,PER33,R2XX
!
      NDIMR = 131074
      CALL MMA_ALLOCATE(RVB,NDIMR,LABEL='RVB')
      CALL MMA_ALLOCATE(YVB,NDIMR,LABEL='YVB')
      CALL MMA_ALLOCATE(DRDY2,NDIMR,LABEL='DRDY2')
      CALL MMA_ALLOCATE(FAS,NDIMR,LABEL='FAS')
      CALL MMA_ALLOCATE(SDRDY,NDIMR,LABEL='SDRDY')
      CALL MMA_ALLOCATE(VBZ,NDIMR,LABEL='VBZ')
!
      CALL MMA_ALLOCATE(P,NDIMR,LABEL='P')
      CALL MMA_ALLOCATE(WF1,NDIMR,LABEL='WF1')
      CALL MMA_ALLOCATE(WF2,NDMIR,LABEL='WF2')
      P0=0
      MV2=0
      LV2=0
      G3=0
      IF(NEND.GT.NDIMR) THEN
          WRITE(6,602) NEND,NDIMR
          RETURN
          ENDIF
      ZTW= 1.D0/12.d0
      YH2 = YH*YH
      DV = YH2*ZTW
      E= EO*YH2/BvWN
      IPASS = 1
      OV01 = 0.D0
      OV02 = 0.D0
      OV03 = 0.D0
      OV11 = 0.D0
      OV22 = 0.D0
      OV12 = 0.D0
      OV33 = 0.D0
      OV23 = 0.D0
      OV13 = 0.D0
      PER01 = 0.D0
      PER02 = 0.D0
      PER03 = 0.D0
      PER11 = 0.D0
      PER12 = 0.D0
      PER13 = 0.D0
      PER22 = 0.D0
      PER23 = 0.D0
      PER33 = 0.D0
!** First, calculate the expectation value of  1/r**2  and hence Bv
      R2IN= 0.5D0*(RM2(NBEG)*WF0(NBEG)**2 + RM2(NEND)*WF0(NEND)**2)
      DO   I= NBEG+1, NEND-1
         R2IN= R2IN+ RM2(I)*WF0(I)**2
         ENDDO
      R2IN = R2IN*YH
      RCNST(1)= R2IN*BvWN
!
!** On First pass  IPASS=1  and calculate first-order wavefx., Dv & Hv
!  On second pass  IPASS=2  and calculate second-order wavefx., Lv & Mv
!  On third pass   IPASS=3  and calculate third-order wavefx., Nv & Ov
!
   10 P1= 0.D0
      P2= 0.D0
!
!     P1= WF0(NEND)
!     P2= WF0(NEND-1)
!
      P(NEND) = P1
      P(NEND-1) = P2
      V1 = V(NEND) - E*DRDY2(NEND)
      V2 = V(NEND-1) - E*DRDY2(NEND-1)
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1)                                       &
     &                   - DV*(RM2(NEND) - R2IN*DRDY2(NEND))*WF0(NEND)
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF0(NEND-1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF0(NEND)                   &
     &                     - (RM2(NEND) - R2IN*DRDY2(NEND))*WF1(NEND))
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF1(NEND-1)           &
     &                                 - DVV*WF0(NEND-1)*DRDY2(NEND-1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF1(NEND) - HVV*WF0(NEND)   &
     &                     - (RM2(NEND) - R2IN*DRDY2(NEND))*WF2(NEND))
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF2(NEND-1)           &
     &             - (DVV*WF1(NEND-1) + HVV*WF0(NEND-1))*DRDY2(NEND-1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      M= NEND-1
!** Now - integrate inward from outer end of range
      DO  I = NBEG+2,NEND
          M = M-1
          Y3 = Y2 + Y2 - Y1 + YH2*G2 + V2*P2
          R2XX= R2IN*DRDY2(M)
          IF(IPASS.EQ.1) G3 = (RM2(M)- R2XX)*WF0(M)
          IF(IPASS.EQ.2) G3 = (RM2(M)-R2XX)*WF1(M) - DVV*WF0(M)*DRDY2(M)
          IF(IPASS.EQ.3) G3 = (RM2(M)- R2XX)*WF2(M)                     &
     &                            - (DVV*WF1(M) + HVV*WF0(M))*DRDY2(M)
          V3 = V(M) - E*DRDY2(M)
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          IF(V3.LT.0.D0)  GO TO 32
          P(M) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          ENDDO
      GO TO 90
!** Escaped loop at outer turning point:  initialize outward integration
   32 PRS = P3
      PRT = P(M+1)
      P1 = 0.D0
      P2 = 0.D0
!
!     P1 = WF0(NBEG)
!     P2 = WF0(NBEG+1)
!
      P(NBEG) = P1
      P(NBEG+1) = P2
      V1 = V(NBEG) - E*DRDY2(NBEG)
      V2 = V(NBEG+1) - E*DRDY2(NBEG+1)
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1)                                       &
     &                   - DV*(RM2(NBEG) - R2IN*DRDY2(NBEG))*WF0(NBEG)
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF0(NEND)*DRDY2(NBEG)       &
     &                     - (RM2(NBEG) - R2IN*DRDY2(NBEG))*WF1(NBEG))
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF1(NBEG+1)           &
     &                                 - DVV*WF0(NBEG+1)*DRDY2(NBEG+1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*((DVV*WF1(NEND) + HVV*WF0(NEND)) &
     &       *DRDY2(NBEG)  - (RM2(NBEG) - R2IN*DRDY2(NBEG))*WF2(NBEG))
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF2(NBEG+1)           &
     &             - (DVV*WF1(NBEG+1) + HVV*WF0(NBEG+1))*DRDY2(NBEG+1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      AR = 0.D0
      M1 = M+1
!** Now ... integrate outward from inner end of range
      DO  I = NBEG+2,M1
          Y3 = Y2 + Y2 - Y1 + YH2*G2 + V2*P2
          P0 = WF0(I)
          R2XX= R2IN*DRDY2(I)
          IF(IPASS.EQ.1) G3 = (RM2(I)-R2XX)*P0
          IF(IPASS.EQ.2) G3 = (RM2(I)-R2XX)*WF1(I) - DVV*P0*DRDY2(I)
          IF(IPASS.EQ.3) G3 = (RM2(I)-R2XX)*WF2(I)                      &
     &                                - (DVV*WF1(I) + HVV*P0)*DRDY2(I)
          V3 = V(I) - E*DRDY2(I)
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          P(I) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          AR = AR + P0*P3*DRDY2(I)
          ENDDO
!** Average for 2 adjacent mesh points to get Joel's "(a-b)"
      AMB2 = (P3-PRT)/P0
      AMB1 = (P(M)-PRS)/WF0(M)
      AMB = (AMB1+AMB2)*0.5D0
      M2 = M+2
!** Find the rest of the overlap with zero-th order solution ...
      DO  I = M2,NEND
          P0 = WF0(I)
          PI = P(I) + AMB*P0
          P(I) = PI
          AR = AR + PI*P0*DRDY2(I)
          ENDDO
      OV = AR*YH
      DO  I = NBEG,NEND
          P0 = WF0(I)
! ... and project out contribution of zero'th-order part of solution
          PI = P(I) - OV*P0
          PIF = PI*RM2(I)
          IF(IPASS.EQ.1) THEN
!** Now - on first pass accumulate integrals for Dv and Hv
              WF1(I) = PI
              OV01 = OV01 + PI*P0 * drdy2(i)
              OV11 = OV11 + PI*PI * drdy2(i)
              PER01 = PER01 + PIF*P0
              PER11 = PER11 + PI*PIF
            ELSEIF(IPASS.EQ.2) THEN
! ... and on next pass, accumulate integrals for Lv and Mv
              WF2(I) = PI
              P1 = WF1(I)
              OV02 = OV02 + PI*P0 * drdy2(i)
              OV12 = OV12 + PI*P1 * drdy2(i)
              OV22 = OV22 + PI*PI * drdy2(i)
              PER02 = PER02 + PIF*P0
              PER12 = PER12 + PIF*P1
              PER22 = PER22 + PI*PIF
            ELSEIF(IPASS.EQ.3) THEN
! ... and on next pass, accumulate integrals for Nv and Ov
              P1 = WF1(I)
              P2 = WF2(I)
              OV03 = OV03 + PI*P0 * drdy2(i)
              OV13 = OV13 + PI*P1 * drdy2(i)
              OV23 = OV23 + PI*P2 * drdy2(i)
              OV33 = OV33 + PI*PI * drdy2(i)
              PER03 = PER03 + PIF*P0
              PER13 = PER13 + PIF*P1
              PER23 = PER23 + PIF*P2
              PER33 = PER33 + PIF*PI
            ENDIF
          ENDDO
      IF(IPASS.EQ.1) THEN
          DVV = YH*PER01
          HVV = YH*(PER11 - R2IN*OV11)
          IPASS = 2
          RCNST(2) = DVV*BvWN
          RCNST(3) = HVV*BvWn
          GO TO 10
        ELSEIF(IPASS.EQ.2) THEN
          HV2 = YH*PER02*BvWN
          LVV = YH*(PER12 - R2IN*OV12 - DVV*OV11)
          MVV = YH*(PER22 - R2IN*OV22 - 2.D0*DVV*OV12 - HVV*OV11)
          IPASS = 3
          RCNST(4) = LVV*BvWN
          RCNST(5) = MVV*BvWN
          GO TO 10
        ELSEIF(IPASS.EQ.3) THEN
          LV2 = YH*PER03*BvWN
          MV2 = YH*(PER13 - R2IN*OV13 - DVV*OV12 - HVV*OV11)*BvWN
          NVV = YH*(PER23 - R2IN*OV23 - DVV*(OV13 + OV22)               &
     &                                     - 2.D0*HVV*OV12 - LVV*OV11)
          OVV = YH*(PER33 - R2IN*OV33 - 2.D0*DVV*OV23                   &
     &             - HVV*(2.D0*OV13+ OV22) - 2.D0*LVV*OV12 - MVV*OV11)
          RCNST(6) = NVV*BvWN
          RCNST(7) = OVV*BvWN
        ENDIF
      IF(WARN.GT.0) THEN
          IF(DMAX1(DABS(OV01),DABS(OV02),DABS(OV01)).GT.1.D-9)          &
     &                                     WRITE(6,604) OV01,OV02,OV03
          TSTHV= dabs(RCNST(3)/HV2-1.D0)
          TSTLV= dabs(RCNST(4)/LV2-1.D0)
          TSTMV= dabs(RCNST(5)/MV2-1.D0)
          IF(DMAX1(TSTHV,TSTLV,TSTMV).GT.1.d-5)                         &
     &                                  WRITE(6,603) TSTHV,TSTLV,TSTMV
          ENDIF
      CALL MMA_DEALLOCATE(RVB)
      CALL MMA_DEALLOCATE(YVB)
      CALL MMA_DEALLOCATE(DRDY2)
      CALL MMA_DEALLOCATE(FAS)
      CALL MMA_DEALLOCATE(SDRDY)
      CALL MMA_DEALLOCATE(VBZ)
!
      CALL MMA_DEALLOCATE(P)
      CALL MMA_DEALLOCATE(WF1)
      CALL MMA_DEALLOCATE(WF2)
      RETURN
   90 WRITE(6,601) EO
      RETURN
  601 FORMAT(' *** ERROR in CDJOEL *** for input energy  E =',f12.4,    &
     &   '  never reach outer turning point')
  602 FORMAT(/' *** Dimensioning PROBLEM in CDJOEL ***   NEND=',i6,     &
     &  ' > NDIMR=',i6)
  603 FORMAT(' ** CAUTION ** Comparison tests for Hv, Lv & Mv give:',   &
     & 3(1Pd9.1))
  604 FORMAT(' ** CAUTION ** CDJOEL orthogonality tests OV01,OV02 & OV03&
     &:',3(1Pd9.1))
      END
