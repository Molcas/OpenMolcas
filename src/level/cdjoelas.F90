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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,V,WF0,RM2,RCNST)
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

use LEVEL_COMMON, only: DRDY2, NDIMR
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Twelve, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NBEG, NEND, WARN
real(kind=wp), intent(in) :: EO, BvWN, YH, V(NEND), WF0(NEND), RM2(NEND)
real(kind=wp), intent(out) :: RCNST(7)
integer(kind=iwp) :: I, IPASS, M, M1, M2
real(kind=wp) :: AMB, AMB1, AMB2, AR, DV, DVV, E, G2, G3, HV2, HVV, LV2, LVV, MV2, MVV, NVV, OV, OV01, OV02, OV03, OV11, OV12, &
                 OV13, OV22, OV23, OV33, OVV, P0, P1, P2, P3, PER01, PER02, PER03, PER11, PER12, PER13, PER22, PER23, PER33, PI, &
                 PIF, PRS, PRT, R2IN, R2XX, TSTHv, TSTLv, TSTMv, V1, V2, V3, Y1, Y2, Y3, YH2, ZTW
logical(kind=iwp) :: Break
real(kind=wp), allocatable :: P(:), WF1(:), WF2(:)

P0 = 0
MV2 = 0
LV2 = 0
G3 = 0
if (NEND > NDIMR) then
  write(u6,602) NEND,NDIMR
  return
end if
call mma_allocate(P,NEND,LABEL='P')
call mma_allocate(WF1,NEND,LABEL='WF1')
call mma_allocate(WF2,NEND,LABEL='WF2')
ZTW = One/Twelve
YH2 = YH*YH
DV = YH2*ZTW
E = EO*YH2/BvWN
OV01 = Zero
OV02 = Zero
OV03 = Zero
OV11 = Zero
OV22 = Zero
OV12 = Zero
OV33 = Zero
OV23 = Zero
OV13 = Zero
PER01 = Zero
PER02 = Zero
PER03 = Zero
PER11 = Zero
PER12 = Zero
PER13 = Zero
PER22 = Zero
PER23 = Zero
PER33 = Zero
DVV = Zero
HVV = Zero
LVV = Zero
MVV = Zero
HV2 = Zero
!** First, calculate the expectation value of  1/r**2  and hence Bv
R2IN = Half*(RM2(NBEG)*WF0(NBEG)**2+RM2(NEND)*WF0(NEND)**2)
R2IN = R2IN+sum(RM2(NBEG+1:NEND-1)*WF0(NBEG+1:NEND-1)**2)*YH
RCNST(1) = R2IN*BvWN

! On First pass  IPASS=1  and calculate first-order wavefx., Dv & Hv
! On second pass IPASS=2  and calculate second-order wavefx., Lv & Mv
! On third pass  IPASS=3  and calculate third-order wavefx., Nv & Ov

do IPASS=1,3
  P1 = Zero
  P2 = Zero

  !P1 = WF0(NEND)
  !P2 = WF0(NEND-1)

  P(NEND) = P1
  P(NEND-1) = P2
  V1 = V(NEND)-E*DRDY2(NEND)
  V2 = V(NEND-1)-E*DRDY2(NEND-1)
  select case (IPASS)
    case (1)
      Y1 = P1*(One-ZTW*V1)-DV*(RM2(NEND)-R2IN*DRDY2(NEND))*WF0(NEND)
      G2 = (RM2(NEND-1)-R2IN*DRDY2(NEND-1))*WF0(NEND-1)
    case (2)
      Y1 = P1*(One-ZTW*V1)+DV*(DVV*WF0(NEND)-(RM2(NEND)-R2IN*DRDY2(NEND))*WF1(NEND))
      G2 = (RM2(NEND-1)-R2IN*DRDY2(NEND-1))*WF1(NEND-1)-DVV*WF0(NEND-1)*DRDY2(NEND-1)
    case (3)
      Y1 = P1*(One-ZTW*V1)+DV*(DVV*WF1(NEND)-HVV*WF0(NEND)-(RM2(NEND)-R2IN*DRDY2(NEND))*WF2(NEND))
      G2 = (RM2(NEND-1)-R2IN*DRDY2(NEND-1))*WF2(NEND-1)-(DVV*WF1(NEND-1)+HVV*WF0(NEND-1))*DRDY2(NEND-1)
  end select
  Y2 = P2*(One-ZTW*V2)-DV*G2
  M = NEND-1
  ! Now - integrate inward from outer end of range
  Break = .false.
  do I=NBEG+2,NEND
    M = M-1
    Y3 = Y2+Y2-Y1+YH2*G2+V2*P2
    R2XX = R2IN*DRDY2(M)
    select case (IPASS)
      case (1)
        G3 = (RM2(M)-R2XX)*WF0(M)
      case (2)
        G3 = (RM2(M)-R2XX)*WF1(M)-DVV*WF0(M)*DRDY2(M)
      case (3)
        G3 = (RM2(M)-R2XX)*WF2(M)-(DVV*WF1(M)+HVV*WF0(M))*DRDY2(M)
    end select
    V3 = V(M)-E*DRDY2(M)
    P3 = (Y3+DV*G3)/(One-ZTW*V3)
    if (V3 < Zero) then
      Break = .true.
      exit
    end if
    P(M) = P3
    Y1 = Y2
    Y2 = Y3
    V2 = V3
    P2 = P3
    G2 = G3
  end do
  ! Escaped loop at outer turning point:  initialize outward integration
  if (Break) then
    PRS = P3
    PRT = P(M+1)
    P1 = Zero
    P2 = Zero

    !P1 = WF0(NBEG)
    !P2 = WF0(NBEG+1)

    P(NBEG) = P1
    P(NBEG+1) = P2
    V1 = V(NBEG)-E*DRDY2(NBEG)
    V2 = V(NBEG+1)-E*DRDY2(NBEG+1)
    select case (IPASS)
      case (1)
        Y1 = P1*(One-ZTW*V1)-DV*(RM2(NBEG)-R2IN*DRDY2(NBEG))*WF0(NBEG)
        G2 = (RM2(NBEG+1)-R2IN*DRDY2(NBEG+1))*WF0(NBEG+1)
      case (2)
        Y1 = P1*(One-ZTW*V1)+DV*(DVV*WF0(NEND)*DRDY2(NBEG)-(RM2(NBEG)-R2IN*DRDY2(NBEG))*WF1(NBEG))
        G2 = (RM2(NBEG+1)-R2IN*DRDY2(NBEG+1))*WF1(NBEG+1)-DVV*WF0(NBEG+1)*DRDY2(NBEG+1)
      case (3)
        Y1 = P1*(One-ZTW*V1)+DV*((DVV*WF1(NEND)+HVV*WF0(NEND))*DRDY2(NBEG)-(RM2(NBEG)-R2IN*DRDY2(NBEG))*WF2(NBEG))
        G2 = (RM2(NBEG+1)-R2IN*DRDY2(NBEG+1))*WF2(NBEG+1)-(DVV*WF1(NBEG+1)+HVV*WF0(NBEG+1))*DRDY2(NBEG+1)
    end select
    Y2 = P2*(One-ZTW*V2)-DV*G2
    AR = Zero
    M1 = M+1
    ! Now ... integrate outward from inner end of range
    do I=NBEG+2,M1
      Y3 = Y2+Y2-Y1+YH2*G2+V2*P2
      P0 = WF0(I)
      R2XX = R2IN*DRDY2(I)
      select case (IPASS)
        case (1)
          G3 = (RM2(I)-R2XX)*P0
        case (2)
          G3 = (RM2(I)-R2XX)*WF1(I)-DVV*P0*DRDY2(I)
        case (3)
          G3 = (RM2(I)-R2XX)*WF2(I)-(DVV*WF1(I)+HVV*P0)*DRDY2(I)
      end select
      V3 = V(I)-E*DRDY2(I)
      P3 = (Y3+DV*G3)/(One-ZTW*V3)
      P(I) = P3
      Y1 = Y2
      Y2 = Y3
      V2 = V3
      P2 = P3
      G2 = G3
      AR = AR+P0*P3*DRDY2(I)
    end do
    ! Average for 2 adjacent mesh points to get Joel's "(a-b)"
    AMB2 = (P3-PRT)/P0
    AMB1 = (P(M)-PRS)/WF0(M)
    AMB = (AMB1+AMB2)*Half
    M2 = M+2
    ! Find the rest of the overlap with zero-th order solution ...
    do I=M2,NEND
      P0 = WF0(I)
      PI = P(I)+AMB*P0
      P(I) = PI
      AR = AR+PI*P0*DRDY2(I)
    end do
    OV = AR*YH
    do I=NBEG,NEND
      P0 = WF0(I)
      ! ... and project out contribution of zero'th-order part of solution
      PI = P(I)-OV*P0
      PIF = PI*RM2(I)
      select case (IPASS)
        case (1)
          ! Now - on first pass accumulate integrals for Dv and Hv
          WF1(I) = PI
          OV01 = OV01+PI*P0*DRDY2(i)
          OV11 = OV11+PI*PI*DRDY2(i)
          PER01 = PER01+PIF*P0
          PER11 = PER11+PI*PIF
        case (2)
          ! ... and on next pass, accumulate integrals for Lv and Mv
          WF2(I) = PI
          P1 = WF1(I)
          OV02 = OV02+PI*P0*DRDY2(i)
          OV12 = OV12+PI*P1*DRDY2(i)
          OV22 = OV22+PI*PI*DRDY2(i)
          PER02 = PER02+PIF*P0
          PER12 = PER12+PIF*P1
          PER22 = PER22+PI*PIF
        case (3)
          ! ... and on next pass, accumulate integrals for Nv and Ov
          P1 = WF1(I)
          P2 = WF2(I)
          OV03 = OV03+PI*P0*DRDY2(i)
          OV13 = OV13+PI*P1*DRDY2(i)
          OV23 = OV23+PI*P2*DRDY2(i)
          OV33 = OV33+PI*PI*DRDY2(i)
          PER03 = PER03+PIF*P0
          PER13 = PER13+PIF*P1
          PER23 = PER23+PIF*P2
          PER33 = PER33+PIF*PI
      end select
    end do
    select case (IPASS)
      case (1)
        DVV = YH*PER01
        HVV = YH*(PER11-R2IN*OV11)
        RCNST(2) = DVV*BvWN
        RCNST(3) = HVV*BvWn
      case (2)
        HV2 = YH*PER02*BvWN
        LVV = YH*(PER12-R2IN*OV12-DVV*OV11)
        MVV = YH*(PER22-R2IN*OV22-Two*DVV*OV12-HVV*OV11)
        RCNST(4) = LVV*BvWN
        RCNST(5) = MVV*BvWN
      case (3)
        LV2 = YH*PER03*BvWN
        MV2 = YH*(PER13-R2IN*OV13-DVV*OV12-HVV*OV11)*BvWN
        NVV = YH*(PER23-R2IN*OV23-DVV*(OV13+OV22)-Two*HVV*OV12-LVV*OV11)
        OVV = YH*(PER33-R2IN*OV33-Two*DVV*OV23-HVV*(Two*OV13+OV22)-Two*LVV*OV12-MVV*OV11)
        RCNST(6) = NVV*BvWN
        RCNST(7) = OVV*BvWN
        if (WARN > 0) then
          if (max(abs(OV01),abs(OV02),abs(OV01)) > 1.0e-9_wp) write(u6,604) OV01,OV02,OV03
          TSTHV = abs(RCNST(3)/HV2-One)
          TSTLV = abs(RCNST(4)/LV2-One)
          TSTMV = abs(RCNST(5)/MV2-One)
          if (max(TSTHV,TSTLV,TSTMV) > 1.0e-5_wp) write(u6,603) TSTHV,TSTLV,TSTMV
        end if
    end select
  else
    write(u6,601) EO
    exit
  end if
end do
call mma_deallocate(P)
call mma_deallocate(WF1)
call mma_deallocate(WF2)

return

601 format(' *** ERROR in CDJOEL *** for input energy  E =',f12.4,'  never reach outer turning point')
602 format(/' *** Dimensioning PROBLEM in CDJOEL ***   NEND=',i6,' > NDIMR=',i6)
603 format(' ** CAUTION ** Comparison tests for Hv, Lv & Mv give:',3(1Pd9.1))
604 format(' ** CAUTION ** CDJOEL orthogonality tests OV01,OV02 & OV03:',3(1Pd9.1))

end subroutine CDJOELas
