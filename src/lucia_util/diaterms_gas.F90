!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DIATERMS_GAS(NAEL,IASTR,NBEL,IBSTR,NORB,VEC,NSMST,H,IDC,XB,RJ,RK,NSSOA,NSSOB,ECORE,LUIN,LUOUT,NTOOB,RJKAA,I12,IBLOCK, &
                        NBLOCK,ITASK,FACTOR,I0CHK,I0BLK)
! Terms from diagonal to specific blocks
!
! Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
! Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
!
! Calculate determinant diagonal
!
! ========================
! General symmetry version
! ========================
!
! Jeppe Olsen, July 1995, GAS version
!
! I12 = 1 => only one-body part
!     = 2 =>      one+two-body part

use Index_Functions, only: nTri_Elem
use lucia_data, only: IDISK
use Constants, only: Zero, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NAEL, NBEL, NORB, NSMST, IDC, NSSOA(NSMST,*), NSSOB(NSMST,*), LUIN, LUOUT, NTOOB, I12, &
                                 IBLOCK(8,*), NBLOCK, ITASK, I0CHK, I0BLK(*)
integer(kind=iwp), intent(inout) :: IASTR(NAEL,*), IBSTR(NBEL,*)
real(kind=wp), intent(_OUT_) :: VEC(*), RJKAA(*)
real(kind=wp), intent(in) :: H(NORB), RJ(NTOOB,NTOOB), ECORE, FACTOR
real(kind=wp), intent(out) :: XB(NORB)
real(kind=wp), intent(inout) :: RK(NTOOB,NTOOB)
integer(kind=iwp) :: IA, IAEL, IAMPACK, IASM, IASTOP, IASTRT, IATP, IB, IBEL, IBSM, IBTP, IDET, IDUM_ARR(1), IEL, IMZERO, IOFF, &
                     IPACK, ITDET, JBLOCK, JEL, LDET, NASTR1, NBSTR1, NIA, NIB
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: NAST
#endif
real(kind=wp) :: EAA, EB, HB, RJBB, X


if (LUIN > 0) IDISK(LUIN) = 0
if (LUOUT > 0) IDISK(LUOUT) = 0

#ifdef _DEBUGPRINT_
write(u6,*) ' ======================'
write(u6,*) ' DIATERMS_GAS in action'
write(u6,*) ' ======================'
write(u6,*)
write(u6,*) ' LUIN,LUOUT = ',LUIN,LUOUT
write(u6,*) ' NBLOCK =',NBLOCK
write(u6,*) ' I0CHK = ',I0CHK

write(u6,*) ' Diagonal one electron integrals'
call WRTMAT(H,1,NORB,1,NORB)
if (I12 == 2) then
  write(u6,*) ' Coulomb and exchange integrals'
  call WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
  write(u6,*)
  call WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
  write(u6,*) ' I12 and ITASK = ',I12,ITASK
end if
write(u6,*) ' FACTOR = ',FACTOR
#endif

!*3 Diagonal elements according to Handys formulae
!   (corrected for error)
!
!   DIAG(IDET) = HII*(NIA+NIB)
!              + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
!              + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
!              +         J(I,J) * NIA*NJB
!
! K goes to J - K
if (I12 == 2) RK(:,:) = RJ(:,:)-RK(:,:)

ITDET = 0
IDET = 0
do JBLOCK=1,NBLOCK
  if (IBLOCK(1,JBLOCK) > 0) then
    IATP = IBLOCK(1,JBLOCK)
    IBTP = IBLOCK(2,JBLOCK)
    IASM = IBLOCK(3,JBLOCK)
    IBSM = IBLOCK(4,JBLOCK)
    IOFF = IBLOCK(6,JBLOCK)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' Block in action : IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
#   endif

    if ((IDC == 2) .and. (IASM == IBSM) .and. (IATP == IBTP)) then
      IPACK = 1
    else
      IPACK = 0
    end if

    ! Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
    !                          0.5*SUM(I,J) (J(I,J) - K(I,J))*N(I)*N(J)

    ! Obtain alpha strings of sym IASM and type IATP
    IDUM_ARR = 0
    call GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,NORB,0,IDUM_ARR,IDUM_ARR)
#   ifdef _DEBUGPRINT_
    write(u6,*) ' After GETSTR for A strings'
    write(u6,*) ' alpha strings obtained'
    NAST = NSSOA(IASM,IATP)
    call IWRTMA(IASTR,NAEL,NAST,NAEL,NAST)
#   endif

    IOFF = 1
    NIA = NSSOA(IASM,IATP)
    do IA=1,NSSOA(IASM,IATP)
      EAA = Zero
      do IEL=1,NAEL
        IAEL = IASTR(IEL,IA)
        EAA = EAA+H(IAEL)
        if (I12 == 2) then
          do JEL=1,NAEL
            EAA = EAA+Half*RK(IASTR(JEL,IA),IAEL)
          end do
        end if
      end do
      RJKAA(IA-IOFF+1) = EAA
    end do
    ! Obtain alpha strings of sym IBSM and type IBTP
    call GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,NORB,0,IDUM_ARR,IDUM_ARR)
    NIB = NSSOB(IBSM,IBTP)

    IMZERO = 0
    if (LUIN > 0) then
      call IDAFILE(LUIN,2,IDUM_ARR,1,IDISK(LUIN))
      LDET = IDUM_ARR(1)
      call IDAFILE(LUIN,2,IDUM_ARR,1,IDISK(LUIN))
      IDET = 0
      call FRMDSC(VEC,LDET,-1,LUIN,IMZERO,IAMPACK)
    end if

    if (I0CHK == 1) then
      IMZERO = I0BLK(JBLOCK)
      if (IMZERO == 1) then
        ! Update offset to next block
        if ((IPACK == 1) .and. (IATP == IBTP)) then
          IDET = IDET+nTri_Elem(NIA)
        else
          IDET = IDET+NIA*NIB
        end if
      end if
    end if
    !write(u6,*) ' DIATERMS_GAS : I0CHK,JBLOCK IMZERO',I0CHK,JBLOCK,IMZERO

    if (IMZERO /= 1) then
      ! Calculate ...

      do IB=1,NIB

        ! Terms depending only on IB

        HB = Zero
        RJBB = Zero
        XB(:) = Zero

        do IEL=1,NBEL
          IBEL = IBSTR(IEL,IB)
          HB = HB+H(IBEL)

          if (I12 == 2) then
            do JEL=1,NBEL
              RJBB = RJBB+RK(IBSTR(JEL,IB),IBEL)
            end do

            XB(:) = XB(:)+RJ(1:NORB,IBEL)
          end if
        end do
        EB = HB+Half*RJBB+ECORE

        if ((IPACK == 1) .and. (IATP == IBTP)) then
          IASTRT = IB
        else
          IASTRT = 1
        end if

        IASTOP = NSSOA(IASM,IATP)
        do IA=IASTRT,IASTOP
          IDET = IDET+1
          ITDET = ITDET+1
          X = EB+RJKAA(IA-IOFF+1)
          do IEL=1,NAEL
            X = X+XB(IASTR(IEL,IA))
          end do
          ! Obtain VEC = (DIAGONAL + FACTOR) ** -1 VEC (ITASK = 1)
          ! Obtain VEC = (DIAGONAL + FACTOR)       VEC (ITASK = 2)
          if (ITASK == 1) then
            if (abs(X+FACTOR) > 1.0e-10_wp) then
              VEC(IDET) = VEC(IDET)/(X+FACTOR)
            else
              VEC(IDET) = Zero
            end if
          else
            VEC(IDET) = VEC(IDET)*(X+FACTOR)
          end if
          !write(u6,*) ' IDET,X,VEC(IDET) ',IDET,X,VEC(IDET)
        end do
      end do
    end if

    if (LUOUT > 0) then
      IDUM_ARR(1) = LDET
      call ITODS(IDUM_ARR,1,-1,LUOUT)
      call TODSC(VEC,LDET,-1,LUOUT)
      !write(u6,*) ' Number of elements transferred to DISC ',LDET
      IDET = 0
    end if

  end if
end do

if (LUOUT > 0) then
  IDUM_ARR(1) = -1
  call ITODS(IDUM_ARR,1,-1,LUOUT)
end if

!write(u6,*) ' Mission DIATERMS finished'

end subroutine DIATERMS_GAS
