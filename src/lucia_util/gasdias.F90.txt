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
!               2015, Lasse Kragh Soerensen                            *
!***********************************************************************

subroutine GASDIAS(NAEL,IASTR,NBEL,IBSTR,NORB,DIAG,NSMST,H,XB,RJ,RK,NSSOA,NSSOB,LUDIA,ECORE,PSSIGN,NTOOB,ICISTR,RJKAA,I12,IBLTP, &
                   NBLOCK,IBLKFO,I_AM_OUT,N_ELIMINATED_BATCHES)
! Calculate determinant diagonal
! Turbo-ras version
!
! Driven by IBLKFO, May 97
!
! ========================
! General symmetry version
! ========================
!
! Jeppe Olsen, July 1995, GAS version
!
! I12 = 1 => only one-body part
!     = 2 =>      one+two-body part
!
! Added the possibility to zero out unwanted diagonal part
! which is needed for highly excited states. Lasse 2015

use lucia_data, only: IDISK
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NAEL, NBEL, NORB, NSMST, NSSOA(NSMST,*), NSSOB(NSMST,*), LUDIA, NTOOB, ICISTR, I12, IBLTP(*), &
                                 NBLOCK, IBLKFO(8,NBLOCK), I_AM_OUT(*), N_ELIMINATED_BATCHES
integer(kind=iwp), intent(inout) :: IASTR(NAEL,*), IBSTR(NBEL,*)
real(kind=wp), intent(_OUT_) :: DIAG(*), RJKAA(*)
real(kind=wp), intent(in) :: H(NORB), RJ(NTOOB,NTOOB), ECORE, PSSIGN
real(kind=wp), intent(out) :: XB(NORB)
real(kind=wp), intent(inout) :: RK(NTOOB,NTOOB)
integer(kind=iwp) :: I, I_AM_NOT_WANTED, IA, IAEL, IASM, IASTOP, IASTRT, IATP, IB, IBEL, IBLK, IBSM, IBSTOP, IBSTRT, IBTP, IDET, &
                     IDUM_ARR(1), IEL, IOFF, IREST1, ITDET, JEL, NASTR1, NBSTR1
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: IBLOCK, II
#endif
real(kind=wp) :: EAA, EB, HB, RJBB, X, XADD

if (PSSIGN == -One) then
  XADD = 1.0e6_wp
else
  XADD = Zero
end if

#ifdef _DEBUGPRINT_
write(u6,*) ' Diagonal one electron integrals'
call WRTMAT(H,1,NORB,1,NORB)
write(u6,*) ' Core energy ',ECORE
if (I12 == 2) then
  write(u6,*) ' Coulomb and exchange integrals'
  call WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
  write(u6,*)
  call WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
end if

write(u6,*) ' TTSS for Blocks'
do IBLOCK=1,NBLOCK
  write(u6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
end do

write(u6,*) ' I12 = ',I12
#endif

! Diagonal elements according to Handys formulae
! (corrected for error)
!
! DIAG(IDET) = HII*(NIA+NIB)
!            + 0.5 * ( J(I,J)-K(I,J) ) * NIA*NJA
!            + 0.5 * ( J(I,J)-K(I,J) ) * NIB*NJB
!            +         J(I,J) * NIA*NJB
!
! K goes to J - K
if (I12 == 2) RK(:,:) = RJ(:,:)-RK(:,:)
IDET = 0
ITDET = 0
if (LUDIA /= 0) IDISK(LUDIA) = 0

do IBLK=1,NBLOCK
  ! Lasse addition
  I_AM_NOT_WANTED = 0
  do I=1,N_ELIMINATED_BATCHES
    if (I_AM_OUT(I) == IBLK) then
      I_AM_NOT_WANTED = 1
      exit
    end if
  end do
  ! Lasse addition end

  IATP = IBLKFO(1,IBLK)
  IBTP = IBLKFO(2,IBLK)
  IASM = IBLKFO(3,IBLK)
  IBSM = IBLKFO(4,IBLK)

  if (IBLTP(IASM) == 2) then
    IREST1 = 1
  else
    IREST1 = 0
  end if

  ! Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
  !                          0.5*SUM(I,J) (J(I,J) - K(I,J))*N(I)*N(J)

  ! Obtain alpha strings of sym IASM and type IATP
  IDUM_ARR = 0
  call GETSTR_TOTSM_SPGP(1,IATP,IASM,NAEL,NASTR1,IASTR,NORB,0,IDUM_ARR,IDUM_ARR)
  IOFF = 1
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
  ! Obtain beta strings of sym IBSM and type IBTP
  call GETSTR_TOTSM_SPGP(2,IBTP,IBSM,NBEL,NBSTR1,IBSTR,NORB,0,IDUM_ARR,IDUM_ARR)
  IBSTRT = 1
  IBSTOP = NSSOB(IBSM,IBTP)
  do IB=IBSTRT,IBSTOP

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

    if ((IREST1 == 1) .and. (IATP == IBTP)) then
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
      ! Lasse addition
      if (I_AM_NOT_WANTED == 0) then
        if (IB == IA) then
          DIAG(IDET) = X+XADD
        else
          DIAG(IDET) = X
        end if
      else
        DIAG(IDET) = Zero
      end if
      ! Lasse addition end
    end do
    ! End of loop over alpha strings|
  end do
  ! End of loop over betastrings
  ! Yet a RAS block of the diagonal has been constructed
  if (ICISTR >= 2) then
#   ifdef _DEBUGPRINT_
    write(u6,*) ' number of diagonal elements to disc ',IDET
    call WRTMAT(DIAG,1,IDET,1,IDET)
#   endif
    IDUM_ARR(1) = IDET
    call ITODS(IDUM_ARR,1,-1,LUDIA)
    call TODSC(DIAG,IDET,-1,LUDIA)
    IDET = 0
  end if
end do
! End of loop over blocks

#ifdef _DEBUGPRINT_
write(u6,*) ' Number of diagonal elements generated (1)',ITDET
#endif

#ifdef _DEBUGPRINT_
if (ICISTR <= 1) then
  write(u6,*) ' CIDIAGONAL'
  call WRTMAT(DIAG,1,IDET,1,IDET)
end if
#endif

if (ICISTR >= 2) then
  IDUM_ARR(1) = -1
  call ITODS(IDUM_ARR,1,-1,LUDIA)
end if

end subroutine GASDIAS
