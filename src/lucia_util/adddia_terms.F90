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
! Copyright (C) 2011, Jeppe Olsen                                      *
!               2011, Giovanni Li Manni                                *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine ADDDIA_TERMS(NAEL,IASTR,NBEL,IBSTR,NORB,CVEC,SVEC,NSMST,H,XB,RJ,RK,NSSOA,NSSOB,ECORE,NTOOB,RJKAA,IASPGP,IASM,IBSPGP, &
                        IBSM,FACTOR)
! Update Sigma vector with diagonal terms for a given block
!     SVEC(IASPGP,IBSPGP) = SVEC(IASPGP,IBSPGP)
!                         + FACTOR*DIAG(IASPGP,IBSPGP)CVEC(IASPGP,IBSPGP)
! ========================
! General symmetry version
! ========================
!
! Jeppe Olsen and Giovanni Li Manni, September 2011
!
! I12 = 1 => only one-body part
!     = 2 =>      one+two-body part

use Constants, only: Zero, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NAEL, NBEL, NORB, NSMST, NSSOA(NSMST,*), NSSOB(NSMST,*), NTOOB, IASPGP, IASM, IBSPGP, IBSM
integer(kind=iwp), intent(inout) :: IASTR(NAEL,*), IBSTR(NBEL,*)
real(kind=wp), intent(in) :: CVEC(*), H(NORB), RJ(NTOOB,NTOOB), ECORE, FACTOR
real(kind=wp), intent(inout) :: SVEC(*), RK(NTOOB,NTOOB)
real(kind=wp), intent(out) :: XB(NORB)
real(kind=wp), intent(_OUT_) :: RJKAA(*)
integer(kind=iwp) :: I12, IA, IAEL, IB, IBEL, IDET, IDUM(1), IEL, JEL, NASTR1, NBSTR1, NIA, NIB
real(kind=wp) :: EAA, EB, HB, RJBB, X

IDUM(1) = 0
I12 = 2
!if (LUIN > 0) rewind(LUIN)
!if (LUOUT > 0) rewind(LUOUT)

#ifdef _DEBUGPRINT_
write(u6,*) ' ======================'
write(u6,*) ' ADDDIA_TERMS in action'
write(u6,*) ' ======================'
write(u6,*)
write(u6,*) ' IASM, IASPGP, IBSM, IBSPGP = ',IASM,IASPGP,IBSM,IBSPGP

write(u6,*) ' Diagonal one electron integrals'
call WRTMAT(H,1,NORB,1,NORB)
if (I12 == 2) then
  write(u6,*) ' Coulomb and exchange integrals'
  call WRTMAT(RJ,NORB,NORB,NTOOB,NTOOB)
  write(u6,*)
  call WRTMAT(RK,NORB,NORB,NTOOB,NTOOB)
end if
write(u6,*) ' FACTOR = ',FACTOR
#endif

!*3 Diagonal elements according to Handys formulae
!   (corrected for error)
!
!   DIAG(IDET) = HII*(NIA+NIB)
!              + 0.5 * ( J(I,J)-K(I,J) ) * N(I,A)*N(J,A)
!              + 0.5 * ( J(I,J)-K(I,J) ) * N(I,B)*N(J,B)
!              +         J(I,J) * N(I,A)*N(J,B)
! N(X) are occupation numbers

! K goes to J - K
if (I12 == 2) RK(:,:) = RJ(:,:)-RK(:,:)

! Construct array RJKAA(*) =   SUM(I) H(I)*N(I) +
!                          0.5*SUM(I,J) (J(I,J) - K(I,J))*N(I)*N(J)
!
! Obtain alpha strings of sym IASM and type IASPGP
call GETSTR_TOTSM_SPGP(1,IASPGP,IASM,NAEL,NASTR1,IASTR,NORB,0,IDUM,IDUM)

NIA = NSSOA(IASM,IASPGP)

#ifdef _DEBUGPRINT_
write(u6,*) ' After GETSTR for A strings'
write(u6,*) ' alpha strings obtained'
call IWRTMA(IASTR,NAEL,NIA,NAEL,NIA)
#endif

do IA=1,NIA
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
  RJKAA(IA) = EAA
end do
! Obtain alpha strings of sym IBSM and type IBTP
call GETSTR_TOTSM_SPGP(2,IBSPGP,IBSM,NBEL,NBSTR1,IBSTR,NORB,0,IDUM,IDUM)
NIB = NSSOB(IBSM,IBSPGP)
IDET = 0
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
  do IA=1,NSSOA(IASM,IASPGP)
    IDET = IDET+1
    X = EB+RJKAA(IA)
    do IEL=1,NAEL
      X = X+XB(IASTR(IEL,IA))
    end do
    SVEC(IDET) = SVEC(IDET)+CVEC(IDET)*(X+FACTOR)
  end do ! IA
end do ! IB

#ifdef _DEBUGPRINT_
write(u6,*) ' Input and output vectord, ADDDIA_TERMS'
call WRTMAT(CVEC,1,IDET,1,IDET)
call WRTMAT(SVEC,1,IDET,1,IDET)
#endif

end subroutine ADDDIA_TERMS
