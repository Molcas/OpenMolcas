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

subroutine FrzDel(NREMO,IREMO,EOCC,E1,NREME,IREME,EEXT,E2,CMO,CMO1,ISEQ)
! MOLPT2(MOLCAS) SUBROUTINE
! CODED AJS, MAR. 15, 1990
!
! THIS SUBROUTINE IS USED TO MOVE THE ADDITIONALLY FROZEN AND
! ADDITIONALLY DELETED ORBITALS TO THE BOTTOM OR TO THE TOP
! RESPECTIVELY, OF THE EIGENVECTOR LIST. THE ORBITAL ENERGIES='PRIN'
! ARE REARRANGED ACCORDINGLY

#include "intent.fh"

use MBPT2_Global, only: nBas
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NREMO(*), IREMO(8,*), NREME(*), IREME(8,*)
real(kind=wp), intent(_OUT_) :: EOCC(*), EEXT(*), CMO(*)
real(kind=wp), intent(in) :: E1(*), E2(*), CMO1(*)
integer(kind=iwp), intent(_OUT_) :: ISEQ(*)
integer(kind=iwp) :: I, IAD, IAD0, IADD, IADE, IADF, IADO, IADR, IEE, IENE, IENO, IEO, ISYM, J
#include "corbinf.fh"

IAD0 = 1
IADR = 1
IEO = 0
IENO = 0
IENE = 0
IEE = 0
do ISYM=1,NSYM
  IADF = IADR
  IADO = IADF+NBAS(ISYM)*(NFRO(ISYM)+NREMO(ISYM))
  IADE = IADF+NBAS(ISYM)*(NFRO(ISYM)+NOCC(ISYM))
  IADD = IADF+NBAS(ISYM)*(NBAS(ISYM)-NDEL(ISYM)-NREME(ISYM))
  do I=1,NBAS(ISYM)
    ISEQ(I) = I
  end do
  do I=1,NFRO(ISYM)
    ISEQ(I) = 0
  end do
  do I=NBAS(ISYM),NBAS(ISYM)-NDEL(ISYM)+1,-1
    ISEQ(I) = 0
  end do
  do I=1,NREMO(ISYM)
    J = IREMO(ISYM,I)
    ISEQ(J) = 0
  end do
  do I=1,NREME(ISYM)
    J = IREME(ISYM,I)+NFRO(ISYM)+NOCC(ISYM)
    ISEQ(J) = 0
  end do

  do I=1,NFRO(ISYM)+NOCC(ISYM)
    IAD = IAD0+(I-1)*NBAS(ISYM)
    if (ISEQ(I) == 0) then
      ! frozen, move to appropriate place in symmetry block
      ! observe the sequence: frozen/occupied/external/deleted
      call DCOPY_(NBAS(ISYM),CMO1(IAD),1,CMO(IADF),1)
      IADF = IADF+NBAS(ISYM)
    else
      ! occupied...
      call DCOPY_(NBAS(ISYM),CMO1(IAD),1,CMO(IADO),1)
      IADO = IADO+NBAS(ISYM)
      IENO = IENO+1
      EOCC(IENO) = E1(IEO+I-NFRO(ISYM))
    end if
  end do
  do I=NFRO(ISYM)+NOCC(ISYM)+1,NBAS(ISYM)
    IAD = IAD0+(I-1)*NBAS(ISYM)
    if (ISEQ(I) == 0) then
      ! delete, move to appropriate place in symmetry block
      ! observe the sequence: frozen/occupied/external/deleted
      call DCOPY_(NBAS(ISYM),CMO1(IAD),1,CMO(IADD),1)
      IADD = IADD+NBAS(ISYM)
    else
      ! external...
      call DCOPY_(NBAS(ISYM),CMO1(IAD),1,CMO(IADE),1)
      IADE = IADE+NBAS(ISYM)
      IENE = IENE+1
      EEXT(IENE) = E2(IEE+I-NFRO(ISYM)-NOCC(ISYM))
    end if
  end do
  IADR = IADR+NBAS(ISYM)**2
  IAD0 = IADR
  IEO = IEO+NOCC(ISYM)
  IEE = IEE+NEXT(ISYM)

  ! UPDATE NFRO, NOCC, NEXT, NDEL, AND NORB

  NFRO(ISYM) = NFRO(ISYM)+NREMO(ISYM)
  NOCC(ISYM) = NOCC(ISYM)-NREMO(ISYM)
  NDEL(ISYM) = NDEL(ISYM)+NREME(ISYM)
  NEXT(ISYM) = NEXT(ISYM)-NREME(ISYM)
  NORB(ISYM) = NBAS(ISYM)-NFRO(ISYM)-NDEL(ISYM)
end do

return

end subroutine FrzDel
