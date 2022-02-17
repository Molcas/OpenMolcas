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

subroutine faibj3(NSIJ,IFT,AIBJ,FSEC,FAC,IIN,INS,IPOA,IPOF)

use mrci_global, only: NSYM, NVIR
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NSIJ, IFT, IPOF(9)
real(kind=wp), intent(in) :: AIBJ(*), FAC
real(kind=wp), intent(_OUT_) :: FSEC(*)
integer(kind=iwp), intent(inout) :: IIN
integer(kind=iwp), intent(out) :: INS, IPOA(9)
integer(kind=iwp) :: IAB, IASYM, IBSYM

call IPO(IPOA,NVIR,MUL,NSYM,NSIJ,IFT)
! INTEGRAL COMBINATION APPROPRIATE FOR SINGLET-COUPLING:
do IASYM=1,NSYM
  IBSYM = MUL(NSIJ,IASYM)
  if (IBSYM > IASYM) cycle
  IAB = IPOA(IASYM+1)-IPOA(IASYM)
  if (IAB == 0) cycle
  if (NSIJ == 1) then
    call SECEQ(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),FSEC(IIN+1),NVIR(IASYM),0,FAC)
  else
    call SECNE(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),FSEC(IIN+1),NVIR(IASYM),NVIR(IBSYM),0)
  end if
  IIN = IIN+IAB
end do
INS = IIN
! INTEGRAL COMBINATION APPROPRIATE FOR TRIPLET-COUPLING:
do IASYM=1,NSYM
  IBSYM = MUL(NSIJ,IASYM)
  if (IBSYM > IASYM) cycle
  IAB = IPOA(IASYM+1)-IPOA(IASYM)
  if (IAB == 0) cycle
  if (NSIJ == 1) then
    call SECEQ(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),FSEC(IIN+1),NVIR(IASYM),1,FAC)
  else
    call SECNE(AIBJ(IPOF(IASYM)+1),AIBJ(IPOF(IBSYM)+1),FSEC(IIN+1),NVIR(IASYM),NVIR(IBSYM),1)
  end if
  IIN = IIN+IAB
end do

return

end subroutine faibj3
