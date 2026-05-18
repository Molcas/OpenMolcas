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

subroutine Get_Cholesky_Vectors(ITK,ITQ,JSYM,Array,mArray,nArray,IBSTA,IBEND)

use definitions, only: iwp, wp
use CHOVEC_IO, only: NPQ_CHOTYPE, NVLOC_CHOBATCH, IDLOC_CHOGROUP
use caspt2_global, only: LUDRA
use caspt2_module, only: NSYM

implicit none
integer(kind=iwp), intent(in) :: ITK, ITQ, JSYM, IBSTA, IBEND
integer(kind=iwp), intent(in) :: mArray
integer(kind=iwp), intent(out) :: nArray
real(kind=wp), intent(out) :: Array(mArray)
integer(kind=iwp) ICASE, LKETSM, ISYK, NQK, IB, NV, NKETSM, IDISK

! ugly hack to convert separate k/q orbital types into a specific case
ICASE = ITK*ITQ
if (ICASE == 3) then
  ICASE = 4
else
  ICASE = ICASE/2
end if

LKETSM = 1
do ISYK=1,NSYM
  NQK = NPQ_CHOTYPE(ICASE,ISYK,JSYM)
  if (NQK == 0) cycle
  do IB=IBSTA,IBEND
    NV = NVLOC_CHOBATCH(IB)
    NKETSM = NQK*NV
    IDISK = IDLOC_CHOGROUP(ICASE,ISYK,JSYM,IB)
    call DDAFILE(LUDRA,2,Array(LKETSM),NKETSM,IDISK)
    LKETSM = LKETSM+NKETSM
  end do
end do
nArray = LKETSM-1

end subroutine Get_Cholesky_Vectors
