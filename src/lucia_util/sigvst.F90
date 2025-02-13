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

subroutine SIGVST(ISGVST,NSMST)
! Obtain ISGVST(ISM) : Symmetry of sigma v on string of symmetry ism

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSMST
integer(kind=iwp), intent(out) :: ISGVST(NSMST)
integer(kind=iwp) :: IML, IPARI, ISM, ISM_, MIML, MISM, NTEST

do ISM=1,NSMST
  ISM_ = ISM
  call MLSM(IML,IPARI,ISM_,'ST',2)
  !    MLSM(IML,IPARI,ISM,TYPE,IWAY)
  MIML = -IML
  call MLSM(MIML,IPARI,MISM,'ST',1)
  ISGVST(ISM) = MISM
end do

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' ISGVST array'
  write(u6,*) ' ============'
  call IWRTMA(ISGVST,1,NSMST,1,NSMST)
end if

end subroutine SIGVST
