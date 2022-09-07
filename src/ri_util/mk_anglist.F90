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
! Copyright (C) 2007,2008, Roland Lindh                                *
!***********************************************************************

subroutine Mk_AngList(iAL,nCompA,nCompB,iD_c,nD_c,List2,nList2,mData,iAng,jAng)

use Definitions, only: iwp

implicit none
integer(kind=iwp) :: nCompA, nCompB, iAL(nCompA,nCompB), nD_c, iD_c(nD_c), nList2, mData, List2(mData,nList2), iAng, jAng
integer(kind=iwp) :: ijSO, jD_c

call IZero(iAL,nCompA*nCompB)
do jD_c=1,nD_c
  ijSO = iD_c(jD_c)
  if ((List2(1,ijSO) == iAng) .and. (List2(2,ijSO) == jAng)) then
    iAL(List2(3,ijSO),List2(4,ijSO)) = 1
  end if
end do

return

end subroutine Mk_AngList
