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

subroutine Union(iU,nU,iV,nV,iR,iM,nM)

implicit real*8(A-H,O-Z)
integer iU(nU), iV(nV), iM(8)
logical RinT_

! M is formed as U union RU

call iCopy(nU,iU,1,iM,1) ! copy the first elements
nM = nU
do i=1,nV
  iRV = ieor(iR,iV(i))
  if (.not. RinT_(iM,nM,iRV)) then
    nM = nM+1
    iM(nM) = iRV
  end if
end do

return

end subroutine Union
