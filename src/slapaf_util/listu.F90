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
! Copyright (C) 1993, Roland Lindh                                     *
!               Giovanni Ghigo                                         *
!***********************************************************************

subroutine ListU(Lu,Lbl,gq,mInt,nIter)
!***********************************************************************
!                                                                      *
! Object: to print last gradient or internal coordinate list           *
!                                                                      *
! Called from: RlxCtl                                                  *
!                                                                      *
!   Author: Giovanni Ghigo                                             *
!   Adapted from List by Roland Lindh, Dept. of Theoretical Chemistry, *
!             University of Lund, SWEDEN                               *
!***********************************************************************

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lu, mInt, nIter
character(len=8), intent(in) :: Lbl(mInt)
real(kind=wp), intent(in) :: gq(mInt,nIter)
integer(kind=iwp) :: igq

write(Lu,*)
write(Lu,*) '****************************'
write(Lu,*) '* Value of internal forces *'
write(Lu,*) '----------------------------'
do igq=1,mInt
  write(Lu,'(1X,A8,1X,F9.5)') Lbl(igq),gq(igq,nIter)
end do
write(Lu,*)

return

end subroutine ListU
