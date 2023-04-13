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

subroutine UrobInt(W,maxa,maxbe,n,lun)
! vyraba subor LunInt so simulovanymi (a"be"|b"ga") integralmi
! kde je za sebou N blokov

use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: W(1)
integer(kind=iwp) :: maxa, maxbe, n, lun
integer(kind=iwp) :: i, length

!open(unit=Lun,file='IntFil',form='unformatted')

!1 def length

length = maxa*maxa*maxbe*maxbe

!2 cycle over N

do i=1,n

  !2.1 full W with random numbers

  call RNFill(length,W(1),1.0e-2_wp)

  !2.2 write block

  write(u6,*) 'Vint',i,length
  call wri_chcc(lun,length,W(1))

end do

rewind(lun)

return

end subroutine UrobInt
