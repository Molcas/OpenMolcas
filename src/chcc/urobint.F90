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

implicit none
#include "chcc1.fh"
#include "o2v4.fh"
integer maxa, maxbe, n, lun
real*8 W(1)
! help variables
integer length, i

!open(unit=Lun,file='IntFil',form='unformatted')

!1 def legth

length = maxa*maxa*maxbe*maxbe

!2 cycle over N

do i=1,n

  !2.1 full W with random numbers

  call RNFill(length,W(1),1.0d-2)

  !2.2 write block

  write(6,*) 'Vint',i,length
  call wri_chcc(lun,length,W(1))

end do

rewind(lun)

return

end subroutine UrobInt
