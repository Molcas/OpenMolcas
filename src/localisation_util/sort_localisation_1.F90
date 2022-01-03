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
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Sort_Localisation_1(CMO,U,nBas,nOcc)
! Thomas Bondo Pedersen, November 2005.
!
! Purpose: sort CMO columns according to U.

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nBas, nOcc
real(kind=wp) :: CMO(nBas,nOcc), U(nOcc,nOcc)
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ip1, ip2, ipC, ipI1, ipI2, j, jmax, kOff, lC, lI1, lI2
real(kind=wp) :: Umax, Utst

! Allocations.
! ------------

lI1 = nOcc
lI2 = nOcc
lC = nBas*nOcc
call GetMem('Sr1I1','Allo','Inte',ipI1,lI1)
call GetMem('Sr1I2','Allo','Inte',ipI2,lI2)
call GetMem('Sr1C','Allo','Real',ipC,lC)

! Find max U element in each row.
! -------------------------------

ip1 = ipI1-1
do i=1,nOcc
  iWork(ip1+i) = i
end do

ip2 = ipI2-1
do i=1,nOcc
  jmax = 0
  Umax = -huge(Umax)
  do j=1,nOcc
    if (iWork(ipI1-1+j) == j) then
      Utst = abs(U(i,j))
      if (Utst > Umax) then
        jmax = j
        Umax = Utst
      end if
    end if
  end do
  if (jmax == 0) then
    call SysAbendMsg('Sort_Localisation_1','Error:','jmax=0')
  else
    iWork(ip1+jmax) = 0
    iWork(ip2+i) = jmax
  end if
end do

! Swap MOs according to I2.
! -------------------------

call dCopy_(nBas*nOcc,CMO,1,Work(ipC),1)
do i=1,nOcc
  kOff = ipC+nBas*(iWork(ipI2-1+i)-1)
  call dCopy_(nBas,Work(kOff),1,CMO(1,i),1)
end do

! De-allocate.
! ------------

call GetMem('Sr1C','Free','Real',ipC,lC)
call GetMem('Sr1I2','Free','Inte',ipI2,lI2)
call GetMem('Sr1I1','Free','Inte',ipI1,lI1)

end subroutine Sort_Localisation_1
