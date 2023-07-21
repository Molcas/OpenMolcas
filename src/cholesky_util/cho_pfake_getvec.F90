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
! Copyright (C) 2007, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_PFake_GetVec(Vec,lVec,IDV,lIDV,InfV,iSym,nRead,iRedC)

implicit none
integer lVec, lIDV
real*8 Vec(lVec)
integer IDV(lIDV)
integer InfV(2,*)
integer iSym, nRead, iRedC
character*16 SecNam
parameter(SecNam='Cho_PFake_GetVec')
integer ipV, Mem, iVec, n, m

nRead = 0
ipV = 1
Mem = lVec
do iVec=1,lIDV
  n = 0
  m = 0
  call Cho_VecRd(Vec(ipV),Mem,IDV(iVec),IDV(iVec),iSym,n,iRedC,m)
  if (n == 1) then
    nRead = nRead+1
    ipV = ipV+m
    Mem = Mem-m
    InfV(1,iVec) = m
  else if (n == 0) then
    return
  else
    call Cho_Quit('Logical error in '//SecNam,103)
  end if
end do

end subroutine Cho_PFake_GetVec
