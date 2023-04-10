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

subroutine DefParo3v3(NvGrp,maxdim)
! This routine does:
! define parameters in o3v3.fh using NvGrp,maxdim
!
! I/O parameter description:
! NvGrp    - # of groups in a,b,be,ga set (I)
! maxdim   - # maximal dimension of (a,b,be,ga)" Groups(O)

implicit none
#include "chcc1.fh"
#include "o3v3.fh"
#include "chcc_files.fh"
integer NvGrp, maxdim
! help variables
real*8 rdim
integer i, j
integer Up(1:MaxGrp), Low(1:MaxGrp)

!1 define parameters of Groups of v set

rdim = 1.0d0*nv/(1.0d0*NvGrp)

do i=1,NvGrp

  if (i == 1) then
    Up(i) = int(rdim*i)
    Low(i) = 1
  else if (i == NvGrp) then
    Up(i) = nv
    Low(i) = Up(i-1)+1
  else
    Up(i) = int(rdim*i)
    Low(i) = Up(i-1)+1
  end if

  DimGrpv(i) = (Up(i)-Low(i))+1

end do

!5 find maximal dimensions of v'

maxdim = DimGrpv(1)
do i=1,NvGrp
  if (DimGrpv(i) > maxdim) maxdim = DimGrpv(i)
end do

!6.1 def L2Name, T2Name, I2Name,I3Name,Tmp1Name,Tmp2Name

do i=1,MaxGrp
  do j=1,MaxGrp
    call DefParo3v3Hlp1(i,j,'L2',L2Name(i,j))
    call DefParo3v3Hlp1(i,j,'T2',T2Name(i,j))
    call DefParo3v3Hlp1(i,j,'I2',I2Name(i,j))
    call DefParo3v3Hlp1(i,j,'I3',I3Name(i,j))
    call DefParo3v3Hlp1(i,j,'X1',Tmp1Name(i,j))
    call DefParo3v3Hlp1(i,j,'X2',Tmp2Name(i,j))
  end do
end do

!6.2 def L1Name,I1Name

do i=1,MaxGrp
  call DefParo3v3Hlp2(i,'L1vc',L1Name(i))
  call DefParo3v3Hlp2(i,'I1in',I1Name(i))
end do

return

end subroutine DefParo3v3
