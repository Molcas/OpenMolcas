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

subroutine DefParReord_t3(NaGrpR,maxdim)
! This routine does:
! define parameters in cht3_reord.fh using NaGrpR,maxdim
!
! I/O parameter description:
! NxGrpR   - # of groups in a (=b) set (I)
! maxdim   - # maximal dimension of a (=b) Groups(O)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: NaGrpR, maxdim
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "files.fh"
integer i, j, Low(MaxGrp), Up(MaxGrp)
real(kind=wp) :: rdim

!1 define parameters of Groups of a set

rdim = real(nv,kind=wp)/real(NaGrpR,kind=wp)

do i=1,NaGrpR

  if (i == 1) then
    Up(i) = int(rdim*i)
    Low(i) = 1
  else if (i == NaGrpR) then
    Up(i) = nv
    Low(i) = Up(i-1)+1
  else
    Up(i) = int(rdim*i)
    Low(i) = Up(i-1)+1
  end if

  DimGrpaR(i) = (Up(i)-Low(i))+1

end do

!2 find maximal dimensions of a'

maxdim = DimGrpaR(1)
do i=1,NaGrpR
  if (DimGrpaR(i) > maxdim) then
    maxdim = DimGrpaR(i)
  end if
end do

!3.1 def L2Name, T2Name, I2Name,I3Name

do i=1,MaxGrp
  do j=1,MaxGrp
    write(L2Name(i,j),'(A2,I0.2,I0.2)') 'L2',i,j
    write(T2Name(i,j),'(A2,I0.2,I0.2)') 'T2',i,j
    write(I2Name(i,j),'(A2,I0.2,I0.2)') 'I2',i,j
    write(I3Name(i,j),'(A2,I0.2,I0.2)') 'I3',i,j
  end do
end do

!3.2 def L1Name,I1Name

do i=1,MaxGrp
  write(L1Name(i),'(A4,I0.2)') 'L1vc',i
  write(I1Name(i),'(A4,I0.2)') 'I1in',i
end do

!3.3 def L0Name,I0Name

L0Name = 'L0vctr'
I0Name = 'I0intg'

return

end subroutine DefParReord_t3
