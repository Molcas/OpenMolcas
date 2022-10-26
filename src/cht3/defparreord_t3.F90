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
! define parameters in cht3_global using NaGrpR,maxdim
!
! I/O parameter description:
! NxGrpR   - # of groups in a (=b) set (I)
! maxdim   - # maximal dimension of a (=b) Groups(O)

use ChT3_global, only: DimGrpaR, L1Name, L2Name, nv, T2Name
use stdalloc, only: mma_allocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NaGrpR
integer(kind=iwp), intent(out) :: maxdim
integer(kind=iwp) :: i, j, Low, Up_prev, Up
real(kind=wp) :: rdim

call mma_allocate(DimGrpaR,NaGrpR,label='DimGrpaR')
call mma_allocate(L1Name,NaGrpR,label='L1Name')
call mma_allocate(L2Name,NaGrpR,NaGrpR,label='L2Name')
call mma_allocate(T2Name,NaGrpR,NaGrpR,label='T2Name')

!1 define parameters of Groups of a set

rdim = real(nv,kind=wp)/real(NaGrpR,kind=wp)

Up_prev = 0
do i=1,NaGrpR

  Low = Up_prev+1
  if (i == NaGrpR) then
    Up = nv
  else
    Up = int(rdim*i)
  end if
  Up_prev = Up

  DimGrpaR(i) = Up-Low+1

end do

!2 find maximal dimensions of a'

maxdim = DimGrpaR(1)
do i=1,NaGrpR
  if (DimGrpaR(i) > maxdim) maxdim = DimGrpaR(i)
end do

!3 def L1Name, L2Name, T2Name

do i=1,NaGrpR
  write(L1Name(i),'(A4,I0.2)') 'L1vc',i
  do j=1,NaGrpR
    write(L2Name(i,j),'(A2,I0.2,I0.2)') 'L2',i,j
    write(T2Name(i,j),'(A2,I0.2,I0.2)') 'T2',i,j
  end do
end do

return

end subroutine DefParReord_t3
