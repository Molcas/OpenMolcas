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

subroutine my_block(vblock,vblock_my)
! this subroutine calculates maximum overlap between juzek's
! vblock segmentation and palo's dimgrp

use ChT3_global, only: DimGrpaR, nv, NvGrp
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: vblock
integer(kind=iwp), intent(out) :: vblock_my
integer(kind=iwp) :: i, i_f, i_l, i_tmp, j, pos, vblock_my_tmp
logical(kind=iwp) :: found

vblock_my = 0
i_l = 0
i_f = 0

do i=1,nv,vblock

  ! - find initial position of the i-th juzek's block

  pos = 0
  found = .false.
  do j=1,NvGrp
    pos = pos+DimGrpaR(j)
    if ((i <= pos) .and. (.not. found)) then
      i_f = j
      found = .true.
      !mp write(u6,'(A,3(i5,x))') 'i,i_f,pos     = ',i,i_f,pos
    end if
  end do

  if ((i+vblock-1) <= nv) then
    i_tmp = i+vblock-1
  else
    i_tmp = nv
  end if
  !mp write(u6,'(A,2(i5,x))') 'i,i_tmp        = ',i,i_tmp

  ! - find terminal position of the i-th juzek's block

  pos = 0
  found = .false.
  do j=1,NvGrp
    pos = pos+DimGrpaR(j)
    if ((i_tmp <= pos) .and. (.not. found)) then
      i_l = j
      found = .true.
      !mp write(u6,'(A,3(i5,x))') 'i_tmp,i_l,pos = ',i_tmp,i_l,pos
    end if
  end do

  vblock_my_tmp = 0
  do j=i_f,i_l
    vblock_my_tmp = vblock_my_tmp+DimGrpaR(j)
  end do

  if (vblock_my_tmp > vblock_my) vblock_my = vblock_my_tmp

  !mp write(u6,'(A,2(i5,x))') 'vblock_my_tmp, vblock_my',vblock_my_tmp,vblock_my
  !mp write(u6,*)
end do

return

end subroutine my_block
