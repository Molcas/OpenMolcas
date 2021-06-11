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

subroutine writedrt(id)

use gugadrt_global, only: ja, jb, jd, jj, jm, js, jt, jv, kk, ludrt, no, norb_inn
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: id
integer(kind=iwp) :: i, idisk, idum(1), idx(2), jbuf(4*(id+1)), nc

nc = 1
do i=0,id
  jbuf(nc:nc+3) = jj(1:4,i)
  nc = nc+4
end do
idisk = 0
! number of nodes
call idafile(ludrt,2,idx,2,idisk)
idisk = idx(2)
idum(1) = id
call idafile(ludrt,1,idum,1,idisk)
call idafile(ludrt,1,ja,id,idisk)
call idafile(ludrt,1,jb,id,idisk)
call idafile(ludrt,1,jm,id,idisk)
call idafile(ludrt,1,jbuf,4*(id+1),idisk)
call idafile(ludrt,1,kk(0),1+id,idisk)
call idafile(ludrt,1,no(0),norb_inn+2,idisk)
idum(1) = jv
call idafile(ludrt,1,idum,1,idisk)
call idafile(ludrt,1,jd,8,idisk)
call idafile(ludrt,1,jt,8,idisk)
call idafile(ludrt,1,js,8,idisk)

return

end subroutine writedrt
