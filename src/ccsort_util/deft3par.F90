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

subroutine DefT3par(noa,nsym)
! this routine does:
! 0) Open t3nam file, with lunt3 lun
! define parameters, required for T3 integral handling, namely
! 1) def T3IndPos(i)
!    address positions for all occupied orbitals in t3nam file
! 2) def T3Off(ii,isym)
!    relative shifts of address for ii-th block of R_i(a,bc)
!    for each symmetry
!
! noa   - array with occupation numbers
! nsym  - actual number of irreps

use ccsort_global, only: daddr, lunt3, mapdri, mapiri, mbas, posri0
use CCT3_global, only: T3IntPos, t3nam, T3Off
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: noa(8), nsym
integer(kind=iwp) :: i, idum(1), ii, iorb, length, post, symi
real(kind=wp) :: dum(1)

!0  open t3nam file
!lunt3= 1
call daname(lunt3,t3nam)

!1 set address poiter to 0
daddr(lunt3) = 0

!2 first record in t3nam file is T3IntPos
!  (emulate writing of T3IntPos)
idum(1) = 0
dum(1) = Zero
call idafile(lunt3,0,idum,mbas,daddr(lunt3))

iorb = 0
!3 cycle over irreps
do symi=1,nsym

  !3.1 make mapd and mapi for  R_i(a,bc)
  call ccsort_t3grc0(3,8,4,4,4,0,symi,posri0,post,mapdri,mapiri)

  !3.2 cycle over occupied orbitals in symi
  do i=1,noa(symi)

    !3.2.1 save initial address for this orbital
    iorb = iorb+1
    T3IntPos(iorb) = daddr(lunt3)

    !3.2.2 emulate writing of mapd and mapp
    call idafile(lunt3,0,idum,513*6,daddr(lunt3))
    call idafile(lunt3,0,idum,8*8*8,daddr(lunt3))

    !3.2.3  cycle over all blocks of R_i(a,bc), which will be stored separately
    do ii=1,mapdri(0,5)

      !3.2.3.1 def T3Off(ii,symi)
      !        note, that iorb is always proper one, since only besides
      !        first occ. orbital in given irrep T3Off is defined
      if (i == 1) T3Off(ii,symi) = daddr(lunt3)-T3IntPos(iorb)

      !3.2.3.2 emulate writing of each block
      length = mapdri(ii,2)
      call ddafile(lunt3,0,dum,length,daddr(lunt3))

    end do
  end do
end do

return

end subroutine DefT3par
