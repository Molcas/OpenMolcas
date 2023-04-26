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

subroutine addinta(wrk,wrksize,syma,ammap)
! this routine does for all a in syma
! 1- reconstruct #2 <_a,m,p,q> from TEMPDA2 file
! 2- prepare corresponding <_am p q> (like <amef>aaaa) to #3
! and write it to open INTA1-4
! N.B.  this routine uses following foreign routines:
! dawrtmap
! dawri

use ccsort_global, only: luna1, luna2, luna3, luna4, map3, mbas, nva, nvb
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize, syma, ammap(mbas,8,8)
real(kind=wp), intent(_OUT_) :: wrk(wrksize)
integer(kind=iwp) :: a, lenefaaaa, lenefabab, lenefbaab, lenefbbbb, lenejaaaa, lenejabab, lenejabba, lenejbaab, lenejbaba, &
                     lenejbbbb, post, rc

! map2 of #2 <_a,m|p,q> are prepared

! make required map3 and write them to INTA1-4
! define lengths of this mediates

!1   to INTA1 <m,_a||ef>aaaa, <m,_a||ef>baab
call ccsort_grc0(3,2,1,3,3,0,syma,post,map3)
call deflength(map3,lenefaaaa)
call dawrtmap(luna1,map3,rc)
call ccsort_grc0(3,0,2,3,4,0,syma,post,map3)
call deflength(map3,lenefbaab)
call dawrtmap(luna1,map3,rc)

!2   to INTA2 <m,_a||ef>bbbb, <m,_a||ef>abab
call ccsort_grc0(3,2,2,4,4,0,syma,post,map3)
call deflength(map3,lenefbbbb)
call dawrtmap(luna2,map3,rc)
call ccsort_grc0(3,0,1,3,4,0,syma,post,map3)
call deflength(map3,lenefabab)
call dawrtmap(luna2,map3,rc)

!3   to INTA3 <m,_a||ej>aaaa, <m,_a||ej>baab, <m,_a||ej>baba
call ccsort_grc0(3,0,1,3,1,0,syma,post,map3)
call deflength(map3,lenejaaaa)
call dawrtmap(luna3,map3,rc)
call ccsort_grc0(3,0,2,3,2,0,syma,post,map3)
call deflength(map3,lenejbaab)
call dawrtmap(luna3,map3,rc)
call ccsort_grc0(3,0,2,4,1,0,syma,post,map3)
call deflength(map3,lenejbaba)
call dawrtmap(luna3,map3,rc)

!4   to INTA4 <m,_a||ej>bbbb, <m,_a||ej>abba, <m,_a||ej>abab
call ccsort_grc0(3,0,2,4,2,0,syma,post,map3)
call deflength(map3,lenejbbbb)
call dawrtmap(luna4,map3,rc)
call ccsort_grc0(3,0,1,4,1,0,syma,post,map3)
call deflength(map3,lenejabba)
call dawrtmap(luna4,map3,rc)
call ccsort_grc0(3,0,1,3,2,0,syma,post,map3)
call deflength(map3,lenejabab)
call dawrtmap(luna4,map3,rc)

! cycle over a

do a=1,nvb(syma)

  ! reconstruct #2 <_a,m,p,q> for given _a
  call mkampq(wrk,wrksize,a,ammap)

  ! get contributions to INTA2 <m,_a||ef>bbbb, <m,_a||ef>abab
  ! and write it there

  if (lenefbbbb > 0) then
    call expmpq(wrk,wrksize,syma,2,2,4,4,1,1)
    call dawri(luna2,lenefbbbb,wrk(map3%d(1,1)))
  end if

  if (lenefabab > 0) then
    call expmpq(wrk,wrksize,syma,0,1,3,4,1,0)
    call dawri(luna2,lenefabab,wrk(map3%d(1,1)))
  end if

  ! get contributions to INTA4 <m,_a||ej>bbbb, <m,_a||ej>abba, <m,_a||ej>abab
  ! and write it there

  if (lenejbbbb > 0) then
    call expmpq(wrk,wrksize,syma,0,2,4,2,1,1)
    call dawri(luna4,lenejbbbb,wrk(map3%d(1,1)))
  end if

  if (lenejabba > 0) then
    call expmpq(wrk,wrksize,syma,0,1,4,1,0,1)
    call dawri(luna4,lenejabba,wrk(map3%d(1,1)))
  end if

  if (lenejabab > 0) then
    call expmpq(wrk,wrksize,syma,0,1,3,2,1,0)
    call dawri(luna4,lenejabab,wrk(map3%d(1,1)))
  end if

  if (a > (nvb(syma)-nva(syma))) then
    ! contributions to INTA1 and INTA3 only for a-alfa

    ! get contributions to INTA1 <m,_a||ef>aaaa, <m,_a||ef>baab if any
    ! and write it there

    if (lenefaaaa > 0) then
      call expmpq(wrk,wrksize,syma,2,1,3,3,1,1)
      call dawri(luna1,lenefaaaa,wrk(map3%d(1,1)))
    end if

    if (lenefbaab > 0) then
      call expmpq(wrk,wrksize,syma,0,2,3,4,0,1)
      call dawri(luna1,lenefbaab,wrk(map3%d(1,1)))
    end if

    ! get contributions to INTA3 <m,_a||ej>aaaa, <m,_a||ej>baab, <m,_a||ej>baba
    ! and write it there

    if (lenejaaaa > 0) then
      call expmpq(wrk,wrksize,syma,0,1,3,1,1,1)
      call dawri(luna3,lenejaaaa,wrk(map3%d(1,1)))
    end if

    if (lenejbaab > 0) then
      call expmpq(wrk,wrksize,syma,0,2,3,2,0,1)
      call dawri(luna3,lenejbaab,wrk(map3%d(1,1)))
    end if

    if (lenejbaba > 0) then
      call expmpq(wrk,wrksize,syma,0,2,4,1,1,0)
      call dawri(luna3,lenejbaba,wrk(map3%d(1,1)))
    end if

  end if

end do

return

end subroutine addinta
