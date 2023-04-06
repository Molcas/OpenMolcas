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

subroutine init(wrk,wrksize,lunabij1,lunabij2,lunabij3)
! this routine does FI1,FII1,FIII1,T11,T21
! 1) def F1(a,e) = fok(a,e)
! 2) def F2(m,i) = fok(m,i)
! 3) def F3(e,m) = fok(e,m)
! 4) def T1n(a,i) = fok(a,i)
! 5) def T2m(abij)= <ab||ij>
!
! lunabij1 - lun of file, where <ab||ij>aaaa is stored (I)
! lunabij2 - lun of file, where <ab||ij>bbbb is stored (I)
! lunabij3 - lun of file, where <ab||ij>abab is stored (I)
!
! N.B. this routine uses and destroys help files : none

use ccsd_global, only: f11, f12, f21, f22, f31, f32, fk1, fk2, fk3, fk4, fk5, fk6, t13, t14, t21, t22, t23
use Para_Info, only: MyRank
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: wrksize
real(kind=wp), intent(inout) :: wrk(wrksize)
integer(kind=iwp), intent(_IN_) :: lunabij1, lunabij2, lunabij3
integer(kind=iwp) :: post, rc

!1.1 map fok(a,b)aa to f1(a,e)aa
call map(wrk,wrksize,2,1,2,0,0,fk1,1,f11,post,rc)

!1.2 map fok(a,b)bb to f1(a,e)bb
call map(wrk,wrksize,2,1,2,0,0,fk2,1,f12,post,rc)

!2.1 map fok(i,j)aa to f2(i,j)aa
call map(wrk,wrksize,2,1,2,0,0,fk5,1,f21,post,rc)

!2.2 map fok(i,j)bb to f2(i,j)bb
call map(wrk,wrksize,2,1,2,0,0,fk6,1,f22,post,rc)

!3.1 map fok(a,i)aa to f3(a,i)aa
call map(wrk,wrksize,2,1,2,0,0,fk3,1,f31,post,rc)

!3.2 map fok(a,i)bb to f3(a,i)bb
call map(wrk,wrksize,2,1,2,0,0,fk4,1,f32,post,rc)

if (myRank == 0) then

  !4.1 map fok(a,i)aa to t1n(a,i)aa
  call map(wrk,wrksize,2,1,2,0,0,fk3,1,t13,post,rc)

  !4.2 map fok(a,i)bb to t1n(a,i)bb
  call map(wrk,wrksize,2,1,2,0,0,fk4,1,t14,post,rc)

else

  !4.3 set t1naa (t13) =0
  call set0(wrk,wrksize,t13)

  !4.4 set t1nbb (t14) =0
  call set0(wrk,wrksize,t14)

end if

if (myRank == 0) then

  !5.1 load <ab||ij>aaaa from lunabij1 to t2n(ab,ij)aaaa
  call filemanager(2,lunabij1,rc)
  call getmediate(wrk,wrksize,lunabij1,t21,rc)

  !5.2 load <ab||ij>bbbb from lunabij2 to t2n(ab,ij)bbbb
  call filemanager(2,lunabij2,rc)
  call getmediate(wrk,wrksize,lunabij2,t22,rc)

  !5.3 load <ab||ij>abab from lunabij3 to t2n(ab,ij)abab
  call filemanager(2,lunabij3,rc)
  call getmediate(wrk,wrksize,lunabij3,t23,rc)

else

  !5.4 set t2naaaa (t21) =0
  call set0(wrk,wrksize,t21)

  !5.5 set t2nbbbb (t22) =0
  call set0(wrk,wrksize,t22)

  !5.6 set t2nabab (t23) =0
  call set0(wrk,wrksize,t23)

end if

return

end subroutine init
