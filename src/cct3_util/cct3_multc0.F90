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

subroutine cct3_multc0(wrk,wrksize,mvec,ix,c,key)
! This routine realizes multiplying according mvec
! for C=A*B
! N.B. if key=0, C file is not vanished (ie can be used for
! adding to some existing file)
!
! If C=A*B process is faster or comparable with C=AT*B then mchntyp should be set to 1.
! If C=AT*B is significantly faster than C=A*T (more than 20%), then mchntyp should be set to 2. (default is 1)
! if mchntyp is 2, then
!   1) processes with scale(A)/scale(B) > scalelim will be calculated as C=A*B
!   2) processes with scale(A)/scale(B) < scalelim will be calculated as C=AT*B
! Note, that for mchntyp =2 more memory is required, due to requirement of
! aditional o2v2 help file posd0

use CCT3_global, only: Map_Type, mchntyp, posd0, slim
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: wrksize, mvec(4096,7), ix, key
real(kind=wp), intent(inout) :: wrk(wrksize)
type(Map_Type), intent(in) :: c
integer(kind=iwp) :: ic, iix, nhelp1, nhelp2, nhelp3, nhelp4, nhelp5, nhelp6
real(kind=wp) :: sc

!1 set C=0

if (key == 1) then

  ! C matrix must be vanished

  do ic=1,c%d(0,5)
    nhelp1 = c%d(ic,1)
    nhelp2 = c%d(ic,2)
    wrk(nhelp1:nhelp1+nhelp2-1) = Zero
  end do

end if

!2 C=C+A*B

do iix=1,ix

  ! skip this summation if yes/no=0
  if (mvec(iix,1) == 0) cycle

  ! realize individual summation

  ! def positions of A,B,C
  nhelp1 = mvec(iix,2)
  nhelp2 = mvec(iix,3)
  nhelp3 = mvec(iix,4)

  ! def rowA(rowC), colA(rowB,sum), colB(colC)
  nhelp4 = mvec(iix,5)
  nhelp5 = mvec(iix,6)
  nhelp6 = mvec(iix,7)

  if (mchntyp == 1) then

    ! Typ 1
    call cct3_mc0c1a3b(nhelp4,nhelp5,nhelp5,nhelp6,nhelp4,nhelp6,nhelp4,nhelp5,nhelp6,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))

  else

    ! Typ2
    sc = real(nhelp4,kind=wp)/real(nhelp6,kind=wp)
    if (sc > slim) then
      call cct3_mc0c1a3b(nhelp4,nhelp5,nhelp5,nhelp6,nhelp4,nhelp6,nhelp4,nhelp5,nhelp6,wrk(nhelp1),wrk(nhelp2),wrk(nhelp3))

    else
      ! map D=AT
      call cct3_map21(wrk(nhelp1),wrk(posd0),nhelp4,nhelp5,2,1,1)
      ! calc C=DT*B
      call cct3_mc0c1at3b(nhelp5,nhelp4,nhelp5,nhelp6,nhelp4,nhelp6,nhelp4,nhelp5,nhelp6,wrk(posd0),wrk(nhelp2),wrk(nhelp3))
    end if

  end if

end do

return

end subroutine cct3_multc0
