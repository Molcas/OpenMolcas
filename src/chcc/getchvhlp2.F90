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

subroutine GetChVHlp2(L2Status,NL2,used,kam)
! this routine does:
! Define, which one of unused L2 arrays allocated in memory
! will be used for placing newly needed L2x
! Criteria of present approach:
! 1) First those positions that are allocated, but never used yet
! 2) if no 1), take first unused in this step, according L2Status
!    N.B. 2) druhy krok moze byt eventuelne vylepseny, zatial
!         takto odflaknute
!
! description of variables:
! L2Status  - L2 Status array (I)
! NL2       - number of L2 arrays reserved in memory (I)
! used      - array indicating, which L2x are already used
!             in this step (those cannot be touched) (I)
! kam       - index of position, where new L2 can be placed (O)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: L2Status(4,3), NL2, used(4)
integer(kind=iwp), intent(out) :: kam
integer(kind=iwp) :: i

!1 search, if there are never used positions

do i=1,NL2
  if (L2Status(i,1) == 0) then
    kam = i
    return
  end if
end do

!2 find first unused in this step

do i=1,NL2
  if (used(i) == 0) then
    kam = i
    return
  end if
end do

! Jaj nieje dobre ak sme sa dostali az sem
write(u6,*) ' Sorry fish getChVHlp2 '
call Abend()

return

end subroutine GetChVHlp2
