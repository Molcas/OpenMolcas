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

subroutine Pos_QLast(Disc)

use definitions, only: iwp, wp, u6, RtoI
use TList_Mod, only: iTskCan, Not_Used, TskQ, QLast

implicit none
real(kind=wp), intent(inout) :: Disc
integer(kind=iwp) iWR(2), mInts
real(kind=wp) Dummy(1), Quad_ijkl, RST_triplet
logical :: Copy = .true., NoCopy = .false.

if (.not. allocated(TskQ)) return

Quad_ijkl = TskQ(1,iTskCan)
RST_triplet = TskQ(2,iTskCan)
if (Quad_ijkl == Not_Used) return

! If already at the right position return

if ((Quad_ijkl == QLast(1)) .and. (RST_triplet == QLast(2))) return

do

  call iRBuf(iWR,2,Copy)
  call dRBuf(QLast,2,Copy)
  mInts = iWR(2)
  if ((QLast(1) == Quad_ijkl) .and. (QLast(2) == RST_triplet)) then
    if (mInts > 0) call dRBuf(Dummy,mInts,NoCopy)
    Disc = Disc+dble(2/RtoI+2+mInts)
    return
  else if (QLast(1) <= Quad_ijkl) then
    if (mInts > 0) call dRBuf(Dummy,mInts,NoCopy)
    Disc = Disc+dble(2/RtoI+2+mInts)
    cycle
  else
    write(u6,*) 'Pos_QLast: batch is lost!'
    write(u6,'(A,2F10.1)') 'Index,1.0:  ',QLast(1),QLast(2)
    write(u6,'(A,2F10.1)') 'Looking for ',Quad_ijkl,RST_triplet
    write(u6,*) ' iTskCan,=',iTskCan
    call RecPrt('TskQ',' ',TskQ,2,iTskCan)
    write(u6,*)
    call Abend()
  end if

end do

write(u6,*) 'Pos_QLast: Fatal problem!'
call Abend()

end subroutine Pos_QLast
