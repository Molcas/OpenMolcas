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

function PiAt(MxBond,IAt,IAn,NBond,IBond)

use Definitions, only: iwp

implicit none
logical(kind=iwp) :: PiAt
integer(kind=iwp), intent(in) :: MxBond, IAt, IAn(*), NBond(*), IBond(MxBond,*)
integer(kind=iwp) :: IAnJ, ISPI, ITpJ, J, JAt, K, KAt, NBnJ, NPiJ
integer(kind=iwp), external :: IColAt

PiAt = .false.
ISPI = -1
do J=1,NBond(IAT)
  JAt = IBond(J,IAt)
  IAnJ = IAn(JAt)
  ITpJ = IColAt(IAnJ)
  NBnJ = NBond(JAt)
  NPiJ = 0
  do K=1,NBnJ
    KAt = IBond(K,JAt)
    if ((IAn(KAt) == 6) .and. (NBond(KAt) == 3)) NPiJ = NPiJ+1
  end do
  if ((IAnJ == 6) .and. (NBnJ == 3)) then
    if (NPiJ >= 2) then
      ISPI = ISPI+2
    else
      ISPI = ISPI+1
    end if
  end if
  if ((ITpJ == 5) .and. (NBnJ == 2)) ISPI = ISPI+1
  if ((ITpJ == 5) .and. (NPiJ >= 2)) ISPI = ISPI+1
end do
if (ISPI > 0) PiAt = .true.

return

end function PiAt
