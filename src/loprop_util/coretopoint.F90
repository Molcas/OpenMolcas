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

subroutine CoreToPoint(nAt,MP,TP)

implicit real*8(a-h,o-z)
#include "WrkSpc.fh"
real*8 MP(*), TP(*)
dimension ByggareBas(6)
logical GoHere
data ByggareBas/2.0d0,8.0d0,8.0d0,18.0d0,18.0d0,32.0d0/

! A crude but to the point algorithm to put core electrons and
! nuclei together, and separating the presumably more diffuse
! part of the charge distribution in a separate chunk.

call GetMem('NucC','Allo','Real',iNucCh,nAt)
call Get_dArray('Nuclear charge',Work(iNucCh),nAt)
kAt = 0
dToPoint = 0.0d0
do iAt=1,nAt
  ByggMeraHus = 0.0d0
  GoHere = .true.
  dScaleOffSave = Work(iNucCh+iAt-1)
  do i=1,6
    dScaleOff = dScaleOffSave-ByggareBas(i)
    if ((dScaleOff <= 0.0d0) .and. GoHere) then
      dToPoint = ByggMeraHus
      GoHere = .false.
    end if
    dScaleOffSave = dScaleOff
    ByggMeraHus = ByggMeraHus+ByggareBas(i)
  end do
  kAt = kAt+iAt
  MP(kAt) = MP(kAt)+dToPoint
  TP(iAt) = Work(iNucCh+iAt-1)-dToPoint
end do
call GetMem('NucC','Free','Real',iNucCh,nAt)

return

end subroutine CoreToPoint
