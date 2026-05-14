!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine SOS(iStabO,nStabO,lOper)
!***********************************************************************
!                                                                      *
! Object: to generate the stabilizer S for the operator O.             *
!                                                                      *
!     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN                               *
!             February '91                                             *
!***********************************************************************

use Symmetry_Info, only: iChTbl, iOper, nIrrep
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: iStabO(8), nStabO
integer(kind=iwp), intent(in) :: lOper
integer(kind=iwp) :: iIrrep, iS

#ifdef _DEBUGPRINT_
write(u6,*) ' In SOS'
write(u6,*) ' lOper=',lOper
do iS=0,nIrrep-1
  write(u6,'(8I5)') (iChTbl(iIrrep,iS),iIrrep=0,nIrrep-1)
end do
#endif
if ((lOper < 0) .or. (lOper > 255)) then
  call WarningMessage(2,'SOS: Symmetry label is corrupted.')
  write(u6,*) 'lOper=',lOper
  call Abend()
end if
nStabO = 0
outer: do iS=0,nIrrep-1
  inner: do iIrrep=0,nIrrep-1
    if (.not. btest(lOper,iIrrep)) cycle inner
    if (iChTbl(iIrrep,iS) /= 1) cycle outer
  end do inner
  nStabO = nStabO+1
  iStabO(nStabO) = iOper(iS)
end do outer

end subroutine SOS
