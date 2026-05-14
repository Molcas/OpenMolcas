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

!#define _DEBUGPRINT_
subroutine MATCG(CIN,COUT,NROWI,NROWO,IROWI1,NGCOL,IGAT,GATSGN)
! Gather columns of CIN with phase
!
! COUT(IR,IC) = GATSGN(IC)*CIN(IR+IROWI1-1,IGAT(IC)) if IGAT(IC) /= 0
! COUT(IR,IC) = 0                                    if IGAT(IC) /= 0

use Constants, only: Zero
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NROWI, NROWO, IROWI1, NGCOL, IGAT(NGCOL)
real(kind=wp), intent(in) :: CIN(NROWI,*), GATSGN(NGCOL)
real(kind=wp), intent(out) :: COUT(NROWO,NGCOL)
integer(kind=iwp) :: IG, IGFRM
real(kind=wp) :: SGN

!write(u6,*) ' MATCG NROWI,NROWO,IROWI1,NGCOL'
!write(u6,*) NROWI,NROWO,IROWI1,NGCOL
do IG=1,NGCOL
  !write(u6,*) ' igat,sign ',IGAT(IG),GATSGN(IG)
  if (IGAT(IG) == 0) then
    COUT(:,IG) = Zero
  else
    IGFRM = IGAT(IG)
    SGN = GATSGN(IG)
    COUT(:,IG) = SGN*CIN(IROWI1:IROWI1+NROWO-1,IGFRM)
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Column gathered matrix'
call WRTMAT(COUT,NROWO,NGCOL,NROWO,NGCOL)
#endif

end subroutine MATCG
