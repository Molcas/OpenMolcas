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
subroutine MATCAS(CIN,COUT,NROWI,NROWO,IROWO1,NGCOL,ISCA,SCASGN)
! COUT(IR+IROWO1-1,ISCA(IC)) =
! COUT(IR+IROWO1-1,ISCA(IC)) + CIN(IR,IC)*SCASGN(IC)
! (if IGAT(IC) /= 0)

use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NROWI, NROWO, IROWO1, NGCOL, ISCA(NGCOL)
real(kind=wp), intent(in) :: CIN(NROWI,*), SCASGN(*)
real(kind=wp), intent(inout) :: COUT(NROWO,*)
integer(kind=iwp) :: IC, ICEXP, MAXCOL
real(kind=wp) :: SGN

MAXCOL = 0
do IC=1,NGCOL
  if (ISCA(IC) /= 0) then
    ICEXP = ISCA(IC)
    MAXCOL = max(MAXCOL,ICEXP)
    SGN = SCASGN(IC)
    COUT(IROWO1:IROWO1+NROWI-1,ICEXP) = COUT(IROWO1:IROWO1+NROWI-1,ICEXP)+SGN*CIN(1:NROWI,IC)
  end if
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Output from MATCAS'
call WRTMAT(COUT,NROWO,MAXCOL,NROWO,MAXCOL)
#endif

end subroutine MATCAS
