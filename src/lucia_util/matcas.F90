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

subroutine MATCAS(CIN,COUT,NROWI,NROWO,IROWO1,NGCOL,ISCA,SCASGN)
! COUT(IR+IROWO1-1,ISCA(IC)) =
! COUT(IR+IROWO1-1,ISCA(IC)) + CIN(IR,IC)*SCASGN(IC)
! (if IGAT(IC) /= 0)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NROWI, NROWO, IROWO1, NGCOL, ISCA(*)
real(kind=wp) :: CIN(NROWI,*), COUT(NROWO,*), SCASGN(*)
integer(kind=iwp) :: IC, ICEXP, IR, MAXCOL, NTEST
real(kind=wp) :: SGN

MAXCOL = 0
do IC=1,NGCOL
  if (ISCA(IC) /= 0) then
    ICEXP = ISCA(IC)
    MAXCOL = max(MAXCOL,ICEXP)
    SGN = SCASGN(IC)
    do IR=1,NROWI
      COUT(IR+IROWO1-1,ICEXP) = COUT(IR+IROWO1-1,ICEXP)+SGN*CIN(IR,IC)
    end do
  end if
end do

NTEST = 0
if (NTEST /= 0) then
  write(u6,*) ' Output from MATCAS'
  call WRTMAT(COUT,NROWO,MAXCOL,NROWO,MAXCOL)
end if

end subroutine MATCAS
