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
subroutine NGETH1(H,ISM,ITP,JSM,JTP)
! One-electron integrals over orbitals belonging to given OS class

use MCLR_Data, only: NTSOB
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: H(*)
integer(kind=iwp), intent(in) :: ISM, ITP, JSM, JTP
integer(kind=iwp) :: I, IJ, J, NI, NJ
real(kind=wp), external :: GTH1ES_MCLR

NI = NTSOB(ITP,ISM)
NJ = NTSOB(JTP,JSM)
IJ = 0
do J=1,NJ
  do I=1,NI
    IJ = IJ+1
    H(IJ) = GTH1ES_MCLR(I,ITP,ISM,J,JTP,JSM)
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' H1 for itp ism jtp jsm ',ITP,ISM,JTP,JSM
call WRTMAT(H,NI,NJ,NI,NJ)
#endif

end subroutine NGETH1
