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
! One-electron integrals over orbitals belonging to
! given OS class

use MCLR_Data, only: NTSOB

implicit none
!.Output
real*8 H(*)
integer ISM, ITP, JSM, JTP
integer NI, NJ, IJ, I, J
real*8, external :: GTH1EN

NI = NTSOB(ITP,ISM)
NJ = NTSOB(JTP,JSM)
IJ = 0
do J=1,NJ
  do I=1,NI
    IJ = IJ+1
    H(IJ) = GTH1EN(I,ITP,ISM,J,JTP,JSM)
  end do
end do

#ifdef _DEBUGPRINT_
write(6,*) ' H1 for itp ism jtp jsm ',ITP,ISM,JTP,JSM
call WRTMAT(H,NI,NJ,NI,NJ)
#endif

end subroutine NGETH1
