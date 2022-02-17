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

subroutine DCORR(JREFX,AREF,ICSPCK,DMO)

use mrci_global, only: ENP, IPRINT, IRC, LN, Lu_27, NREF
use Constants, only: One
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JREFX(*), ICSPCK(*)
real(kind=wp), intent(_OUT_) :: AREF(*), DMO(*)
integer(kind=iwp) :: I, IAD27, II1, IJ, IK, INDA, IOC
real(kind=wp) :: FAC, TSUM
integer(kind=iwp), external :: ICUNP

! CORRECTION TO DENSITY MATRIX IN ACPF CASE.
if (IPRINT >= 7) write(u6,*) ' ENP IN DENS =',ENP
FAC = One-(One/ENP)
IAD27 = 0
call dDAFILE(Lu_27,2,AREF,NREF,IAD27)
IK = 0
do INDA=1,IRC(1)
  II1 = (INDA-1)*LN
  if (JREFX(INDA) /= 0) then
    IK = IK+1
    TSUM = AREF(IK)*AREF(IK)*FAC
    IJ = 0
    do I=1,LN
      IOC = (1+ICUNP(ICSPCK,II1+I))/2
      IJ = IJ+I
      DMO(IJ) = DMO(IJ)+IOC*TSUM
    end do
  end if
end do

return

end subroutine DCORR
