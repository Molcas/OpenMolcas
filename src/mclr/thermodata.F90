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

subroutine ThermoData(in_Freq,in_nFreq)
! Compute thermodynamic data at different temperatures.

use Constants, only: auTocm
use Temperatures, only: DefTemp

implicit real*8(a-h,o-z)
#include "Molcas.fh"
real*8 in_Freq(in_nFreq), Freq(MxAtom*3-6)

! Remove translational and rotational frequencies

nFreq = 0
do i=1,in_nFreq
  if (in_Freq(i) > 20.0d0) then
    nFreq = nFreq+1
    Freq(nFreq) = in_Freq(i)
  end if
end do

! Is the system linear?

nAtom = (nFreq+6)/3
nTR = 3*nAtom-nFreq ! Number of trans and rot fg

do i=1,nFreq
  ! Convert frequecnies from cm-1 to hartree
  !Freq(i) = Freq(i)*4.55633538D-06
  Freq(i) = Freq(i)/auTocm
end do

do i=1,size(DefTemp)
  call Thermo_Vib(nFreq,Freq,DefTemp(i),nTR,i)
end do

end subroutine ThermoData
