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

function memtra(npam)
! -------------------------------------------------------------------
! The following function returns the size needed by ptrans for
! temporaries. It also puts into common some offsets and stuff.
! -------------------------------------------------------------------

use etwas, only: mIrrep, nAsh, nCred, nScr1, nScr2
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: memtra
integer(kind=iwp), intent(in) :: nPam(4,0:7)
integer(kind=iwp) :: isym, mxa2, mxa3, mxa4, mxact, mxS, mxS1, mxS2, mxS234, mxS3, mxS34, mxS4, na, nscr3, nscr4, nscr5

!iQ = 1
mxact = 0
mxS1 = 0
mxS2 = 0
mxS3 = 0
mxS4 = 0
do isym=0,mirrep-1
  na = nash(isym)
  if (na == 0) cycle
  mxact = max(mxact,na)
  mxS1 = max(mxS1,npam(1,isym))
  mxS2 = max(mxS2,npam(2,isym))
  mxS3 = max(mxS3,npam(3,isym))
  mxS4 = max(mxS4,npam(4,isym))
end do
mxS = max(mxS1,mxS2,mxS3,mxS4)
mxa2 = mxact*mxact
mxa3 = mxa2*mxact
mxa4 = mxa3*mxact
mxS34 = mxS3*mxS4
mxS234 = mxS2*mxS34
! Max sizes, in common, needed for certain temporaries:
ncred = max(1,mxS*mxact)

nscr1 = mxa4
nscr2 = mxa3*mxS4
nscr3 = mxa2*mxS34
nscr4 = mxact*mxS234
nscr5 = mxS1*mxS234
nScr1 = max(1,nscr1,nscr3,nscr5)
nScr2 = max(1,nscr2,nscr4)
memtra = nCred+2*nScr1+nScr2+3

return

end function memtra
