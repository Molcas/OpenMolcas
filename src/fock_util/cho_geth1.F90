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

subroutine CHO_GETH1(nBtri,H1,RFpert,ERFNuc)

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
integer nBTri
real*8 H1(nBTri)
logical RFpert
real*8 ERFNuc
character*8 OneLbl
real*8, allocatable :: Tmp(:)

iRc = -1
iOpt = 6
iCmp = 1
iSyLab = 1
OneLbl = 'OneHam  '

call RdOne(iRc,iOpt,OneLbl,iCmp,H1,iSyLab)

if (IRC /= 0) then
  write(6,*)
  write(6,*) '    *** ERROR IN SUBROUTINE  CHO_GETH1 *** '
  write(6,*) '   BARE NUCLEI HAMILTONIAN IS NOT AVAILABLE'
  write(6,*)
  call Abend()
end if

ERFNuc = Zero
if (RFpert) then
  call mma_allocate(Tmp,nBTri,Label='Tmp')
  call Get_dScalar('RF Self Energy',ERFNuc)
  call Get_dArray('Reaction field',Tmp,nBtri)

  H1(:) = H1(:)+Tmp(:)

  call mma_deallocate(Tmp)
end if

return

end subroutine CHO_GETH1
