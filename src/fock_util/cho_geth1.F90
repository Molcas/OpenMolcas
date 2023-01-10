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

subroutine CHO_GETH1(nBTri,H1,RFpert,ERFNuc)

use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBTri
real(kind=wp), intent(out) :: H1(nBTri), ERFNuc
logical(kind=iwp), intent(in) :: RFpert
integer(kind=iwp) :: iCmp, iOpt, iRc, iSyLab
character(len=8) :: OneLbl
real(kind=wp), allocatable :: Tmp(:)

iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iCmp = 1
iSyLab = 1
OneLbl = 'OneHam  '

call RdOne(iRc,iOpt,OneLbl,iCmp,H1,iSyLab)

if (IRC /= 0) then
  write(u6,*)
  write(u6,*) '    *** ERROR IN SUBROUTINE  CHO_GETH1 ***'
  write(u6,*) '   BARE NUCLEI HAMILTONIAN IS NOT AVAILABLE'
  write(u6,*)
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
