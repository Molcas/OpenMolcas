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

subroutine HssPrt_MCLR(ideg,Hess,ldisp)

use Index_Functions, only: iTri, nTri_Elem
use input_mclr, only: ChIrr, nIrrep, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ideg(*), ldisp(nsym)
real(kind=wp), intent(in) :: Hess(*)
integer(kind=iwp) :: i, iaa, idisp, ii, iIrrep, j, jj, kDisp(8)
character(len=39) :: Title
real(kind=wp), allocatable :: Temp(:)

iDisp = 0
do iIrrep=1,nIrrep
  kDisp(iIrrep) = iDisp
  iDisp = iDisp+lDisp(iIrrep)
  write(u6,*) lDisp(iIrrep)
end do

call mma_allocate(Temp,iDisp**2,Label='Temp')
iaa = 0
do iIrrep=1,nIrrep
  if (ldisp(iirrep) /= 0) then
    write(title,'(A,A)') 'Hessian in Irrep ',chirr(iIrrep)

    do i=1,lDisp(iirrep)
      do j=1,i
        ii = iTri(i,j)
        jj = iaa+iTri(i,j)
        Temp(ii) = Hess(jj)*sqrt(real(ideg(i+kdisp(iirrep))*ideg(j+kdisp(iirrep)),kind=wp))
      end do
    end do
    call TriPrt(title,' ',Temp,ldisp(iirrep))
    iaa = iaa+nTri_Elem(ldisp(iirrep))
  end if
end do
call mma_deallocate(Temp)

end subroutine HssPrt_MCLR
