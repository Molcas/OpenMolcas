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

subroutine Compute_V12_Inner(V,V12,VTri,Vec,nDim)

implicit real*8(A-H,O-Z)
#include "real.fh"
real*8 V(nDim,nDim), V12(nDim,nDim), VTri(nDim*(nDim+1)/2), Vec(nDim,nDim)

call FZero(Vec,nDim**2)
call dcopy_(nDim,[One],0,Vec,nDim+1)

do i=1,nDim
  do j=1,i
    VTri(i*(i-1)/2+j) = V(i,j)
  end do
end do

call JACOB(VTri,Vec,nDim,nDim)

call FZero(V12,nDim**2)
do i=1,nDim
# ifdef _DEBUGPRINT_
  tmp = VTri(i*(i+1)/2)
  write(6,*) 'i,tmp=',i,tmp
# endif
  tmp = sqrt(VTri(i*(i+1)/2))
  if (tmp < 1.0D-90) then
    V12(i,i) = 1.0d90
  else
    V12(i,i) = One/tmp
  end if
end do

call DGEMM_('N','T',nDim,nDim,nDim,1.0d0,V12,nDim,Vec,nDim,0.0d0,V,nDim)
call DGEMM_('N','N',nDim,nDim,nDim,1.0d0,Vec,nDim,V,nDim,0.0d0,V12,nDim)

return

end subroutine Compute_V12_Inner
