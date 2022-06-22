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

subroutine Get_OCOF(nPrim,nBas,Work,nVec_p,OCOF)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nPrim, nBas, nVec_p
real(kind=wp), intent(in) :: Work(nVec_p)
real(kind=wp), intent(out) :: oCof(nBas,nPrim)
integer(kind=iwp) :: i, iStdOut, iVec_p, j

iStdOut = u6 ! Added by EB
iVec_p = 0
do i=1,nBas
  do j=1,nPrim
    OCOF(i,j) = Work((i-1)*nPrim+j)
    iVec_p = iVec_p+1
    if (iVec_p > nVec_p) then
      write(iStdOut,*) 'iVec_p > nVec_p'
      write(iStdOut,*) iVec_p,' > ',nVec_p
      write(iStdOut,*) 'nPrim=',nPrim
      call Abend()
    end if
  end do
end do

return

end subroutine Get_OCOF
