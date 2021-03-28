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

implicit real*8(a-h,o-z)
dimension oCof(nBas,nPrim)
dimension Work(nVec_p)

iStdOut = 6 ! Added by EB
iVec_p = 0
do i=1,nBas
  do j=1,nPrim
    OCOF(i,j) = Work((i-1)*nPrim+j)
    iVec_p = iVec_p+1
    if (iVec_p > nVec_p) then
      write(iStdOut,*) 'iVec_p.gt.nVec_p'
      write(iStdOut,*) iVec_p,'.gt.',nVec_p
      write(iStdOut,*) 'nPrim=',nPrim
      call ABEND()
    end if
  end do
end do

return

end subroutine Get_OCOF
