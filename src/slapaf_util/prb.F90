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

subroutine PrB(dB,nq,nDim)

implicit real*8(a-h,o-z)
real*8 dB(nq,nDim,nDim)

do iq=1,nQ
  write(6,*) ' iq=',iq
  do iDim=1,nDim
    write(6,'(9F10.6)') (dB(iq,iDim,jDim),jDim=1,nDim)
  end do
end do

return

end subroutine PrB
