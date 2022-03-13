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

! Here we construct the density matrix given the orbital coefficients.
subroutine DENSI_MO(DENS,ORBCO,IS,IA,NBAS,IDIM)

implicit real*8(A-H,O-Z)
dimension DENS(*), ORBCO(IDIM,*)

IJ = 0
do I=1,NBAS
  do J=1,I
    IJ = IJ+1
    DENS(IJ) = 0.0d0
  end do
end do
do I=IS,IS+IA-1
  IJ = 0
  do J=1,NBAS
    do K=1,J
      IJ = IJ+1
      DENS(IJ) = DENS(IJ)+4.d0*ORBCO(J,I)*ORBCO(K,I)
    end do
    DENS(IJ) = DENS(IJ)-ORBCO(J,I)*ORBCO(J,I)*2.d0
  end do
end do

return

end subroutine DENSI_MO
