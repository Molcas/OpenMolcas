!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2003, Giovanni Ghigo                                   *
!***********************************************************************

subroutine Alter_MO(CMO)
!***********************************************************************
!                                                                      *
!    purpose:                                                          *
!    The new keyword ALTEr exchanges pair of MO taken from the files   *
!    INPORB or JOBOLD before to start the RASSCF calculation.          *
!    The keyword must be followed by the number of pairs to exchange   *
!    NAlter and, for each pair, by the symmetry specie MAlter(iAlter   *
!    ,1) and the two indices of the MO to exchange MAlter(iAlter,2-3). *
!    NALTEr and MAlter(MaxAlter,3) (the dimension is fixed) are passed *
!    through module general_data.                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     G. Ghigo                                                         *
!     University of Lund, Sweden, September 2003                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use output_ras, only: LF
use general_data, only: NALTER, MALTER, NBAS

implicit none
real*8 CMO(*)
real*8 CMOex
integer iAlter, iAlteri, iAlterj, iCMO, iSym

! Local print level (if any)
write(LF,*)
write(LF,'(6X,A)') 'Molecular orbitals exchanged:'

do iAlter=1,NAlter
  write(LF,'(8X,A,I2,A,2I4)') 'In symmetry ',MAlter(iAlter,1),' :',MAlter(iAlter,2),MAlter(iAlter,3)
  iAlterI = 0
  iAlterJ = 0
  if (MAlter(iAlter,1) > 1) then
    do iSym=1,MAlter(iAlter,1)-1
      iAlterI = iAlterI+nBas(iSym)**2
      iAlterJ = iAlterJ+nBas(iSym)**2
    end do
  end if
  iAlterI = iAlterI+nBas(MAlter(iAlter,1))*(MAlter(iAlter,2)-1)
  iAlterJ = iAlterJ+nBas(MAlter(iAlter,1))*(MAlter(iAlter,3)-1)
  do iCMO=1,nBas(MAlter(iAlter,1))
    CMOex = CMO(iAlterI+iCMO)
    CMO(iAlterI+iCMO) = CMO(iAlterJ+iCMO)
    CMO(iAlterJ+iCMO) = CMOex
  end do
end do
write(LF,*)

end subroutine Alter_MO
