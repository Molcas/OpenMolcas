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

use general_data, only: MALTER, NALTER, NBAS
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CMO(*)
integer(kind=iwp) :: iAlter, iAlteri, iAlterj, iOrb, iSym, jOrb
real(kind=wp), allocatable :: CMOex(:)

write(u6,*)
write(u6,'(6X,A)') 'Molecular orbitals exchanged:'

call mma_allocate(CMOex,maxval(nBas(:)),Label='CMOex')

do iAlter=1,NAlter
  iSym = MAlter(iAlter,1)
  iOrb = MAlter(iAlter,2)
  jOrb = MAlter(iAlter,3)
  write(u6,'(8X,A,I2,A,2I4)') 'In symmetry ',iSym,' :',iOrb,jOrb
  iAlterI = sum(nBas(1:iSym-1)**2)+nBas(iSym)*(iOrb-1)
  iAlterJ = sum(nBas(1:iSym-1)**2)+nBas(iSym)*(jOrb-1)
  CMOex(1:nBas(iSym)) = CMO(iAlterI+1:iAlterI+nBas(iSym))
  CMO(iAlterI+1:iAlterI+nBas(iSym)) = CMO(iAlterJ+1:iAlterJ+nBas(iSym))
  CMO(iAlterJ+1:iAlterJ+nBas(iSym)) = CMOex(1:nBas(iSym))
end do
write(u6,*)

call mma_deallocate(CMOex)

end subroutine Alter_MO
