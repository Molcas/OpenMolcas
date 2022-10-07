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

subroutine defcommon(no,nv)
! this routine does:
!
! define commons needed in DIRCC routines

use ChT3_global, only: IOPT, IT, NNOAB, NNUAB, NOAB, NUAB
use Index_Functions, only: nTri_Elem
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: no, nv

! ----  UHF -----

noab(1) = no
noab(2) = no

nnoab(1) = nTri_Elem(noab(1)-1)
nnoab(2) = nTri_Elem(noab(2)-1)
nnoab(3) = noab(1)*noab(2)

nuab(1) = nv
nuab(2) = nv

nnuab(1) = nTri_Elem(nuab(1)-1)
nnuab(2) = nTri_Elem(nuab(2)-1)
nnuab(3) = nuab(1)*nuab(2)

! ----  PARAM ----

!??????????????????
it = 1

! ------ my_mpi_world_com --------

! zatial pre sekvencny chod

!mp !me = 0
!mp !nprocs = 1
!mp !llmpi = .false.
!mp !nws = 1
!mp !iws(1) = 1
!mp !lws(1) = .true.

! ------ pre Get3DM -----------

! ------    IOPT    -----------

IOPT(1) = 0

! sets IOPT(2) to an extreme number to force one file - temporary !!!  preskumat !!!
! in order to be compatible with RHF i=2^31-1
iopt(2) = 2147483647

! zatial tolko

return

end subroutine defcommon
