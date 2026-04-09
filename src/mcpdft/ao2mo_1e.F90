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
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************
!  ao2mo_1e
!
!> @brief transforms 1 particle integrals from AO to MO
!>
!> @author Matthew R. Hennefarth
!>
!> @param[in] cmo  MO coefficients
!> @param[in] d_ao unscaled, folded 1 particle integral matrix in AO basis
!> @param[in] nSym number of symmetry groups
!> @param[in] nBas number of basis functions per symmetry group
!> @param[in] nOrb number of orbitals per symmetry group
!> @param[in] nFro number of frozen orbitals per symmetry group
!>
!> @param[out] d_mo unscaled, folded 1 particle integral matrix in MO basis
!***********************************************************************

subroutine ao2mo_1e(cmo,d_ao,d_mo,nSym,nBas,nOrb,nFro)

use Index_Functions, only: nTri_Elem
use general_data, only: ntot1, ntot2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: cmo(ntot2), d_ao(ntot1)
real(kind=wp), intent(out) :: d_mo(ntot1)
integer(kind=iwp), intent(in) :: nsym, nbas(nsym), norb(nsym), nfro(nsym)
integer(kind=iwp) :: ibas, ifro, ioff1, ioff2, ioff3, iorb, isym
real(kind=wp), allocatable :: tmp1(:), tmp2(:)

iOff1 = 1
iOff2 = 1
iOff3 = 1
do iSym=1,nSym
  iBas = nBas(iSym)
  iOrb = nOrb(iSym)
  if ((iBas == 0) .or. (iOrb == 0)) cycle
  iFro = nFro(iSym)
  call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
  call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
  call Square(d_ao(iOff1),Tmp1,1,iBas,iBas)
  call DGEMM_('N','N',iBas,iOrb,iBas, &
              One,Tmp1,iBas, &
              CMO(iOff2+(iFro*iBas)),iBas, &
              Zero,Tmp2,iBas)
  call DGEMM_Tri('T','N',iOrb,iOrb,iBas, &
                 One,Tmp2,iBas, &
                 CMO(iOff2+(iFro*iBas)),iBas, &
                 Zero,d_mo(iOff3),iOrb)
  call mma_deallocate(Tmp2)
  call mma_deallocate(Tmp1)
  iOff1 = iOff1+nTri_Elem(iBas)
  iOff2 = iOff2+iBas**2
  iOff3 = iOff3+nTri_Elem(iOrb)
end do

end subroutine ao2mo_1e
