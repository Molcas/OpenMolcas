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

!> @brief transforms 1 particle integrals from AO to MO
!>
!> @author Matthew R. Hennefarth
!>
!> @param[in] cmo MO coefficients
!> @param[in] d_ao unscaled, folded 1 particle integral matrix in AO basis
!> @param[in] nSym number of symmetry groups
!> @param[in] nBas number of basis functions per symmetry group
!> @param[in] nOrb number of orbitals per symmetry group
!> @param[in] nFro number of frozen orbitals per symmetry group
!>
!> @param[out] d_mo unscaled, folded 1 particle integral matrix in MO basis
subroutine ao2mo_1e(cmo,d_ao,d_mo,nSym,nBas,nOrb,nFro)
  use definitions,only:iwp,wp
  use constants,only:one,zero
  use stdalloc,only:mma_allocate,mma_deallocate
  implicit none

  integer(kind=iwp),intent(in) :: nsym, nbas(*),norb(*),nfro(*)
  real(kind=wp),intent(in) :: cmo(*),d_ao(*)
  real(kind=wp),intent(out) :: d_mo(*)

  integer(kind=iwp) :: ioff1,ioff2,ioff3,isym,ibas,iorb,ifro
  real(kind=wp),allocatable :: tmp1(:),tmp2(:)

  iOff1 = 1
  iOff2 = 1
  iOff3 = 1
  Do iSym = 1,nSym
    iBas = nBas(iSym)
    iOrb = nOrb(iSym)
    If(iBas == 0 .or. iOrb == 0) then
      Cycle
    endif
    iFro = nFro(iSym)
    Call mma_allocate(Tmp1,iBas*iBas,Label='Tmp1')
    Call mma_allocate(Tmp2,iOrb*iBas,Label='Tmp2')
    Call Square(d_ao(iOff1),Tmp1,1,iBas,iBas)
    Call DGEMM_('N','N',iBas,iOrb,iBas, &
                one,Tmp1,iBas, &
                CMO(iOff2+(iFro*iBas)),iBas, &
                zero,Tmp2,iBas)
    Call DGEMM_Tri('T','N',iOrb,iOrb,iBas, &
                   one,Tmp2,iBas, &
                   CMO(iOff2+(iFro*iBas)),iBas, &
                   zero,d_mo(iOff3),iOrb)
    Call mma_deallocate(Tmp2)
    Call mma_deallocate(Tmp1)
    iOff1 = iOff1+(iBas*iBas+iBas)/2
    iOff2 = iOff2+iBas**2
    iOff3 = iOff3+(iOrb*iOrb+iOrb)/2
  EndDo
endsubroutine
