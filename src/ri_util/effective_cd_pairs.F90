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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Effective_CD_Pairs(ij2,nij_Eff)

use Index_Functions, only: nTri_Elem
use Basis_Info, only: dbsc, nBas, nBas_Aux, nCnttp, Shells
use Symmetry_Info, only: nIrrep
use Cholesky, only: iSOShl, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: nij_Eff
integer(kind=iwp), allocatable, intent(out) :: ij2(:,:)
integer(kind=iwp) :: i, iAng, iCnttp, iIrrep, ij, ij_Eff, iOff, iShll, iSym, j, jOff, nAux_Tot, nij, nSkal_Valence, nVal_Tot
integer(kind=iwp), allocatable :: ij3(:), SO_ab(:)

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the max number of auxiliary shells in case of CD.
! Hence we do not have any explicit auxiliary basis set!

nSkal_Valence = 0
do iCnttp=1,nCnttp
  if (.not. dbsc(iCnttp)%Aux) then
    do iAng=0,dbsc(iCnttp)%nVal-1
      iShll = dbsc(iCnttp)%iVal+iAng
      if (.not. Shells(iShll)%Aux) nSkal_Valence = nSkal_Valence+dbsc(iCnttp)%nCntr
    end do
  end if
end do

! Max number of pairs

nij = nTri_Elem(nSkal_Valence)
call mma_allocate(ij3,nij,Label='ij3')
ij3(:) = 0
!write(u6,*) 'nij3=',nij
!                                                                      *
!***********************************************************************
!                                                                      *
nAux_Tot = 0
nVal_Tot = 0
do iIrrep=0,nIrrep-1
  nAux_Tot = nAux_Tot+nBas_Aux(iIrrep)
  nVal_Tot = nVal_Tot+nBas(iIrrep)
end do

call mma_allocate(SO_ab,2*nAux_Tot,Label='SO_ab')
SO_ab(:) = 0

iOff = 1
jOff = 0
nSym = nIrrep
do iSym=1,nSym
  iIrrep = iSym-1
  call CHO_X_GET_PARDIAG(iSym,SO_ab(iOff))

  call Get_Auxiliary_Shells(SO_ab(iOff),nBas_Aux(iIrrep),jOff,iSOShl,nVal_Tot,ij3,nij)

  jOff = jOff+nBas_Aux(iIrrep)
  iOff = iOff+2*nBas_Aux(iIrrep)
end do
call mma_deallocate(SO_ab)
!                                                                      *
!***********************************************************************
!                                                                      *
nij_Eff = 0
do i=1,nij
  nij_Eff = nij_Eff+ij3(i)
end do
!write(u6,*) 'nij_Eff=',nij_Eff
if (nij_Eff > nij) then
  call WarningMessage(2,'Effective_CD_Pairs: nij_Eff > nij')
  call Abend()
end if

call mma_allocate(ij2,2,nij_Eff,Label='ij2')
ij = 0
ij_Eff = 0
do i=1,nSkal_Valence
  do j=1,i
    ij = ij+1
    !write (u6,*) 'i,j,ij=',i,j,ij
    if (ij3(ij) == 1) then
      ij_Eff = ij_Eff+1
      !write(u6,*) 'ij_Eff=',ij_Eff
      ij2(1,ij_Eff) = i
      ij2(2,ij_Eff) = j
    end if
  end do
end do
if (ij_Eff /= nij_Eff) then
  call WarningMessage(2,'Effective_CD_Pairs: ij_Eff /= nij_Eff')
  call Abend()
end if
call mma_deallocate(ij3)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Effective_CD_Pairs
