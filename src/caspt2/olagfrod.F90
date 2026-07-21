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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine OLagFroD(NBSQT,NASHT,DIA,DI,RDMSA,Trf)

use Index_Functions, only: nTri_Elem
use caspt2_global, only: CMOPT2
use general_data, only: NASH
use caspt2_module, only: NBAS, NFRO, NISH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBSQT, nAshT
real(kind=wp), intent(out) :: DIA(NBSQT), DI(NBSQT)
real(kind=wp), intent(in) :: RDMSA(nAshT**2), Trf(NBSQT)
integer(kind=iwp) :: iAOsq, iAOtr, iSym, nAshI, nBasI, nCorI, nFroI, nIshI
real(kind=wp), allocatable :: WRK1(:), WRK2(:)

call mma_allocate(WRK1,NBSQT,Label='WRK1')
call mma_allocate(WRK2,NBSQT,Label='WRK2')

iAOtr = 0
iAOsq = 1
do iSym=1,nSym
  nFroI = nFro(iSym)
  nIshI = nIsh(iSym)
  nAshI = nAsh(iSym)
  nBasI = nBas(iSym)
  nCorI = nFroI+nIshI

  !! full density matrix
  !call SQUARE(WRK1(1+iAOtr),DIA(iAOsq),1,nBasI,nBasI)
  ! !! off-diagonal elements have to be halved
  !do Mu=1,nBasI
  !  do Nu=1,nBasI
  !    if (Mu == Nu) cycle
  !    DIA(iAOsq+Mu-1+nBasI*(Nu-1)) = Half*DIA(iAOsq+Mu-1+nBasI*(Nu-1))
  !  end do
  !end do

  !! inactive density matrix
  call DGEMM_('N','T',nBasI,nBasI,nCorI,Two,CMOPT2,nBasI,CMOPT2,nBasI,Zero,DI(iAOsq),nBasI)

  !! inactive+active density matrix
  !! Somehow, the above density matrix obtained by calling
  !! Get_D1AO is incorrect... at least, cannot be used.
  ! 1) inactive part
  DIA(1:nBasI**2) = DI(1:nBasI**2)
  ! 2)  RDMSA is defined in CASSCF orbitals, so transform RDMSA to
  !     CASPT2 orbital basis
  call DGemm_('T','N',nAshI,nAshI,nAshI,One,Trf(1+nCorI+nBasI*nCorI),nBasI,RDMSA,nAshI,Zero,WRK2,nAshI)
  call DGemm_('N','N',nAshI,nAshI,nAshI,One,WRK2,nAshI,Trf(1+nCorI+nBasI*nCorI),nBasI,Zero,WRK1,nAshI)
  ! 3) Finally, add the active part
  call DGemm_('N','N',nBasI,nAshI,nAshI,One,CMOPT2(1+nBasI*nCorI),nBasI,WRK1,nAshI,Zero,WRK2,nBasI)
  call DGemm_('N','T',nBasI,nBasI,nAshI,One,WRK2,nBasI,CMOPT2(1+nBasI*nCorI),nBasI,One,DIA,nBasI)

  iAOtr = iAOtr+nTri_Elem(nBasI)
  iAOsq = iAOsq+nBasI*nBasI
end do

call mma_deallocate(WRK1)
call mma_deallocate(WRK2)

return

end subroutine OLagFroD
