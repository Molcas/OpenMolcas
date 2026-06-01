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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine DEPSATrf(NBSQT,nAshT,DEPSA,FPT2,WRK1,WRK2)

use caspt2_global, only: CMOPT2
use caspt2_module, only: IfChol, NASH, NBAS, NBAST, NFRO, NISH, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBSQT, nAshT
real(kind=wp), intent(in) :: DEPSA(nAshT,nAshT)
real(kind=wp), intent(out) :: FPT2(NBSQT), WRK1(NBSQT), WRK2(NBSQT)
integer(kind=iwp) :: iAshI, iOrb, iSym, iSymA, iSymB, iSymI, iSymJ, jAsh, jAshI, jOrb, nBasI, nCorI
real(kind=wp), allocatable :: DAO(:), DMO(:)

FPT2(1:nBasT**2) = Zero

iSym = 1
iSymA = 1
iSymI = 1
iSymB = 1
iSymJ = 1

!if ((nFroT /= 0) .and. IfChol) then
if (IfChol) then
  !! DEPSA(MO) -> DEPSA(AO) -> G(D) in AO -> G(D) in MO
  !! The Cholesky vectors do not contain frozen orbitals...
  call mma_allocate(DAO,NBSQT,Label='DAO')
  call mma_allocate(DMO,NBSQT,Label='DMO')
  !! First, MO-> AO transformation of DEPSA
  do iSym=1,nSym
    DMO(:) = Zero
    nCorI = nFro(iSym)+nIsh(iSym)
    nBasI = nBas(iSym)
    do jAsh=1,nAsh(iSym)
      DMO(nCorI+nBasI*(nCorI+jAsh-1)+1:nCorI+nBasI*(nCorI+jAsh-1)+nAsh(iSym)) = DEPSA(1:nAsh(iSym),jAsh)
    end do
    call OLagTrf(1,iSym,NBSQT,CMOPT2,DMO,DAO,WRK1)
  end do
  !! Compute G(D)
  WRK1(1:NBSQT) = Zero
  DMO(:) = Zero
  !! it's very inefficient
  call OLagFro4(NBSQT,1,1,1,1,1,DAO,WRK1,DMO,WRK1,WRK2)
  !! G(D) in AO -> G(D) in MO
  do iSym=1,nSym
    call OLagTrf(2,iSym,NBSQT,CMOPT2,FPT2,DMO,WRK1)
  end do
  call mma_deallocate(DAO)
  call mma_deallocate(DMO)
else
  nCorI = nFro(iSym)+nIsh(iSym)
  do iAshI=1,nAsh(iSym)
    iOrb = nCorI+iAshI
    do jAshI=1,nAsh(iSym)
      jOrb = nCorI+jAshI

      call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
      FPT2(1:nBasT**2) = FPT2(1:nBasT**2)+DEPSA(iAshI,jAshI)*WRK1(1:nBasT**2)

      call Exch(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
      FPT2(1:nBasT**2) = FPT2(1:nBasT**2)-Half*DEPSA(iAshI,jAshI)*WRK1(1:nBasT**2)
    end do
  end do
end if

return

end subroutine DEPSATrf
