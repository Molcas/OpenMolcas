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

subroutine Flip_Flop(Primitive)

use Basis_Info, only: DBSC, iCnttp_Dummy, nCnttp, Shells
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(in) :: Primitive
integer(kind=iwp) :: iAng, iCnttp, iShll, nExpi, nTest

S%MaxBas(:) = 0
S%MaxPrm(:) = 0

do iCnttp=1,nCnttp
  nTest = dbsc(iCnttp)%nVal-1
  if (DBSC(iCnttp)%Aux .and. (iCnttp == iCnttp_Dummy)) nTest = -1

  do iAng=0,S%iAngMx
    if (iAng > nTest) cycle
    iShll = dbsc(iCnttp)%iVal+iAng
    nExpi = Shells(iShll)%nExp
    if (nExpi == 0) cycle
    if (Shells(iShll)%nBasis_C == 0) cycle

    ! Decontract only the ordinary basis sets!

    call mma_deallocate(Shells(iShll)%pCff)
    if (Primitive .and. (.not. Shells(iShll)%Aux) .and. (.not. Shells(iShll)%Frag)) then
      Shells(iShll)%nBasis = nExpi
      call mma_allocate(Shells(iShll)%pCff,nExpi,Shells(iShll)%nBasis,Label='pCff')
      Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_p(:,:,1)
    else
      Shells(iShll)%nBasis = Shells(iShll)%nBasis_C
      call mma_allocate(Shells(iShll)%pCff,nExpi,Shells(iShll)%nBasis,Label='pCff')
      Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
    end if
    S%MaxPrm(iAng) = max(S%MaxPrm(iAng),nExpi)
    S%MaxBas(iAng) = max(S%MaxBas(iAng),Shells(iShll)%nBasis)

  end do ! iAng
end do   ! iCnttp

return

end subroutine Flip_Flop
