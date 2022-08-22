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

subroutine RToSph(F,nBeta,ishll,lb,iAng,nveccb)
!***********************************************************************
!
!   Transform <core|B> from Cartesian components to Sperical harmonics
!
!***********************************************************************
! @parameter F The cartesian components of <core|B>(in)
!              The spherical components of <core|B>(out)
! @parameter nBeta Number of exponents
! @parameter ishll Shell number for ECP
! @parameter lb angular momenta Ket
! @parameter iAng angular momenta core
! @parameter Number of derivatives
!***********************************************************************

use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: Shells
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: F(*)
integer(kind=iwp) :: nBeta, ishll, lb, iAng, nveccb
integer(kind=iwp) :: ncb, nExpi
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)
! Statement function
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

!***********************************************************************

ncb = nelem(lb)*nelem(iang)
nExpi = Shells(iShll)%nExp
call mma_allocate(TMP1,nExpi*ncb*nVecCB*nBeta,Label='Tmp1')
call mma_allocate(TMP2,nExpi*ncb*nVecCB*nBeta,Label='Tmp2')

! ) kj,cbx -> cbx,kj

call DgeTMo(F,nBeta*nExpi,nBeta*nExpi,ncb*nVecCB,Tmp1,ncb*nVecCB)

! 2) bxkj,C = c,bxkj * c,C

call DGEMM_('T','N',nElem(lb)*nVecCB*nExpi*nBeta,(2*iAng+1),nElem(iAng),One,Tmp1,nElem(iAng),RSph(ipSph(iAng)),nElem(iAng),Zero, &
            Tmp2,nElem(lb)*nVecCB*nExpi*nBeta)

! 3) bx,kjC -> kjC,bx

call DgeTMo(Tmp2,nElem(lb)*nVecCB,nElem(lb)*nVecCB,nExpi*nBeta*(2*iAng+1),Tmp1,nExpi*nBeta*(2*iAng+1))

call dcopy_(nExpi*nBeta*(2*iAng+1)*nElem(lb)*nVecCB,Tmp1,1,F,1)

call mma_deallocate(Tmp1)
call mma_deallocate(Tmp2)

return

end subroutine RToSph
