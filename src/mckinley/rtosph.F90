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

use Index_Functions, only: nTri_Elem1
use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: Shells
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: F(*)
integer(kind=iwp), intent(in) :: nBeta, ishll, lb, iAng, nveccb
integer(kind=iwp) :: ncb, nExpi
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)

!***********************************************************************

ncb = nTri_Elem1(lb)*nTri_Elem1(iang)
nExpi = Shells(iShll)%nExp
call mma_allocate(TMP1,nExpi*ncb*nVecCB*nBeta,Label='Tmp1')
call mma_allocate(TMP2,nExpi*ncb*nVecCB*nBeta,Label='Tmp2')

! ) kj,cbx -> cbx,kj

call DgeTMo(F,nBeta*nExpi,nBeta*nExpi,ncb*nVecCB,Tmp1,ncb*nVecCB)

! 2) bxkj,C = c,bxkj * c,C

call DGEMM_('T','N',nTri_Elem1(lb)*nVecCB*nExpi*nBeta,(2*iAng+1),nTri_Elem1(iAng),One,Tmp1,nTri_Elem1(iAng),RSph(ipSph(iAng)), &
            nTri_Elem1(iAng),Zero,Tmp2,nTri_Elem1(lb)*nVecCB*nExpi*nBeta)

! 3) bx,kjC -> kjC,bx

call DgeTMo(Tmp2,nTri_Elem1(lb)*nVecCB,nTri_Elem1(lb)*nVecCB,nExpi*nBeta*(2*iAng+1),F,nExpi*nBeta*(2*iAng+1))

call mma_deallocate(Tmp1)
call mma_deallocate(Tmp2)

return

end subroutine RToSph
