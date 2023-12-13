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

subroutine RToCore(F,nBeta,ishll,lb,iAng,nveccb)

!***********************************************************************
!
! Transformation kernel to atomic orbitals in normalized spherical harmonics
!
!***********************************************************************
! @parameter F The cartesian components of <A|core>
! @parameter nBeta Number of exponents
! @parameter ishll Shell number for ECP
! @parameter lb angular momenta Ket
! @parameter iAng angular momenta core
! @parameter nveccb Number of derivatives
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
integer(kind=iwp) :: nBasisi, ncb, nExpi
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)

!***********************************************************************
ncb = nTri_Elem1(lb)*nTri_Elem1(iang)
nExpi = Shells(iShll)%nExp
nBasisi = Shells(iShll)%nBasis
call mma_allocate(TMP1,nExpi*ncb*nVecCB*nBeta,Label='Tmp1')
call mma_allocate(TMP2,nExpi*ncb*nVecCB*nBeta,Label='Tmp2')

! And (almost) the same thing for the righthand side, form
! KjCb from kjcb
! 1) jcb,K = k,jcb * k,K

call DGEMM_('T','N',nBeta*ncb*nVecCB,nBasisi,nExpi,One,F,nExpi,Shells(iShll)%pCff,nExpi,Zero,Tmp1,nBeta*ncb*nVecCB)

! 2) j,cbK -> cbK,j

call DgeTMo(Tmp1,nBeta,nBeta,ncb*nVecCB*nBasisi,Tmp2,ncb*nVecCB*nBasisi)

! 3) bKj,C = c,bKj * c,C

call DGEMM_('T','N',nTri_Elem1(lb)*nVecCB*nBasisi*nBeta,(2*iAng+1),nTri_Elem1(iAng),One,Tmp2,nTri_Elem1(iAng),RSph(ipSph(iAng)), &
            nTri_Elem1(iAng),Zero,Tmp1,nTri_Elem1(lb)*nVecCB*nBasisi*nBeta)

! 4) b,KjC -> KjC,b

call DgeTMo(Tmp1,nTri_Elem1(lb)*nVecCB,nTri_Elem1(lb)*nVecCB,nBasisi*nBeta*(2*iAng+1),F,nBasisi*nBeta*(2*iAng+1))

call mma_deallocate(Tmp2)
call mma_deallocate(Tmp1)

return

end subroutine RToCore
