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

subroutine LToCore(F,nAlpha,iShll,la,iAng,nvecac)
!******************************************************************************
!
! Transformation kernel to atomic orbials in normalized spherical harmonics
!
!******************************************************************************
! @parameter F  The cartesian components of <A|core>
! @parameter nAlpha Number of exponents
! @parameter ishll Shell number for ECP
! @parameter la angular momenta LS
! @parameter iAng angular momenta core
! @parameter nVecAC Number of derivatives

use Index_Functions, only: nTri_Elem1
use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: Shells
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: F(*)
integer(kind=iwp), intent(in) :: nAlpha, iShll, la, iAng, nVecAC
integer(kind=iwp) :: iBk, n, nac, nBasisi, nExpi
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)

nac = nTri_Elem1(la)*nTri_Elem1(iang)
nExpi = Shells(iShll)%nExp
nBasisi = Shells(iShll)%nBasis
call mma_allocate(Tmp1,nExpi*nac*nVecAC*nalpha,Label='Tmp1')
call mma_allocate(Tmp2,nExpi*nac*nVecAC*nalpha,Label='Tmp2')
! From the lefthandside overlap, form iKaC from ikac by
! 1) i,kac -> k,aci

n = nExpi*nac*nVecAC
call DgeTMo(F,nAlpha,nAlpha,n,Tmp1,n)

! 2) aciK =  k,aci * k,K (Contract over core orbital)

n = nac*nVecAC*nAlpha
call DGEMM_('T','N',n,nBasisi,nExpi,One,Tmp1,nExpi,Shells(iShll)%pCff,nExpi,Zero,Tmp2,n)

! 3) Mult by shiftoperators aci,K -> Bk(K) * aci,K

do iBk=1,nBasisi
  Tmp1((iBk-1)*n+1:iBk*n) = Shells(iShll)%Bk(iBk)*Tmp2((iBk-1)*n+1:iBk*n)
end do

! 4) a,ciK -> ciKa

call DgeTMo(Tmp1,nTri_Elem1(la),nTri_Elem1(la),nTri_Elem1(iAng)*nVecAC*nAlpha*nBasisi,Tmp2,nTri_Elem1(iAng)*nVecAC*nAlpha*nBasisi)

! 5) iKa,C = c,iKa * c,C

call DGEMM_('T','N',nVecAC*nAlpha*nBasisi*nTri_Elem1(la),(2*iAng+1),nTri_Elem1(iAng),One,Tmp2,nTri_Elem1(iAng),RSph(ipSph(iAng)), &
            nTri_Elem1(iAng),Zero,Tmp1,nVecAC*nAlpha*nBasisi*nTri_Elem1(la))

call DgeTMo(Tmp1,nVecAC,nVecAC,nAlpha*nBasisi*nTri_Elem1(la)*(2*iAng+1),F,nAlpha*nBasisi*nTri_Elem1(la)*(2*iAng+1))

call mma_deallocate(Tmp2)
call mma_deallocate(Tmp1)

return

end subroutine LToCore
