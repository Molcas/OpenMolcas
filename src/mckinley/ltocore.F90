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

use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: Shells
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: F(*)
integer(kind=iwp) :: nAlpha, iShll, la, iAng, nVecAC
integer(kind=iwp) :: iBk, n, nac, nBasisi, nExpi
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)
! Statement function
integer(kind=iwp) :: nElem, ixyz
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

nac = nelem(la)*nelem(iang)
nExpi = Shells(iShll)%nExp
nBasisi = Shells(iShll)%nBasis
call mma_allocate(Tmp1,nExpi*nac*nVecAC*nalpha,Label='Tmp1')
call mma_allocate(Tmp2,nExpi*nac*nVecAC*nalpha,Label='Tmp2')
! From the lefthandside overlap, form iKaC from ikac by
! 1) i,kac -> k,aci

n = nExpi*nac*nVecAC
call DgeTMo(F,nAlpha,nAlpha,n,Tmp1,n)

! 2) aciK =  k,aci * k,K (Contract over core orbital)

call DGEMM_('T','N',nac*nVecAC*nAlpha,nBasisi,nExpi,One,Tmp1,nExpi,Shells(iShll)%pCff,nExpi,Zero,Tmp2,nac*nVecAC*nAlpha)

! 3) Mult by shiftoperators aci,K -> Bk(K) * aci,K

do iBk=1,nBasisi
  call DYaX(nac*nVecAC*nAlpha,Shells(iShll)%Bk(iBk),Tmp2((iBk-1)*nac*nVecAC*nAlpha+1),1,Tmp1((iBk-1)*nac*nVecAC*nAlpha+1),1)
end do

! 4) a,ciK -> ciKa

call DgeTMo(Tmp1,nElem(la),nElem(la),nElem(iAng)*nVecAC*nAlpha*nBasisi,Tmp2,nElem(iAng)*nVecAC*nAlpha*nBasisi)

! 5) iKa,C = c,iKa * c,C

call DGEMM_('T','N',nVecAC*nAlpha*nBasisi*nElem(la),(2*iAng+1),nElem(iAng),One,Tmp2,nElem(iAng),RSph(ipSph(iAng)),nElem(iAng), &
            Zero,Tmp1,nVecAC*nAlpha*nBasisi*nElem(la))

call DgeTMo(Tmp1,nVecAC,nVecAC,nAlpha*nBasisi*nElem(la)*(2*iAng+1),F,nAlpha*nBasisi*nElem(la)*(2*iAng+1))

call mma_deallocate(Tmp2)
call mma_deallocate(Tmp1)

return

end subroutine LToCore
