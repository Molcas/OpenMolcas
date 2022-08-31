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

subroutine LToSph(F,nAlpha,iShll,la,iAng,nVecAC)
!***********************************************************************
!
!   Transform <A|core> from Cartesian components to Sperical harmonics
!
!           Observe that as opposed to the projection operator that this
!           contraction is done in the primitive basis.
!
!***********************************************************************
! @parameter F  The cartesian components of <A|core>(in)
!               The spherical components of <A|core>(out)
! @parameter nAlpha Number of exponents
! @parameter iShll Shell number for ECP
! @parameter la angular momenta LS
! @parameter iAng angular momenta core
! @parameter nVecAC Number of derivatives
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Real_Spherical, only: ipSph, RSph
use Basis_Info, only: Shells
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: F(*)
integer(kind=iwp), intent(in) :: nAlpha, iShll, la, iAng, nVecAC
integer(kind=iwp) :: nac, nExpi
real(kind=wp), allocatable :: Tmp1(:), Tmp2(:)

!***********************************************************************

nExpi = Shells(iShll)%nExp
nac = nTri_Elem1(la)*nTri_Elem1(iang)
call mma_allocate(Tmp1,nExpi*nac*nVecAC*nalpha,Label='Tmp1')
call mma_allocate(Tmp2,nExpi*nac*nVecAC*nalpha,Label='Tmp2')

call DgeTMo(F,nAlpha*nExpi*nTri_Elem1(la),nAlpha*nExpi*nTri_Elem1(la),nTri_Elem1(iAng)*nVecAC,Tmp1,nTri_Elem1(iAng)*nVecAC)

! 2) xika,C = c,xika * c,C

call DGEMM_('T','N',nVecAC*nAlpha*nExpi*nTri_Elem1(la),(2*iAng+1),nTri_Elem1(iAng),One,Tmp1,nTri_Elem1(iAng),RSph(ipSph(iAng)), &
            nTri_Elem1(iAng),Zero,Tmp2,nVecAC*nAlpha*nExpi*nTri_Elem1(la))

! 3) x,ikaC -> ikaC,x

call DGetMo(Tmp2,nVecAC,nVecAC,nAlpha*nExpi*nTri_Elem1(la)*(2*iAng+1),F,nAlpha*nExpi*nTri_Elem1(la)*(2*iAng+1))

call mma_deallocate(Tmp2)
call mma_deallocate(Tmp1)

return

end subroutine LToSph
