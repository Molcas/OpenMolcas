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

subroutine LToSph(F,nalpha,ishll,la,iAng,nvecac)
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
! @parameter ishll Shell number for ECP
! @parameter la angular momenta LS
! @parameter iAng angular momenta core
! @parameter Number of derivatives
!***********************************************************************

use Basis_Info
use Real_Spherical

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
real*8 F(*)
real*8, allocatable :: Tmp1(:), Tmp2(:)
!Statement function
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

!***********************************************************************

nExpi = Shells(iShll)%nExp
nac = nelem(la)*nelem(iang)
call mma_allocate(Tmp1,nExpi*nac*nVecAC*nalpha,Label='Tmp1')
call mma_allocate(Tmp2,nExpi*nac*nVecAC*nalpha,Label='Tmp2')

call DgeTMo(F,nAlpha*nExpi*nElem(la),nAlpha*nExpi*nElem(la),nElem(iAng)*nVecAC,Tmp1,nElem(iAng)*nVecAC)

! 2) xika,C = c,xika * c,C

call DGEMM_('T','N',nVecAC*nAlpha*nExpi*nElem(la),(2*iAng+1),nElem(iAng),One,Tmp1,nElem(iAng),RSph(ipSph(iAng)),nElem(iAng),Zero, &
            Tmp2,nVecAC*nAlpha*nExpi*nElem(la))

! 3) x,ikaC -> ikaC,x

call DGetMo(Tmp2,nVecAC,nVecAC,nAlpha*nExpi*nElem(la)*(2*iAng+1),Tmp1,nAlpha*nExpi*nElem(la)*(2*iAng+1))
call dcopy_(nVecAC*nAlpha*nExpi*nElem(la)*(2*iAng+1),Tmp1,1,F,1)

call mma_deallocate(Tmp2)
call mma_deallocate(Tmp1)

return

end subroutine LToSph
