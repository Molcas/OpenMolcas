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

subroutine One_Int(Kernel,Array,nArray,A,iAng,iComp,nOrdOp,Scr1,nScr1,Scr2,nScr2,naa,SAR,nSAR,iShll_a,nPrim_a,Exp_a,nCntrc_a, &
                   Cff_a,iCmp_a,iShll_r,nPrim_r,Exp_r,nCntrc_r,Cff_r,iCmp_r)

use Basis_Info, only: Shells
use Real_Spherical, only: ipSph, RSph
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
external :: Kernel
integer(kind=iwp), intent(in) :: nArray, iAng, iComp, nOrdOp, nScr1, nScr2, naa, nSAR, iShll_a, nPrim_a, nCntrc_a, iCmp_a, &
                                 iShll_r, nPrim_r, nCntrc_r, iCmp_r
real(kind=wp), intent(out) :: Array(nArray), Scr1(nScr1), Scr2(nScr2), SAR(nSAR)
real(kind=wp), intent(in) :: A(3), Exp_a(nPrim_a), Cff_a(nPrim_a,nCntrc_a), Exp_r(nPrim_r), Cff_r(nPrim_r,nCntrc_r)
integer(kind=iwp) :: mArr, mSar, nHer
real(kind=wp), allocatable :: KAR(:), PAR(:,:), pSAR(:), ZAR(:), ZIAR(:)

mArr = nArray/(nPrim_a*nPrim_r)
call mma_allocate(ZAR,nPrim_a*nPrim_r,Label='ZAR')
call mma_allocate(ZIAR,nPrim_a*nPrim_r,Label='ZIAR')
call mma_allocate(KAR,nPrim_a*nPrim_r,Label='KAR')
call mma_allocate(PAR,nPrim_a*nPrim_r,3,Label='PAR')

mSAR = nPrim_a*nPrim_r*naa
call mma_allocate(pSAR,mSAR,Label='pSAR')

call ZXia(ZAR,ZIAR,nPrim_a,nPrim_r,Exp_a,Exp_r)
call SetUp1(Exp_a,nPrim_a,Exp_r,nPrim_r,A,A,KAR,PAR,ZIAR)

nHer = (2*iAng+2+nOrdOp)/2
call Kernel(Exp_a,nPrim_a,Exp_r,nPrim_r,ZAR,ZIAR,KAR,PAR,pSAR,nPrim_a*nPrim_r,iComp,iAng,iAng,A,A,nHer,Array,mArr,A,nOrdOp)

call mma_deallocate(ZAR)
call mma_deallocate(ZIAR)
call mma_deallocate(KAR)
call mma_deallocate(PAR)

call DGEMM_('T','N',nPrim_r*naa,nCntrc_a,nPrim_a,One,pSAR,nPrim_a,Cff_a,nPrim_a,Zero,Scr1,nPrim_r*naa)
call DGEMM_('T','N',naa*nCntrc_a,nCntrc_r,nPrim_r,One,Scr1,nPrim_r,Cff_r,nPrim_r,Zero,Scr2,naa*nCntrc_a)
#ifdef _DEBUGPRINT_
call RecPrt('S_AR in Cartesian',' ',Scr2,naa,nCntrc_a*nCntrc_r)
#endif

! Transform to spherical Gaussian if needed!

if (Shells(iShll_a)%Transf .or. Shells(iShll_r)%Transf) then

  call CarSph(Scr2,naa,nCntrc_a*nCntrc_r,pSAR,nScr2,RSph(ipSph(iAng)),iAng,Shells(iShll_a)%Transf,Shells(iShll_a)%Prjct, &
              RSph(ipSph(iAng)),iAng,Shells(iShll_r)%Transf,Shells(iShll_r)%Prjct,SAR,iCmp_a*iCmp_r)
else
  call DGeTmO(Scr2,naa,naa,nCntrc_a*nCntrc_r,SAR,nCntrc_a*nCntrc_r)
end if
call mma_deallocate(pSAR)
#ifdef _DEBUGPRINT_
call RecPrt('S_AR in Sphericals',' ',SAR,iCmp_a*iCmp_r,nCntrc_a*nCntrc_r)
#endif

return

end subroutine One_Int
