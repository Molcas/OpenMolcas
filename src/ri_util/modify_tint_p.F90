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

subroutine Modify_TInt_p(TInt,nTheta_All,List2,mData)

use Basis_Info, only: Shells
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nTheta_All, mData, List2(mData,nTheta_All)
real(kind=wp), intent(inout) :: TInt(nTheta_All,nTheta_All)
integer(kind=iwp) :: iShll, iTheta_All, jShll, jTheta_All, kShll, lShll, nConti, nContj, nContk, nContl, nPrimi, nPrimj, nPrimk, &
                     nPriml
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: iPrim, jPrim, kPrim, lPrim
#endif
real(kind=wp) :: Coeff_i, Coeff_j, Coeff_k, Coeff_l
real(kind=wp), external :: DDot_

#ifdef _DEBUGPRINT_
call RecPrt('Modify_TInt_p: TInt',' ',TInt,nTheta_All,nTheta_All)
#endif
do iTheta_All=1,nTheta_All
# ifdef _DEBUGPRINT_
  iPrim = List2(5,iTheta_All)
# endif
  iShll = List2(7,iTheta_All)
  nConti = Shells(iShll)%nBasis_C
  nPrimi = Shells(iShll)%nExp
  Coeff_i = DDot_(nConti,Shells(iShll)%Cff_c(1,1,1),nPrimi,Shells(iShll)%Cff_c(1,1,1),nPrimi)
  Coeff_i = sqrt(Coeff_i)

# ifdef _DEBUGPRINT_
  jPrim = List2(6,iTheta_All)
# endif
  jShll = List2(8,iTheta_All)
  nContj = Shells(jShll)%nBasis_C
  nPrimj = Shells(jShll)%nExp
  Coeff_j = DDot_(nContj,Shells(jShll)%Cff_c(1,1,1),nPrimj,Shells(jShll)%Cff_c(1,1,1),nPrimj)
  Coeff_j = sqrt(Coeff_j)

  do jTheta_All=1,nTheta_All
#   ifdef _DEBUGPRINT_
    kPrim = List2(5,jTheta_All)
#   endif
    kShll = List2(7,jTheta_All)
    nContk = Shells(kShll)%nBasis_C
    nPrimk = Shells(kShll)%nExp
    Coeff_k = DDot_(nContk,Shells(kShll)%Cff_c(1,1,1),nPrimk,Shells(kShll)%Cff_c(1,1,1),nPrimk)
    Coeff_k = sqrt(Coeff_k)

#   ifdef _DEBUGPRINT_
    lPrim = List2(6,jTheta_All)
#   endif
    lShll = List2(8,jTheta_All)
    nContl = Shells(lShll)%nBasis_C
    nPriml = Shells(lShll)%nExp
    Coeff_l = DDot_(nContl,Shells(lShll)%Cff_c(1,1,1),nPriml,Shells(lShll)%Cff_c(1,1,1),nPriml)
    Coeff_l = sqrt(Coeff_l)

    TInt(iTheta_All,jTheta_All) = TInt(iTheta_All,jTheta_All)*Coeff_i*Coeff_j*Coeff_k*Coeff_l
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) Coeff_i,Coeff_j,Coeff_k,Coeff_l
    write(u6,*) iPrim,jPrim,kPrim,lPrim
    write(u6,*) nPrimi,nPrimj,nPrimk,nPriml
#   endif

  end do

end do
#ifdef _DEBUGPRINT_
call RecPrt('Modify_TInt_p: TInt',' ',TInt,nTheta_All,nTheta_All)
#endif

return

end subroutine Modify_TInt_p
