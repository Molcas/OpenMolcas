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

subroutine GF(nX,nDoF,nInter,EVec,EVal,RedM,iNeg,dDipM,mTR,nAtom,DipM)

use Slapaf_Info, only: nDimBC, Smmtrc
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nX, nDoF, nInter, mTR, nAtom
real(kind=wp), intent(out) :: EVec(2,nDoF,nDoF), EVal(2,nDoF), RedM(nDoF)
integer(kind=iwp), intent(out) :: iNeg
real(kind=wp), intent(inout) :: dDipM(3,nInter+mTR)
real(kind=wp), intent(in) :: DipM(3)
integer(kind=iwp) :: i, iAtom, iNC, ix, ixyz, iy, iz
real(kind=wp), allocatable :: F(:), G(:), GInv(:), Tmp1(:), Tmp2(:)

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('GF: dDipM',' ',dDipM,3,nInter)
call RecPrt('GF: DipM',' ',DipM,3,1)
#endif
call mma_allocate(Tmp1,nX**2,Label='Tmp1')
call mma_allocate(Tmp2,nX**2,Label='Tmp2')
!                                                                      *
!***********************************************************************
!                                                                      *
! Note that all calculations will be done in the Cartesian basis!  *
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute harmonic frequencies
!                                                                      *
!***********************************************************************
!                                                                      *
! Generate the G matrix (mass tensor in cartesians)

call mma_allocate(G,nDimBC**2,Label='G')
call mma_allocate(GInv,nDimBC**2,Label='GInv')
call Mk_G(G,GInv,nDimBC)

! Get the force constant matrix in cartesians

call mma_allocate(F,nX**2,Label='F')
call Get_H(F,nX)
!                                                                      *
!***********************************************************************
!                                                                      *
! Form the GF-matrix (actually G^(1/2)FG^(1/2))

call GF_Mult(G,F,Tmp2,nDoF)  ! Result in Tmp2
call mma_deallocate(F)

! Compute the frequencies and harmonic eigenfunctions in Cartesians.

call GF_Harmonic_Frequencies(G,GInv,Tmp1,Tmp2,EVec,EVal,RedM,iNeg,nX,nDoF)

call mma_deallocate(G)
call mma_deallocate(GInv)
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Compute the dipole moment derivative in Cartesians.

call Get_dDipM(dDipM,DipM,nDoF,nInter)
!                                                                      *
!***********************************************************************
!                                                                      *
! Transform from cartesian to normal coordinates

do iNC=1,nDoF
  Tmp2(1:nDoF) = EVec(1,:,iNC)
  ix = (iNC-1)*3+1
  iy = (iNC-1)*3+2
  iz = (iNC-1)*3+3
  Tmp1(ix) = Zero
  Tmp1(iy) = Zero
  Tmp1(iz) = Zero
  i = 0
  do iAtom=1,nAtom
    do ixyz=1,3
      if (Smmtrc(ixyz,iAtom)) then
        i = i+1
        Tmp1(ix) = Tmp1(ix)+dDipM(1,i)*Tmp2(i)
        Tmp1(iy) = Tmp1(iy)+dDipM(2,i)*Tmp2(i)
        Tmp1(iz) = Tmp1(iz)+dDipM(3,i)*Tmp2(i)
      end if
    end do
  end do
end do
dDipM(:,:) = reshape(Tmp1(1:3*nDoF),[3,nDoF])
#ifdef _DEBUGPRINT_
call RecPrt('dDipM(normal coord.)',' ',dDipM,3,nDoF)
#endif
call mma_deallocate(Tmp2)
call mma_deallocate(Tmp1)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine GF
