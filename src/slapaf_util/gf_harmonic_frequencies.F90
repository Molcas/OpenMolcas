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

subroutine GF_Harmonic_Frequencies(G,GInv,Tmp1,Tmp2,EVec,EVal,RedM,iNeg,nX,nDoF)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, One, autocm
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nX, nDoF
real(kind=wp), intent(in) :: G(nX**2), GInv(nX**2)
real(kind=wp), intent(out) :: Tmp1(nDoF,nDoF), EVec(2,nDoF,nDoF), EVal(2,nDoF), RedM(nDoF)
real(kind=wp), intent(inout) :: Tmp2(nX**2)
integer(kind=iwp), intent(out) :: iNeg
integer(kind=iwp) :: iHarm, iiT, iX, jHarm, jj, jX
real(kind=wp) :: r2, rlow, temp, Test_i, Test_j, tmp
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute eigenvectors and eigenfunctions

call unitmat(Tmp1,nDoF)
call NIDiag_new(Tmp2,Tmp1,nDoF,nDoF)

EVal(:,:) = Zero
EVec(:,:,:) = Zero

! Move over eigenvalue and eigenvectors, note that the eigenvectors
! are transformed back to Cartesian coordinates from mass-weighted
! Cartesians.

do iX=1,nDoF
  iiT = nTri_Elem(iX)
  EVal(1,iX) = Tmp2(iiT)
  EVec(1,:,iX) = Tmp1(:,iX)

  r2 = Zero
  do jX=1,nDoF
    jj = (jX-1)*nDoF+jX
    tmp = EVec(1,jX,iX)
    tmp = tmp*sqrt(G(jj))
    EVec(1,jX,iX) = tmp
    r2 = r2+tmp**2
  end do
  EVec(1,:,iX) = EVec(1,:,iX)/sqrt(r2)
end do
!#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
call RecPrt('EVal',' ',EVal,2,nDoF)
call RecPrt('EVec',' ',EVec,nDoF*2,nDoF)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the harmonic frequencies
!                                                                      *
!***********************************************************************
!                                                                      *
iNeg = 0
do iHarm=1,nDoF
  temp = EVal(1,iHarm)

  ! Fix imaginary frequencies

  if (temp >= Zero) then
    EVal(1,iHarm) = sqrt(temp)*autocm
  else
    iNeg = iNeg+1
    EVal(1,iHarm) = -sqrt(abs(temp))*autocm
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('Converted EVal',' ',EVal(1,:),1,nDoF)
#endif

! Normalize over the metric of the masses

do iHarm=1,nDoF
  Tmp1(:,1) = EVec(1,:,iHarm)
  call DGEMM_('N','N',nDoF,1,nDoF,One,GInv,nDoF,Tmp1,nDoF,Zero,Tmp2,nDoF)
  r2 = DDot_(nDoF,Tmp1,1,Tmp2,1)
  RedM(iHarm) = r2
  EVec(1,:,iHarm) = EVec(1,:,iHarm)/sqrt(r2)
end do
#ifdef _DEBUGPRINT_
call RecPrt('Normal coordinates (Q)',' ',EVec,nDoF*2,nDoF)
#endif

! Order, from low to high. Put translations and rotations last.

do iHarm=1,nDoF-1
  Test_i = EVal(1,iHarm)
  if (abs(Test_i) < 1.0e-3_wp) Test_i = 1.0e5_wp
  do jHarm=iHarm+1,nDoF
    Test_j = EVal(1,jHarm)
    if (abs(Test_j) < 1.0e-3_wp) Test_j = 1.0e5_wp
    if (Test_j < Test_i) then
      rlow = Test_i
      Test_i = Test_j
      Test_j = rlow
      rlow = EVal(1,iHarm)
      EVal(1,iHarm) = EVal(1,jHarm)
      EVal(1,jHarm) = rLow
      rlow = RedM(iHarm)
      RedM(iHarm) = RedM(jHarm)
      RedM(jHarm) = rLow
      do iX=1,nDoF
        tmp = EVec(1,iX,iHarm)
        EVec(1,iX,iHarm) = EVec(1,iX,jHarm)
        EVec(1,iX,jHarm) = tmp
      end do
    end if
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Frequencies (cm-1)',' ',EVal(1,:),1,nDoF)
call RecPrt('Reduced masses (u)',' ',RedM,1,nDoF)
call RecPrt('Normal Coordinates',' ',EVec,nDoF*2,nDoF)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine GF_Harmonic_Frequencies
