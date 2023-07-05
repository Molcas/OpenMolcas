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

use Constants, only: Zero, One, autocm
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: iNeg, nX, nDoF
real(kind=wp) :: G(nX**2), GInv(nX**2), Tmp1(nDoF,nDoF), Tmp2(nX**2), EVec(2*nDoF,nDoF), EVal(2*nDoF), RedM(nDoF)
integer(kind=iwp) :: iHarm, iiT, iX, jHarm, jj, jX
real(kind=wp) :: r2, rlow, temp, Test_i, Test_j, tmp
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
! Compute eigenvectors and eigenfunctions

call FZero(Tmp1,nDoF**2)
call dcopy_(nX,[One],0,Tmp1,nDoF+1)
call NIDiag_new(Tmp2,Tmp1,nDoF,nDoF)

call FZero(EVal,2*nDoF)
call FZero(EVec,2*nDoF**2)

! Move over eigenvalue and eigenvectors, note that the eigenvectors
! are transformed back to Cartesian coordinates from mass-weighted
! Cartesians.

do iX=1,nDoF
  iiT = iX*(iX+1)/2
  EVal((iX-1)*2+1) = Tmp2(iiT)
  call dcopy_(nDoF,Tmp1(1,iX),1,EVec(1,iX),2)

  r2 = Zero
  do jX=1,nDoF
    jj = (jX-1)*nDoF+jX
    tmp = EVec((jX-1)*2+1,iX)
    tmp = tmp*sqrt(G(jj))
    EVec((jX-1)*2+1,iX) = tmp
    r2 = r2+tmp**2
  end do
  call DScal_(nDoF,One/sqrt(r2),EVec(1,iX),2)
end do
!define _DEBUGPRINT_
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
do iHarm=1,2*nDoF,2
  jHarm = (iHarm+1)/2
  temp = EVal(iHarm)

  ! Fix imaginary frequencies

  if (temp >= Zero) then
    EVal(jHarm) = sqrt(temp)*autocm
  else
    iNeg = iNeg+1
    EVal(jHarm) = -sqrt(abs(temp))*autocm
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('Converted EVal',' ',EVal,1,nDoF)
#endif

! Normalize over the metric of the masses

do iHarm=1,nDoF
  call dcopy_(nDoF,EVec(1,iHarm),2,Tmp1,1)
  call DGEMM_('N','N',nDoF,1,nDoF,One,GInv,nDoF,Tmp1,nDoF,Zero,Tmp2,nDoF)
  r2 = DDot_(nDoF,Tmp1,1,Tmp2,1)
  RedM(iHarm) = r2
  r2 = One/sqrt(r2)
  call DScal_(nDoF,r2,EVec(1,iHarm),2)
end do
#ifdef _DEBUGPRINT_
call RecPrt('Normal coordinates (Q)',' ',EVec,nDoF*2,nDoF)
#endif

! Order, from low to high. Put translations and rotations last.

do iHarm=1,nDoF-1
  Test_i = EVal(iHarm)
  if (abs(Test_i) < 1.0e-3_wp) Test_i = 1.0e5_wp
  do jHarm=iHarm+1,nDoF
    Test_j = EVal(jHarm)
    if (abs(Test_j) < 1.0e-3_wp) Test_j = 1.0e5_wp
    if (Test_j < Test_i) then
      rlow = Test_i
      Test_i = Test_j
      Test_j = rlow
      rlow = EVal(iHarm)
      EVal(iHarm) = EVal(jHarm)
      EVal(jHarm) = rLow
      rlow = RedM(iHarm)
      RedM(iHarm) = RedM(jHarm)
      RedM(jHarm) = rLow
      call DSwap_(nDoF,EVec(1,iHarm),2,EVec(1,jHarm),2)
    end if
  end do
end do
#ifdef _DEBUGPRINT_
call RecPrt('Frequencies (cm-1)',' ',EVal,1,nDoF)
call RecPrt('Reduced masses (u)',' ',RedM,1,nDoF)
call RecPrt('Normal Coordinates',' ',EVec,nDoF*2,nDoF)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine GF_Harmonic_Frequencies
