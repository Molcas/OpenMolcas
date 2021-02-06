!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2004, Per-Olof Widmark                                 *
!***********************************************************************
!***********************************************************************
!                                                                      *
! This routine will populate according to the aufbau principle.        *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: Oct 2004                                                    *
!                                                                      *
!***********************************************************************

subroutine GoPop(Eps,Occ,Scr,n,PrtEor,PrThr,GapThr)

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: Eps(n), PrThr, GapThr
real(kind=wp), intent(inout) :: Occ(n)
real(kind=wp), intent(out) :: Scr(n)
logical(kind=iwp), intent(in) :: PrtEor
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
!real(kind=wp) ::  eFermi
real(kind=wp) :: TotNucChg, eGap, OccNo, eLo, eHi, tmp
integer(kind=iwp) :: nElec, nAlpha, nBeta, nOcc, nAct, kLo, kHi, i, j, k, m
!----------------------------------------------------------------------*
! Some setup                                                           *
!----------------------------------------------------------------------*
!GapThr = 1.0e-2_wp
!----------------------------------------------------------------------*
! Sort orbital energies                                                *
!----------------------------------------------------------------------*
do i=1,n
  Scr(i) = Eps(i)
end do
do i=1,n
  j = i
  do k=i,n
    if (Scr(k) < Scr(j)) j = k
  end do
  tmp = Scr(i)
  Scr(i) = Scr(j)
  Scr(j) = tmp
end do
!----------------------------------------------------------------------*
! How many alpha/beta                                                  *
!----------------------------------------------------------------------*
call Get_dScalar('Total nuclear Charge',TotNucChg)
nElec = int(TotNucChg+Half)
nBeta = int(Half*(TotNucChg+Half))
nAlpha = nElec-nBeta
!write(u6,'(a,i5)') 'nElec . . . . . . . . . .',nElec
!write(u6,'(a,i5)') 'nAlpha  . . . . . . . . .',nAlpha
!write(u6,'(a,i5)') 'nBeta . . . . . . . . . .',nBeta
!----------------------------------------------------------------------*
! Optonally print sorted orbitals energies                             *
!----------------------------------------------------------------------*
if (PrtEor) then
  m = 0
  do i=1,n
    if (Scr(i) <= PrThr) m = i
  end do
  write(u6,*)
  write(u6,'(a)') 'Sorted orbital energies'
  write(u6,'(a)') '-----------------------'
  write(u6,*)
  write(u6,'(a,i5,a,i5)') 'Printing',m,' out of',n
  write(u6,'(a,f6.1)') 'Filled orbitals:',Half*TotNucChg
  write(u6,*)
  write(u6,'(i5,1h-,i5,2x,10f12.4)') (i,min(i+9,m),(Scr(j),j=i,min(i+9,m)),i=1,m,10)
  write(u6,*)
end if
!----------------------------------------------------------------------*
! Populate alpha                                                       *
!----------------------------------------------------------------------*
if (nAlpha >= n) then
  !write(u6,'(a)') 'Alpha is MB'
  eLo = min(Eps(n)+1.0e-6_wp,Zero)
  eHi = eLo
  OccNo = Zero
else if (nAlpha <= 0) then
  !write(u6,'(a)') 'Alpha is empty'
  eLo = Eps(1)-One
  eHi = eLo
  OccNo = Zero
else
  !write(u6,'(a)') 'Alpha is not MB'
  if (nAlpha > 0) then
    !eFermi = (Scr(nAlpha+1)+Scr(nAlpha))/Two
    eGap = (Scr(nAlpha+1)-Scr(nAlpha))
  else
    !eFermi = Zero
    eGap = Zero
  end if
  !write(u6,'(a,f12.6)') 'eFermi (alpha)  . . . . .',eFermi
  !write(u6,'(a,f12.6)') 'eGap (alpha)  . . . . . .',eGap
  if (eGap > GapThr) then
    !write(u6,'(a)') 'Alpha have large gap'
    OccNo = Zero
    eLo = 0.25_wp*Scr(nAlpha)+0.75_wp*Scr(nAlpha+1)
    eHi = 0.75_wp*Scr(nAlpha)+0.25_wp*Scr(nAlpha+1)
  else
    !write(u6,'(a)') 'Alpha have small gap'
    kLo = 1
    do i=2,nAlpha
      if (Scr(i)-Scr(i-1) > GapThr) kLo = i
    end do
    kHi = n
    do i=n-1,max(1,nAlpha),-1
      if (Scr(i+1)-Scr(i) > GapThr) kHi = i
    end do
    if (kLo > 1) then
      eLo = (Scr(kLo)+Scr(kLo-1))/Two
    else
      eLo = Scr(1)-One
    end if
    if (kHi < n) then
      eHi = (Scr(kHi)+Scr(kHi+1))/Two
    else
      eHi = Scr(n)+One
    end if
    !write(u6,'(a,i5)') 'kLo (alpha) . . . . . . .',kLo
    !write(u6,'(a,i5)') 'kHi (alpha) . . . . . . .',kHi
    !write(u6,'(a,f12.6)') 'eLo (alpha) . . . . . . .',eLo
    !write(u6,'(a,f12.6)') 'eHi (alpha) . . . . . . .',eHi
    nOcc = kLo-1
    nAct = kHi-kLo+1
    OccNo = One*(nAlpha-nOcc)/nAct
  end if
  !write(u6,'(a,f12.6)') 'OccNo (alpha) . . . . . .',OccNo
end if
do i=1,n
  if (Eps(i) < eLo) then
    Occ(i) = Occ(i)+One
  else if (Eps(i) < eHi) then
    Occ(i) = Occ(i)+OccNo
  end if
end do
!----------------------------------------------------------------------*
! Populate beta                                                        *
!----------------------------------------------------------------------*
if (nBeta >= n) then
  !write(u6,'(a)') 'Beta is MB'
  eLo = min(Eps(n)+1.0e-6_wp,Zero)
  eHi = eLo
  OccNo = Zero
else if (nBeta <= 0) then
  !write(u6,'(a)') 'Beta is empty'
  eLo = Eps(1)-One
  eHi = eLo
  OccNo = Zero
else
  !write(u6,'(a)') 'Beta is not MB'
  if (nBeta > 0) then
    !eFermi = (Scr(nBeta+1)+Scr(nBeta))/Two
    eGap = (Scr(nBeta+1)-Scr(nBeta))
  else
    !eFermi = Zero
    eGap = Zero
  end if
  !write(u6,'(a,f12.6)') 'eFermi (beta) . . . . . .',eFermi
  !write(u6,'(a,f12.6)') 'eGap (beta) . . . . . . .',eGap
  if (eGap > GapThr) then
    !write(u6,'(a)') 'Beta have large gap'
    OccNo = Zero
    eLo = 0.25_wp*Scr(nBeta)+0.75_wp*Scr(nBeta+1)
    eHi = 0.75_wp*Scr(nBeta)+0.25_wp*Scr(nBeta+1)
  else
    !write(u6,'(a)') 'Beta have small gap'
    kLo = 1
    do i=2,nBeta
      if (Scr(i)-Scr(i-1) > GapThr) kLo = i
    end do
    kHi = n
    do i=n-1,max(1,nBeta),-1
      if (Scr(i+1)-Scr(i) > GapThr) kHi = i
    end do
    if (kLo > 1) then
      eLo = (Scr(kLo)+Scr(kLo-1))/Two
    else
      eLo = Scr(1)-One
    end if
    if (kHi < n) then
      eHi = (Scr(kHi)+Scr(kHi+1))/Two
    else
      eHi = Scr(n)+One
    end if
    !write(u6,'(a,i5)') 'kLo (beta)  . . . . . . .',kLo
    !write(u6,'(a,i5)') 'kHi (beta)  . . . . . . .',kHi
    !write(u6,'(a,f12.6)') 'eLo (beta)  . . . . . . .',eLo
    !write(u6,'(a,f12.6)') 'eHi (beta)  . . . . . . .',eHi
    nOcc = kLo-1
    nAct = kHi-kLo+1
    OccNo = One*(nBeta-nOcc)/nAct
  end if
  !write(u6,'(a,f12.6)') 'OccNo (beta)  . . . . . .',OccNo
end if
do i=1,n
  if (Eps(i) < eLo) then
    Occ(i) = Occ(i)+One
  else if (Eps(i) < eHi) then
    Occ(i) = Occ(i)+OccNo
  end if
end do
!----------------------------------------------------------------------*
! Print population (debug)                                             *
!----------------------------------------------------------------------*
!write(u6,*)
!write(u6,'(a)') 'Occupation of orbitals'
!write(u6,'(a)') '----------------------'
!write(u6,*)
!write(u6,'(a,i5,a,i5)') 'Printing',m,' out of',n
!write(u6,*)
!write(u6,'(i5,1h-,i5,2x,10f12.4)') (i,Min(i+9,m),(Occ(j),j=i,Min(i+9,m)),i=1,m,10)
!----------------------------------------------------------------------*
! Done!                                                                *
!----------------------------------------------------------------------*
return

end subroutine GoPop
