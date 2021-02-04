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

implicit none
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
real*8 Eps(*)
real*8 Occ(*)
real*8 Scr(*)
integer n
logical PrtEor
real*8 PrThr
real*8 GapThr
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
real*8 TotNucChg
!real*8 eFermi
real*8 eGap
real*8 OccNo
real*8 eLo, eHi
real*8 tmp
integer nElec
integer nAlpha
integer nBeta
integer nOcc
integer nAct
integer kLo, kHi
integer i, j, k, m
!----------------------------------------------------------------------*
! Some setup                                                           *
!----------------------------------------------------------------------*
!GapThr = 1.0d-2
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
nElec = int(TotNucChg+0.5d0)
nBeta = int(0.5d0*(TotNucChg+0.5d0))
nAlpha = nElec-nBeta
!write(6,'(a,i5)') 'nElec . . . . . . . . . .',nElec
!write(6,'(a,i5)') 'nAlpha  . . . . . . . . .',nAlpha
!write(6,'(a,i5)') 'nBeta . . . . . . . . . .',nBeta
!----------------------------------------------------------------------*
! Optonally print sorted orbitals energies                             *
!----------------------------------------------------------------------*
if (PrtEor) then
  m = 0
  do i=1,n
    if (Scr(i) <= PrThr) m = i
  end do
  write(6,*)
  write(6,'(a)') 'Sorted orbital energies'
  write(6,'(a)') '-----------------------'
  write(6,*)
  write(6,'(a,i5,a,i5)') 'Printing',m,' out of',n
  write(6,'(a,f6.1)') 'Filled orbitals:',0.5d0*TotNucChg
  write(6,*)
  write(6,'(i5,1h-,i5,2x,10f12.4)') (i,min(i+9,m),(Scr(j),j=i,min(i+9,m)),i=1,m,10)
  write(6,*)
end if
!----------------------------------------------------------------------*
! Populate alpha                                                       *
!----------------------------------------------------------------------*
if (nAlpha >= n) then
  !write(6,'(a)') 'Alpha is MB'
  eLo = min(Eps(n)+1.0d-6,0.0d0)
  eHi = eLo
  OccNo = 0.0d0
else if (nAlpha <= 0) then
  !write(6,'(a)') 'Alpha is empty'
  eLo = Eps(1)-1.0d0
  eHi = eLo
  OccNo = 0.0d0
else
  !write(6,'(a)') 'Alpha is not MB'
  if (nAlpha > 0) then
    !eFermi = (Scr(nAlpha+1)+Scr(nAlpha))/2.0d0
    eGap = (Scr(nAlpha+1)-Scr(nAlpha))
  else
    !eFermi = 0.0D0
    eGap = 0.0d0
  end if
  !write(6,'(a,f12.6)') 'eFermi (alpha)  . . . . .',eFermi
  !write(6,'(a,f12.6)') 'eGap (alpha)  . . . . . .',eGap
  if (eGap > GapThr) then
    !write(6,'(a)') 'Alpha have large gap'
    OccNo = 0.0d0
    eLo = 0.25d0*Scr(nAlpha)+0.75d0*Scr(nAlpha+1)
    eHi = 0.75d0*Scr(nAlpha)+0.25d0*Scr(nAlpha+1)
  else
    !write(6,'(a)') 'Alpha have small gap'
    kLo = 1
    do i=2,nAlpha
      if (Scr(i)-Scr(i-1) > GapThr) kLo = i
    end do
    kHi = n
    do i=n-1,max(1,nAlpha),-1
      if (Scr(i+1)-Scr(i) > GapThr) kHi = i
    end do
    if (kLo > 1) then
      eLo = (Scr(kLo)+Scr(kLo-1))/2.0d0
    else
      eLo = Scr(1)-1.0d0
    end if
    if (kHi < n) then
      eHi = (Scr(kHi)+Scr(kHi+1))/2.0d0
    else
      eHi = Scr(n)+1.0d0
    end if
    !write(6,'(a,i5)') 'kLo (alpha) . . . . . . .',kLo
    !write(6,'(a,i5)') 'kHi (alpha) . . . . . . .',kHi
    !write(6,'(a,f12.6)') 'eLo (alpha) . . . . . . .',eLo
    !write(6,'(a,f12.6)') 'eHi (alpha) . . . . . . .',eHi
    nOcc = kLo-1
    nAct = kHi-kLo+1
    OccNo = 1.0d0*(nAlpha-nOcc)/nAct
  end if
  !write(6,'(a,f12.6)') 'OccNo (alpha) . . . . . .',OccNo
end if
do i=1,n
  if (Eps(i) < eLo) then
    Occ(i) = Occ(i)+1.0d0
  else if (Eps(i) < eHi) then
    Occ(i) = Occ(i)+OccNo
  end if
end do
!----------------------------------------------------------------------*
! Populate beta                                                        *
!----------------------------------------------------------------------*
if (nBeta >= n) then
  !write(6,'(a)') 'Beta is MB'
  eLo = min(Eps(n)+1.0d-6,0.0d0)
  eHi = eLo
  OccNo = 0.0d0
else if (nBeta <= 0) then
  !write(6,'(a)') 'Beta is empty'
  eLo = Eps(1)-1.0d0
  eHi = eLo
  OccNo = 0.0d0
else
  !write(6,'(a)') 'Beta is not MB'
  if (nBeta > 0) then
    !eFermi = (Scr(nBeta+1)+Scr(nBeta))/2.0d0
    eGap = (Scr(nBeta+1)-Scr(nBeta))
  else
    !eFermi = 0.0D0
    eGap = 0.0d0
  end if
  !write(6,'(a,f12.6)') 'eFermi (beta) . . . . . .',eFermi
  !write(6,'(a,f12.6)') 'eGap (beta) . . . . . . .',eGap
  if (eGap > GapThr) then
    !write(6,'(a)') 'Beta have large gap'
    OccNo = 0.0d0
    eLo = 0.25d0*Scr(nBeta)+0.75d0*Scr(nBeta+1)
    eHi = 0.75d0*Scr(nBeta)+0.25d0*Scr(nBeta+1)
  else
    !write(6,'(a)') 'Beta have small gap'
    kLo = 1
    do i=2,nBeta
      if (Scr(i)-Scr(i-1) > GapThr) kLo = i
    end do
    kHi = n
    do i=n-1,max(1,nBeta),-1
      if (Scr(i+1)-Scr(i) > GapThr) kHi = i
    end do
    if (kLo > 1) then
      eLo = (Scr(kLo)+Scr(kLo-1))/2.0d0
    else
      eLo = Scr(1)-1.0d0
    end if
    if (kHi < n) then
      eHi = (Scr(kHi)+Scr(kHi+1))/2.0d0
    else
      eHi = Scr(n)+1.0d0
    end if
    !write(6,'(a,i5)') 'kLo (beta)  . . . . . . .',kLo
    !write(6,'(a,i5)') 'kHi (beta)  . . . . . . .',kHi
    !write(6,'(a,f12.6)') 'eLo (beta)  . . . . . . .',eLo
    !write(6,'(a,f12.6)') 'eHi (beta)  . . . . . . .',eHi
    nOcc = kLo-1
    nAct = kHi-kLo+1
    OccNo = 1.0d0*(nBeta-nOcc)/nAct
  end if
  !write(6,'(a,f12.6)') 'OccNo (beta)  . . . . . .',OccNo
end if
do i=1,n
  if (Eps(i) < eLo) then
    Occ(i) = Occ(i)+1.0d0
  else if (Eps(i) < eHi) then
    Occ(i) = Occ(i)+OccNo
  end if
end do
!----------------------------------------------------------------------*
! Print population (debug)                                             *
!----------------------------------------------------------------------*
!write(6,*)
!write(6,'(a)') 'Occupation of orbitals'
!write(6,'(a)') '----------------------'
!write(6,*)
!write(6,'(a,i5,a,i5)') 'Printing',m,' out of',n
!write(6,*)
!write(6,'(i5,1h-,i5,2x,10f12.4)') (i,Min(i+9,m),(Occ(j),j=i,Min(i+9,m)),i=1,m,10)
!----------------------------------------------------------------------*
! Done!                                                                *
!----------------------------------------------------------------------*
  return

end subroutine GoPop
