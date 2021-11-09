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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine WriteInt(IntMat,TermMat,mMat,nMat,OccNumMat1,OccNumMat2,MatEl,ForceField,E1,E2,T0,harmfreq1,harmfreq2,x_anharm1, &
                    x_anharm2,l_IntMat_1,l_IntMat_2,l_TermMat_1,l_TermMat_2,nDimTot,nOsc)
!  Purpose:
!    Write vibrational levels, intensities to log.
!
!  Input:
!    FC       : Real*8 two dimensional array - Franck-Condon
!               factors.
!    Int_Mat  : Real*8 two dimensional array - Intensities.
!    Term_Mat : Real*8 two dimensional array - Vibronic levels.
!    mMat     : Integer two dimensional array - oscillator quanta for
!               ground state.
!    nMat     : Integer two dimensional array - oscillator quanta for
!               excited state.
!    m_plot,
!    n_plot   : Integer array - transitions wanted in output.
!
!  Uses
!    TabMod
!    VibMod
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, u6

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
#include "dims.fh"
#include "indims.fh"
#include "inputdata.fh"
parameter(nfreq=2000)
real*8 Intensity, max_Intensity
real*8 IntMat(0:l_IntMat_1,0:l_IntMat_2)
real*8 TermMat(0:l_TermMat_1,0:l_TermMat_2)
real*8 OccNumMat1(0:nDimTot-1,nOsc), OccNumMat2(0:nDimTot-1,nOsc)
integer mMat(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2)
character*23 mMatChar(0:mdim1+1)
character*23 nMatChar(0:ndim1+1)
character*23 CharTemp
integer TermMin, TermMax, ivee_cm, ivee_nm
real*8 E1(nDimTot), E2(nDimTot)
real*8 harmfreq1(nOsc), harmfreq2(nOsc)
real*8 x_anharm1(nOsc,nOsc), x_anharm2(nOsc,nOsc)
logical MatEl, ForceField
character*12 format
integer nvTabDim
#include "inout.fh"
#include "WrkSpc.fh"

write(u6,*)
write(u6,*)
write(u6,*)
write(u6,*)
write(u6,'(a27,a)') ' ',' ================================================='
write(u6,'(a27,a)') ' ','|                                                 |'
write(u6,'(a27,a)') ' ','|       Results from intensity calculations       |'
write(u6,'(a27,a)') ' ','|                                                 |'
write(u6,'(a27,a)') ' ',' ================================================='
write(u6,*)
write(u6,*)

! Dimensions.
max_mOrd = l_IntMat_1
max_nOrd = l_IntMat_2
nvar = nOsc

if (nvar == 1) format = '(a1,  i3,a2)'
if (nvar == 2) format = '(a1, 2i3,a2)'
if (nvar == 3) format = '(a1, 3i3,a2)'
if (nvar == 4) format = '(a1, 4i3,a2)'
if (nvar == 5) format = '(a1, 5i3,a2)'
if (nvar == 6) format = '(a1, 6i2,a2)'
if (nvar == 7) format = '(a1, 7i2,a2)'
if (nvar == 8) format = '(a1, 8i1,a2)'
if (nvar == 9) format = '(a1, 9i1,a2)'
if (nvar == 10) format = '(a1,10i1,a2)'
if (nvar == 11) format = '(a1,11i1,a2)'
if (nvar == 12) format = '(a1,12i1,a2)'
if (nvar == 13) format = '(a1,13i1,a2)'
if (nvar == 14) format = '(a1,14i1,a2)'
if (nvar == 15) format = '(a1,15i1,a2)'
if (nvar == 16) format = '(a1,16i1,a2)'
if (nvar == 17) format = '(a1,17i1,a2)'
if (nvar == 18) format = '(a1,18i1,a2)'
if (nvar == 19) format = '(a1,19i1,a2)'
if (nvar == 20) format = '(a1,20i1,a2)'
if (nvar < 21) then
  do i=0,mdim1
    write(CharTemp,fmt=format) '(',(mMat(i,j),j=1,nvar),')'
    mMatChar(i) = CharTemp
  end do
  do i=0,ndim1
    write(CharTemp,fmt=format) '(',(nMat(i,j),j=1,nvar),')'
    nMatChar(i) = CharTemp
  end do
else
  do i=0,max_mOrd
    write(CharTemp,'(i5)') i
    mMatChar(i) = CharTemp
  end do
  do i=0,max_nOrd
    write(CharTemp,'(i5)') i
    nMatChar(i) = CharTemp
  end do
  !vv
  iprintLevel = 0
  if (iprintLevel == 10) then
    write(u6,*)
    write(u6,*)
    write(u6,*) ' ','The meaning of n and nprime in FC-table below:'
    write(u6,fmt='(a1,a,40a3,a)') ' ','=============',('===',i=1,nvar),'='
    !write(u6,fmt='(40a3)',advance='no')  ('===',i=1,nvar)
    !write(u6,fmt='(   a)',advance='yes') '='
    write(u6,*)
    write(u6,*) ' ','     n                                  Oscillator quanta'
    write(u6,fmt='(a1,a,40a3,a)') ' ','-------------',('---',i=1,nvar),'-'
    !write(u6,fmt='(40a3)',advance='no')  ('---',i=1,nvar)
    !write(u6,fmt='(   a)',advance='yes') '-'
    if (max_mOrd > max_nOrd) then
      do i=0,max_mOrd
        write(u6,'(a4,i4,a5,40i3)') ' ',i,' ',(mMat(i,j),j=1,nvar)
      end do
    else
      do i=0,max_nOrd
        write(u6,'(a4,i4,a5,40i3)') ' ',i,' ',(nMat(i,j),j=1,nvar)
      end do
    end if
    write(u6,fmt='(a1,a,40a3,a)') ' ','=============',('===',i=1,nvar),'='
    !write(u6,fmt='(40a3)',advance='no')  ('===',i=1,nvar)
    !write(u6,fmt='(   a)',advance='yes') '='
  end if
end if

! Write vibrational levels.
if (WriteVibLevels) then
  call GetMem('VibLevel1','Allo','Real',ipVibLevel1,max_mOrd+1)
  call GetMem('VibLevel2','Allo','Real',ipVibLevel2,max_nOrd+1)
  if (MatEl) then
    k = 0
    do iOrd=1,max_mOrd+1
      Work(ipVibLevel1+k) = (E1(iOrd)-E1(1))
      Work(ipVibLevel2+k) = T0+(E2(iOrd)-E2(1))
      k = k+1
    end do
  else
    call GetMem('level1','Allo','Inte',iplevel1,nvar)
    call GetMem('level2','Allo','Inte',iplevel2,nvar)

    G1 = Zero
    G2 = G1+T0
    k = 0
    Work(ipVibLevel1+k) = Zero
    do iv=1,nvar
      iWork(iplevel1+iv-1) = mMat(0,iv)
    end do
    if (max_mOrd > 0) then
      do iOrd=1,max_mOrd
        k = k+1
        do iv=1,nvar
          iWork(iplevel2+iv-1) = mMat(iOrd,iv)
        end do
        l_harm = nOsc
        call TransEnergy(G1,x_anharm1,harmfreq1,iWork(iplevel1),G1,x_anharm1,harmfreq1,iWork(iplevel2),Work(ipVibLevel1+k),l_harm)
      end do
    end if
    k = 0
    do iOrd=0,max_nOrd
      do iv=1,nvar
        iWork(iplevel2+iv-1) = nMat(iOrd,iv)
      end do
      l_harm = nOsc
      call TransEnergy(G1,x_anharm2,harmfreq2,iWork(iplevel1),G2,x_anharm2,harmfreq2,iWork(iplevel2),Work(ipVibLevel2+k),l_harm)
      k = k+1
    end do
    call GetMem('level1','Free','Inte',iplevel1,nvar)
    call GetMem('level2','Free','Inte',iplevel2,nvar)
  end if

  l_m = max_mOrd+1
  l_n = max_nOrd+1
  call GetMem('VibSort1','Allo','Inte',ipVibSort1,l_m*2)
  call GetMem('VibSort2','Allo','Inte',ipVibSort2,l_n*2)

  do iOrd=0,max_mOrd
    iWork(ipVibSort1+iOrd) = int(Work(ipVibLevel1+iOrd)*HarToRcm)
    iWork(ipVibSort1+iOrd+l_m) = iOrd
  end do
  do iOrd=0,max_nOrd
    iwork(ipVibSort2+iOrd) = int(Work(ipVibLevel2+iOrd)*HarToRcm)
    iWork(ipVibSort2+iOrd+l_n) = iOrd
  end do
  call GetMem('VibLevel1','Free','Real',ipVibLevel1,max_mOrd+1)
  call GetMem('VibLevel2','Free','Real',ipVibLevel2,max_nOrd+1)

  if (max_mOrd > 0) then
    do i=0,max_mOrd
      do j=i+1,max_mOrd
        if (iWork(ipVibSort1+j) < iWork(ipVibSort1+i)) then
          itemp = iWork(ipVibSort1+j)
          iwork(ipVibSort1+j) = iWork(ipVibSort1+i)
          iWork(ipVibSort1+i) = itemp
          itemp = iWork(ipVibSort1+j+l_m)
          iwork(ipVibSort1+j+l_m) = iWork(ipVibSort1+i+l_m)
          iwork(ipVibSort1+i+l_m) = itemp
        end if
      end do
    end do
  end if
  if (max_nOrd > 0) then
    do i=0,max_nOrd
      do j=i+1,max_nOrd
        if (iwork(ipVibSort2+j) < iWork(ipVibSort2+i)) then
          itemp = iwork(ipVibSort2+j)
          iWork(ipVibSort2+j) = iwork(ipVibSort2+i)
          iWork(ipVibSort2+i) = itemp
          itemp = iWork(ipVibSort2+j+l_n)
          iWork(ipVibSort2+j+l_n) = iWork(ipVibSort2+i+l_n)
          iWork(ipVibSort2+i+l_n) = itemp
        end if
      end do
    end do
  end if

  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Vibrational levels for ground state'
  write(u6,*) ' ','=============================================='
  write(u6,*)
  write(u6,*) ' ','     Level                     Energy(cm-1)'
  write(u6,*) ' ','----------------------------------------------'
  do i=0,max_mOrd
    iOrd = iwork(ipVibSort1+i+l_m)
    write(u6,'(a4,a15,a15,i7)') ' ',mMatChar(iOrd),' ',iWork(ipVibSort1+i)
  end do
  write(u6,*) ' ','=============================================='

  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Vibrational levels for excited state'
  write(u6,*) ' ','=============================================='
  write(u6,*)
  write(u6,*) ' ','     Level                     Energy(cm-1)'
  write(u6,*) ' ','----------------------------------------------'
  do i=0,max_nOrd
    iOrd = iWork(ipVibSort2+i+l_n)
    write(u6,'(a4,a15,a15,i7)') ' ',nMatChar(iOrd),' ',iWork(ipVibSort2+i)
  end do
  write(u6,*) ' ','=============================================='
  call GetMem('VibSort1','Free','Inte',ipVibSort1,l_m*2)
  call GetMem('VibSort2','Free','Inte',ipVibSort2,l_n*2)

end if
!!
!!---- Find start and end points for a given total quanta in mMat and nMat.
m_plot_max = l_m_plot
n_plot_max = l_n_plot
max_mQuanta = iWork(ipm_plot)
if (m_plot_max > 1) then
  do i=2,m_plot_max
    if (iWork(ipm_plot+i-1) > max_mQuanta) then
      max_mQuanta = iWork(ipm_plot+i-1)
    end if
  end do
end if
max_nQuanta = iWork(ipn_plot)
if (n_plot_max > 1) then
  do i=2,n_plot_max
    if (iWork(ipn_plot+i-1) > max_nQuanta) then
      max_nQuanta = iWork(ipn_plot+i-1)
    end if
  end do
end if
maxQuanta = max(max_mQuanta,max_nQuanta)
call GetMem('mMatStart','Allo','Inte',ipmMatStart,maxQuanta+1)
call GetMem('mMatStop','Allo','Inte',ipmMatStop,maxQuanta+1)
iWork(ipmMatStart) = 0
iWork(ipmMatStop) = 0
do i=1,maxQuanta
  iWork(ipmMatStart+i) = iWork(ipmMatStop+i-1)+1
  call TabDim2_drv(i,nvar,nvTabDim)
  iWork(ipmMatStop+i) = nvTabDim-1
end do

! Intensities.
n = (l_TermMat_1+2)*(l_TermMat_2+2)
l_TermSort = n+1
call GetMem('TermSort','Allo','Inte',ipTermSort,l_TermSort*3)
nval = 0
max_Intensity = Zero

if (MatEl .and. (.not. ForceField)) then
  if (max_mQuanta <= max_nQuanta) then
    do i=1,m_plot_max
      do iOrd=iWork(ipmMatStart+iwork(ipm_plot+i-1)),iWork(ipmMatStop+iWork(ipm_plot+i-1))
        do jOrd=0,3*max_mOrd/4
          iWork(ipTermSort+nval) = int(TermMat(iOrd,jOrd)*HarToRcm)
          iWork(ipTermSort+nval+l_TermSort) = iOrd
          iWork(ipTermSort+nval+l_TermSort*2) = jOrd
          nval = nval+1
        end do
      end do
    end do
  else
    do iOrd=0,3*max_mOrd/4
      do j=1,n_plot_max
        do jOrd=iWork(ipmMatStart+iWork(ipn_plot+j-1)),iWork(ipmMatStop+iWork(ipn_plot+j-1))
          if (IntMat(iOrd,jOrd) > max_Intensity) max_Intensity = IntMat(iOrd,jOrd)
          iWork(ipTermSort+nval) = int(abs(TermMat(iOrd,jOrd))*HarToRcm)
          iWork(ipTermSort+nval+l_TermSort) = iOrd
          iWork(ipTermSort+nval+l_TermSort*2) = jOrd
          nval = nval+1
        end do
      end do
    end do
  end if
else
  do i=1,m_plot_max
    do j=1,n_plot_max
      do iOrd=iWork(ipmMatStart+iWork(ipm_plot+i-1)),iWork(ipmMatStop+iWork(ipm_plot+i-1))
        do jOrd=iwork(ipmMatStart+iWork(ipn_plot+j-1)),iWork(ipmMatStop+iwork(ipn_plot+j-1))
          if (IntMat(iOrd,jOrd) > max_Intensity) max_Intensity = IntMat(iOrd,jOrd)
        end do
      end do
    end do
  end do
  do i=1,m_plot_max
    do j=1,n_plot_max
      do iOrd=iWork(ipmMatStart+iWork(ipm_plot+i-1)),iWork(ipmMatStop+iWork(ipm_plot+i-1))
        do jOrd=iWork(ipmMatStart+iWork(ipn_plot+j-1)),iWork(ipmMatStop+iWork(ipn_plot+j-1))
          !if ((IntMat(iOrd,jOrd) > int_thrs*max_Intensity) .or. &
          !if ((IntMat(iOrd,jOrd) > 1.0e-3_wp*max_Intensity) .or. ((iOrd == 0 ) .and. (jOrd == 0))) then
          iWork(ipTermSort+nval) = int(abs(TermMat(iOrd,jOrd))*HarToRcm)
          iWork(ipTermSort+nval+l_TermSort) = iOrd
          iWork(ipTermSort+nval+l_TermSort*2) = jOrd
          nval = nval+1
          !end if
        end do
      end do
    end do
  end do
end if

nval = nval-1
do i=0,nval
  do j=i+1,nval
    if (iWork(ipTermSort+j) < iWork(ipTermSort+i)) then
      itemp = iWork(ipTermSort+j)
      iWork(ipTermSort+j) = iWork(ipTermSort+i)
      iWork(ipTermSort+i) = itemp
      itemp = iWork(ipTermSort+j+l_TermSort)
      iWork(ipTermSort+j+l_TermSort) = iWork(ipTermSort+i+l_TermSort)
      iWork(ipTermSort+i+l_TermSort) = itemp
      itemp = iWork(ipTermSort+j+l_TermSort*2)
      iWork(ipTermSort+j+l_TermSort*2) = iWork(ipTermSort+i+l_TermSort*2)
      iWork(ipTermSort+i+l_TermSort*2) = itemp
    end if
  end do
end do

! MatEl,ForceField,m_plot_max,n_plot_max== F T 1 2/5

write(u6,*)
write(u6,'(90A)') ' ',('=',iv=1,89)
if (OscStr) then
  write(u6,'(A)') '  Ground                        Excited                         Energy         Oscillator'
  write(u6,'(A)') '  State                         State                      cm-1 /  eV / nm      Strength'
else
  write(u6,'(A)') '  Ground                        Excited                         Energy        Intensities'
  write(u6,'(A)') '  State                         State                      cm-1 /  eV / nm'
end if
write(u6,'(90A)') ' ',('-',iv=1,89)
const = One

do i=0,nval
  iOrd = iWork(ipTermSort+i+l_TermSort)
  jOrd = iWork(ipTermSort+i+l_TermSort*2)
  vee_cm = iWork(ipTermSort+i)
  vee_nm = 1.0e7_wp/vee_cm
  vee_eV = vee_cm*auToeV/auTocm
  vee = vee_eV
  if (Use_cm) vee = vee_cm
  if (Use_nm) vee = vee_nm
  ivee_cm = int(vee_cm+Half)
  ivee_nm = int(vee_nm+Half)
  Intensity = IntMat(iOrd,jOrd)/const
  if ((MatEl) .and. (.not. ForceField)) then
    if (m_plot_max < n_plot_max) then
      write(u6,'(a1,a15,a5,a1,f5.2,a1,f5.2,a1,f5.2,a1,a5,i7,a10,e12.3)') ' ',mMatChar(iOrd),'---> ','(',OccNumMat2(jOrd,1),',', &
                                                                         OccNumMat2(jOrd,2),',',OccNumMat2(jOrd,3),')',' ', &
                                                                         ivee_cm,' ',Intensity
      if (Intensity > 1.0e-6_wp) then
        call Add_Info('Energy',[vee_cm],1,5)
        call Add_Info('Intensity',[Intensity],1,5)
      end if
    else
      write(u6,'(a1,a15,a5,a1,f5.2,a1,f5.2,a1,f5.2,a1,a5,i7,a10,e12.3)') ' ',mMatChar(iOrd),'<--- ','(',OccNumMat2(jOrd,1),',', &
                                                                         OccNumMat2(jOrd,2),',',OccNumMat2(jOrd,3),')',' ', &
                                                                         ivee_cm,' ',Intensity
      if (Intensity > 1.0e-6_wp) then
        call Add_Info('Energy',[vee_cm],1,5)
        call Add_Info('Intensity',[Intensity],1,5)
      end if
    end if
  else
    if (m_plot_max < n_plot_max) then
      write(u6,'(a1,a23,a7,a23,a4,i6,a1,f5.2,a1,i4,a2,e12.3)') ' ',mMatChar(iOrd),' --->  ',nMatChar(jOrd),' ',ivee_cm,'/',vee_eV, &
                                                               '/',ivee_nm,'  ',Intensity
      if (Intensity > 1.0e-6_wp) then
        call Add_Info('Energy',[vee],1,5)
        call Add_Info('Intensity',[Intensity],1,5)
      end if
    else
      write(u6,'(a1,a23,a7,a23,a4,i6,a1,f5.2,a1,i4,a2,e12.3)') ' ',mMatChar(iOrd),' <---  ',nMatChar(jOrd),' ',ivee_cm,'/',vee_eV, &
                                                               '/',ivee_nm,'  ',Intensity
      if (Intensity > 1.0e-6_wp) then
        call Add_Info('Energy',[vee],1,5)
        call Add_Info('Intensity',[Intensity],1,5)
      end if
    end if
  end if
end do
write(u6,'(90A)') ' ',('=',iv=1,89)
write(u6,*)
write(u6,*)

!***********************************************************************

! Write to plot file.

conv = Zero
if ((.not. Use_nm) .and. (.not. Use_cm)) conv = 1.239842e-4_wp
if (Use_cm) conv = One
if (Use_nm) conv = 1.0e7_wp

! TermMin/Max must be in cm-1
if (plotwindow) then
  if (Use_nm) then                         ! nm -> cm-1
    TermMin = int(conv/cmend+0.999999_wp)   ! Note the inversion
    TermMax = int(conv/cmstart-0.999999_wp) ! cmstart <=> cmend
  else if (Use_cm) then                    ! already cm-1
    TermMin = int(cmstart)
    TermMax = int(cmend+0.999999_wp)
  else                                     ! eV -> cm-1
    TermMin = int(cmstart/conv+0.999999_wp)
    TermMax = int(cmend/conv-0.999999_wp)
  end if
else
  TermMin = int(iWork(ipTermSort))-nfreq
  TermMax = int(iWork(ipTermSort+nval))+nfreq
end if
l_plotvec = TermMax-TermMin+1
call GetMem('plotvec','Allo','Real',ipplotvec,l_plotvec)
call dcopy_(l_plotvec,[Zero],0,Work(ipplotvec),1)

call molcas_open(plotunit,'plot.intensity')
!open(unit=plotUnit,file='plot.intensity')
if (broadplot) then
  FWHM = hbarcm/(Two*LifeTime)
  G2 = FWHM*Half
  do iTrans=0,nval
    iOrd = iWork(ipTermSort+iTrans+l_TermSort)
    jOrd = iWork(ipTermSort+iTrans+l_TermSort*2)
    ifreq = iWork(ipTermSort+iTrans)
    do ipoint=TermMin,TermMax
      Work(ipplotVec+ipoint-TermMin) = Work(ipplotVec+ipoint-TermMin)+IntMat(iOrd,jOrd)*(G2**2)/((G2**2)+(ipoint-ifreq)**2)
    end do
  end do
  do ipoint=TermMin,TermMax
    if (Work(ipplotVec+ipoint-TermMin) > Zero) then
      if (.not. Use_nm) then
        write(plotUnit,'(f12.6,e15.6)') ipoint*conv,Work(ipplotVec+ipoint-TermMin)
      else
        write(plotUnit,'(f12.6,e15.6)') conv/ipoint,Work(ipplotVec+ipoint-TermMin)
      end if
    end if
  end do
else
  do iTrans=0,nval
    iOrd = iWork(ipTermSort+iTrans+l_TermSort)
    jOrd = iWork(ipTermSort+iTrans+l_TermSort*2)
    ifreq = iWork(ipTermSort+iTrans)
    if ((TermMin <= ifreq) .and. (ifreq <= TermMax)) then
      Work(ipplotVec+ifreq-TermMin) = Work(ipplotVec+ifreq-TermMin)+IntMat(iOrd,jOrd)
    end if
  end do
  do ipoint=TermMin,TermMax
    if (Work(ipplotVec+ipoint-TermMin) > Zero) then
      if (.not. Use_nm) then
        write(plotUnit,'(f12.6,e15.6)') ipoint*conv,Zero
        write(plotUnit,'(f12.6,e15.6)') ipoint*conv,Work(ipplotVec+ipoint-TermMin)
        write(plotUnit,'(f12.6,e15.6)') ipoint*conv,Zero
      else
        write(plotUnit,'(f12.6,e15.6)') conv/ipoint,Zero
        write(plotUnit,'(f12.6,e15.6)') conv/ipoint,Work(ipplotVec+ipoint-TermMin)
        write(plotUnit,'(f12.6,e15.6)') conv/ipoint,Zero
      end if
    end if
  end do
  !call GetMem('plotvec','Free','Real',ipplotvec,l_plotvec)
end if
call GetMem('plotvec','Free','Real',ipplotvec,l_plotvec)
close(plotUnit)

call GetMem('TermSort','Free','Inte',ipTermSort,l_TermSort*3)
call GetMem('mMatStart','Free','Inte',ipmMatStart,maxQuanta+1)
call GetMem('mMatStop','Free','Inte',ipmMatStop,maxQuanta+1)

! Avoid unused argument warnings
if (.false.) call Unused_real_array(OccNumMat1)

end subroutine WriteInt
