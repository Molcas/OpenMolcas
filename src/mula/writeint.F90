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

subroutine WriteInt(IntMat,TermMat,mMat,nMat,OccNumMat2,MatEl,ForceField,E1,E2,T0,harmfreq1,harmfreq2,x_anharm1,x_anharm2, &
                    l_IntMat_1,l_IntMat_2,l_TermMat_1,l_TermMat_2,nDimTot,nOsc)
!  Purpose:
!    Write vibrational levels, intensities to log.
!
!  Input:
!    FC       : Real two dimensional array - Franck-Condon factors.
!    Int_Mat  : Real two dimensional array - Intensities.
!    Term_Mat : Real two dimensional array - Vibronic levels.
!    mMat     : Integer two dimensional array - oscillator quanta for ground state.
!    nMat     : Integer two dimensional array - oscillator quanta for excited state.
!    m_plot,
!    n_plot   : Integer array - transitions wanted in output.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use mula_global, only: broadplot, cmstart, cmend, hbarcm, LifeTime, m_plot, mdim1, mdim2, n_plot, ndim1, ndim2, OscStr, &
                       plotwindow, Use_cm, Use_nm, WriteVibLevels
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, auTocm, auToeV
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mMat(0:mdim1,mdim2), nMat(0:ndim1,ndim2), l_IntMat_1, l_IntMat_2, l_TermMat_1, l_TermMat_2, &
                                 nDimTot, nOsc
real(kind=wp), intent(in) :: IntMat(0:l_IntMat_1,0:l_IntMat_2), TermMat(0:l_TermMat_1,0:l_TermMat_2), &
                             OccNumMat2(0:nDimTot-1,nOsc), E1(nDimTot), E2(nDimTot), T0, harmfreq1(nOsc), harmfreq2(nOsc), &
                             x_anharm1(nOsc,nOsc), x_anharm2(nOsc,nOsc)
logical(kind=iwp), intent(in) :: MatEl, ForceField
integer(kind=iwp) :: i, ifreq, iOrd, ipoint, iprintLevel, itemp(3), iTrans, iv, ivee_cm, ivee_nm, j, jOrd, k, l_harm, m_plot_max, &
                     max_mOrd, max_mQuanta, max_nOrd, max_nQuanta, maxQuanta, n, n_plot_max, nval, nvar, plotunit, TermMax, &
                     TermMin, nvTabDim
real(kind=wp) :: const, conv, FWHM, G1, G2, Intensity, max_Intensity, vee, vee_cm, vee_eV, vee_nm
character(len=23) :: CharTemp
character(len=12) :: frmt
integer(kind=iwp), allocatable :: level1(:), level2(:), mMatStart(:), mMatStop(:), TermSort(:,:), VibSort1(:,:), VibSort2(:,:)
real(kind=wp), allocatable :: plotvec(:), VibLevel1(:), VibLevel2(:)
character(len=23), allocatable :: mMatChar(:), nMatChar(:)
integer(kind=iwp), parameter :: nfreq = 2000
integer(kind=iwp), external :: isfreeunit

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

call mma_allocate(mMatChar,[0,mdim1+1],label='mMatChar')
call mma_allocate(nMatChar,[0,ndim1+1],label='nMatChar')

if (nvar < 21) then
  if (nvar > 0) write(frmt,'(a,i2,a)') '(a1,',nvar,'i1,a2)'
  do i=0,mdim1
    write(CharTemp,fmt=frmt) '(',(mMat(i,j),j=1,nvar),')'
    mMatChar(i) = CharTemp
  end do
  do i=0,ndim1
    write(CharTemp,fmt=frmt) '(',(nMat(i,j),j=1,nvar),')'
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
  call mma_allocate(VibLevel1,[0,max_mOrd],label='VibLevel1')
  call mma_allocate(VibLevel2,[0,max_nOrd],label='VibLevel2')
  if (MatEl) then
    k = 0
    do iOrd=1,max_mOrd+1
      VibLevel1(k) = E1(iOrd)-E1(1)
      VibLevel2(k) = T0+(E2(iOrd)-E2(1))
      k = k+1
    end do
  else
    call mma_allocate(level1,nvar,label='level1')
    call mma_allocate(level2,nvar,label='level2')

    G1 = Zero
    G2 = G1+T0
    k = 0
    VibLevel1(k) = Zero
    level1(:) = mMat(0,:)
    if (max_mOrd > 0) then
      do iOrd=1,max_mOrd
        k = k+1
        level2(:) = mMat(iOrd,:)
        l_harm = nOsc
        call TransEnergy(G1,x_anharm1,harmfreq1,level1,G1,x_anharm1,harmfreq1,level2,VibLevel1(k),l_harm)
      end do
    end if
    k = 0
    do iOrd=0,max_nOrd
      level2(:) = nMat(iOrd,:)
      l_harm = nOsc
      call TransEnergy(G1,x_anharm2,harmfreq2,level1,G2,x_anharm2,harmfreq2,level2,VibLevel2(k),l_harm)
      k = k+1
    end do
    call mma_deallocate(level1)
    call mma_deallocate(level2)
  end if

  call mma_allocate(VibSort1,[1,2],[0,max_mOrd],label='VibSort1')
  call mma_allocate(VibSort2,[1,2],[0,max_nOrd],label='VibSort2')

  do iOrd=0,max_mOrd
    VibSort1(1,iOrd) = int(VibLevel1(iOrd)*auTocm)
    VibSort1(2,iOrd) = iOrd
  end do
  do iOrd=0,max_nOrd
    VibSort2(1,iOrd) = int(VibLevel2(iOrd)*auTocm)
    VibSort2(2,iOrd) = iOrd
  end do
  call mma_deallocate(VibLevel1)
  call mma_deallocate(VibLevel2)

  if (max_mOrd > 0) then
    do i=0,max_mOrd
      do j=i+1,max_mOrd
        if (VibSort1(1,j) < VibSort1(1,i)) then
          itemp(1:2) = VibSort1(:,j)
          VibSort1(:,j) = VibSort1(:,i)
          VibSort1(:,i) = itemp(1:2)
        end if
      end do
    end do
  end if
  if (max_nOrd > 0) then
    do i=0,max_nOrd
      do j=i+1,max_nOrd
        if (VibSort2(1,j) < VibSort2(1,i)) then
          itemp(1:2) = VibSort2(:,j)
          VibSort2(:,j) = VibSort2(:,i)
          VibSort2(:,i) = itemp(1:2)
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
    iOrd = VibSort1(2,i)
    write(u6,'(a4,a15,a15,i7)') ' ',mMatChar(iOrd),' ',VibSort1(1,i)
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
    iOrd = VibSort2(2,i)
    write(u6,'(a4,a15,a15,i7)') ' ',nMatChar(iOrd),' ',VibSort2(1,i)
  end do
  write(u6,*) ' ','=============================================='
  call mma_deallocate(VibSort1)
  call mma_deallocate(VibSort2)

end if

! Find start and end points for a given total quanta in mMat and nMat.
m_plot_max = size(m_plot)
n_plot_max = size(n_plot)
max_mQuanta = maxval(m_plot)
max_nQuanta = maxval(n_plot)
maxQuanta = max(max_mQuanta,max_nQuanta)
call mma_allocate(mMatStart,[0,maxQuanta],label='mMatStart')
call mma_allocate(mMatStop,[0,maxQuanta],label='mMatStop')
mMatStart(0) = 0
mMatStop(0) = 0
do i=1,maxQuanta
  mMatStart(i) = mMatStop(i-1)+1
  call TabDim(i,nvar,nvTabDim)
  mMatStop(i) = nvTabDim-1
end do

! Intensities.
n = (l_TermMat_1+2)*(l_TermMat_2+2)
call mma_allocate(TermSort,[1,3],[0,n],label='TermSort')
nval = 0
max_Intensity = Zero

if (MatEl .and. (.not. ForceField)) then
  if (max_mQuanta <= max_nQuanta) then
    do i=1,m_plot_max
      do iOrd=mMatStart(m_plot(i)),mMatStop(m_plot(i))
        do jOrd=0,3*max_mOrd/4
          TermSort(1,nval) = int(TermMat(iOrd,jOrd)*auTocm)
          TermSort(2,nval) = iOrd
          TermSort(3,nval) = jOrd
          nval = nval+1
        end do
      end do
    end do
  else
    do iOrd=0,3*max_mOrd/4
      do j=1,n_plot_max
        do jOrd=mMatStart(n_plot(j)),mMatStop(n_plot(j))
          if (IntMat(iOrd,jOrd) > max_Intensity) max_Intensity = IntMat(iOrd,jOrd)
          TermSort(1,nval) = int(abs(TermMat(iOrd,jOrd))*auTocm)
          TermSort(2,nval) = iOrd
          TermSort(3,nval) = jOrd
          nval = nval+1
        end do
      end do
    end do
  end if
else
  do i=1,m_plot_max
    do j=1,n_plot_max
      do iOrd=mMatStart(m_plot(i)),mMatStop(m_plot(i))
        do jOrd=mMatStart(n_plot(j)),mMatStop(n_plot(j))
          if (IntMat(iOrd,jOrd) > max_Intensity) max_Intensity = IntMat(iOrd,jOrd)
        end do
      end do
    end do
  end do
  do i=1,m_plot_max
    do j=1,n_plot_max
      do iOrd=mMatStart(m_plot(i)),mMatStop(m_plot(i))
        do jOrd=mMatStart(n_plot(j)),mMatStop(n_plot(j))
          !if ((IntMat(iOrd,jOrd) > int_thrs*max_Intensity) .or. &
          !if ((IntMat(iOrd,jOrd) > 1.0e-3_wp*max_Intensity) .or. ((iOrd == 0 ) .and. (jOrd == 0))) then
          TermSort(1,nval) = int(abs(TermMat(iOrd,jOrd))*auTocm)
          TermSort(2,nval) = iOrd
          TermSort(3,nval) = jOrd
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
    if (TermSort(1,j) < TermSort(1,i)) then
      itemp(:) = TermSort(:,j)
      TermSort(:,j) = TermSort(:,i)
      TermSort(:,i) = itemp
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
  iOrd = TermSort(2,i)
  jOrd = TermSort(3,i)
  vee_cm = TermSort(1,i)
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

call mma_deallocate(mMatChar)
call mma_deallocate(nMatChar)

!***********************************************************************

! Write to plot file.

conv = Zero
if ((.not. Use_nm) .and. (.not. Use_cm)) conv = 1.239842e-4_wp
if (Use_cm) conv = One
if (Use_nm) conv = 1.0e7_wp

! TermMin/Max must be in cm-1
if (plotwindow) then
  if (Use_nm) then                          ! nm -> cm-1
    TermMin = int(conv/cmend+0.999999_wp)   ! Note the inversion
    TermMax = int(conv/cmstart-0.999999_wp) ! cmstart <=> cmend
  else if (Use_cm) then                     ! already cm-1
    TermMin = int(cmstart)
    TermMax = int(cmend+0.999999_wp)
  else                                      ! eV -> cm-1
    TermMin = int(cmstart/conv+0.999999_wp)
    TermMax = int(cmend/conv-0.999999_wp)
  end if
else
  TermMin = int(TermSort(1,0))-nfreq
  TermMax = int(TermSort(1,nval))+nfreq
end if
call mma_allocate(plotvec,[TermMin,TermMax],label='plotvec')
plotvec(:) = Zero

plotUnit = isfreeunit(20)
call molcas_open(plotUnit,'plot.intensity')
!open(plotUnit,'plot.intensity')
if (broadplot) then
  FWHM = hbarcm/(Two*LifeTime)
  G2 = FWHM*Half
  do iTrans=0,nval
    iOrd = TermSort(2,iTrans)
    jOrd = TermSort(3,iTrans)
    ifreq = TermSort(1,iTrans)
    do ipoint=TermMin,TermMax
      plotVec(ipoint) = plotVec(ipoint)+IntMat(iOrd,jOrd)*(G2**2)/((G2**2)+(ipoint-ifreq)**2)
    end do
  end do
  do ipoint=TermMin,TermMax
    if (plotVec(ipoint) > Zero) then
      if (.not. Use_nm) then
        write(plotUnit,'(f12.6,e15.6)') ipoint*conv,plotVec(ipoint)
      else
        write(plotUnit,'(f12.6,e15.6)') conv/ipoint,plotVec(ipoint)
      end if
    end if
  end do
else
  do iTrans=0,nval
    iOrd = TermSort(2,iTrans)
    jOrd = TermSort(3,iTrans)
    ifreq = TermSort(1,iTrans)
    if ((TermMin <= ifreq) .and. (ifreq <= TermMax)) then
      plotVec(ifreq) = plotVec(ifreq)+IntMat(iOrd,jOrd)
    end if
  end do
  do ipoint=TermMin,TermMax
    if (plotVec(ipoint) > Zero) then
      if (.not. Use_nm) then
        write(plotUnit,'(f12.6,e15.6)') ipoint*conv,Zero
        write(plotUnit,'(f12.6,e15.6)') ipoint*conv,plotVec(ipoint)
        write(plotUnit,'(f12.6,e15.6)') ipoint*conv,Zero
      else
        write(plotUnit,'(f12.6,e15.6)') conv/ipoint,Zero
        write(plotUnit,'(f12.6,e15.6)') conv/ipoint,plotVec(ipoint)
        write(plotUnit,'(f12.6,e15.6)') conv/ipoint,Zero
      end if
    end if
  end do
end if
call mma_deallocate(plotvec)
close(plotUnit)

call mma_deallocate(TermSort)
call mma_deallocate(mMatStart)
call mma_deallocate(mMatStop)

end subroutine WriteInt
