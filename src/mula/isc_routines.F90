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
! Copyright (C) 2008,2009, Giovanni Ghigo                              *
!***********************************************************************

! Contains:
!   ISC_Rho
!   ISC_Ene
!   ISC_Rate
!   ISC_FCval

!  InterSystem Crossing rate evaluation: "Engine" routines
!  Author: Giovanni Ghigo
!          Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
!          28 Dec-08 - 06 Jan-09

subroutine ISC_Rho(iPrint,nOsc,new_n_max,dRho,energy1,energy2,minQ,dMinWind0,nMaxQ,harmfreq1,harmfreq2)
! Calculate State Density  dRho  GG 30-Dec-08
! Formula (86) taken from  M. Bixon, J. Jortner  JCP,48,715 (1969)

use Constants, only: Zero, One, Two, Six, Twelve, Half, Pi, auTocm, auToeV
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, nOsc
integer(kind=iwp), intent(out) :: new_n_max, minQ, nMaxQ(nOsc)
real(kind=wp), intent(out) :: dRho
real(kind=wp), intent(in) :: energy1, energy2, dMinWind0, harmfreq1(nOsc), harmfreq2(nOsc)
integer(kind=iwp) :: i, iOsc, jOsc
real(kind=wp) :: avFreq, avFreqSq, dAlpha, dBeta, dDn, dEtha, dFE, dLambda, dMaxFreq2, dMinFreq2, dMinWind, dZPE1, dZPE2, GE1, &
                 GE2, T0

if (iPrint >= 2) then
  write(u6,*)
  write(u6,*) ' State Density data:'
  write(u6,*) ' ============================================'
end if

dMinWind = dMinWind0
if (dMinWind == Zero) dMinWind = One
minQ = 0
dRho = Two/Pi
dRho = sqrt(dRho*nOsc)
dRho = dRho*(One-One/(Twelve*nOsc))
avFreq = Zero
avFreqSq = Zero
dLambda = One
dZPE1 = Zero
dZPE2 = Zero
dMinFreq2 = 1.0e99_wp
dMaxFreq2 = Zero
do iOsc=1,nOsc
  avFreq = avFreq+harmfreq2(iOsc)
  avFreqSq = avFreqSq+(harmfreq2(iOsc))**2
  dZPE1 = dZPE1+Half*harmfreq1(iOsc)
  dZPE2 = dZPE2+Half*harmfreq2(iOsc)
  dMinFreq2 = min(dMinFreq2,harmfreq2(iOsc))
  dMaxFreq2 = max(dMaxFreq2,harmfreq2(iOsc))
end do
GE1 = energy1+dZPE1
GE2 = energy2+dZPE2
T0 = abs(GE2-GE1)
avFreq = avFreq/nOsc
avFreqSq = avFreqSq/nOsc

do jOsc=1,nOsc
  dLambda = dLambda*harmfreq2(jOsc)/avFreq
end do
new_n_max = int(Half+(GE1-GE2)/dMinFreq2)
dAlpha = avFreqSq/(avFreq**2)
dBeta = ((nOsc-1)*(nOsc-2)*dAlpha-nOsc**2)/(Six*nOsc)
dEtha = abs(GE1-GE2)/dZPE2
dRho = dRho/avFreq
dRho = dRho*dAlpha
dRho = dRho/(1+dEtha)
dDn = (One+(Two/dEtha))
dDn = dDn**(dEtha/Two)
dDn = dDn*(One+(dEtha/Two))
dDn = dDn**nOsc
dRho = dRho*dDn
dFE = (One+dEtha)**2
dFE = (One-One/dFE)**dBeta
dRho = dRho*dFE

do jOsc=1,nOsc
  nMaxQ(jOsc) = int(Half+(T0+dMinWind/dRho)/harmfreq2(jOsc))
end do
minQ = int(Half+(T0-dMinWind/dRho)/dMaxFreq2)
if (minQ < 0) then
  write(u6,*)
  write(u6,*) ' ***** ERROR ******'
  write(u6,*) ' Window too large !'
  write(u6,*) ' ******************'
  call Quit_OnUserError()
end if

if (iPrint >= 2) then
  write(u6,'(a,f11.6,a)') '  T_0  = ',T0,' (au)'
  write(u6,'(a,f11.3,a)') '  T_0  = ',T0*auTocm,' (cm-1)'
  write(u6,'(a,f11.3,a)') '  T_0  = ',T0*auToeV,' (eV)'
  write(u6,'(a,d14.3,a)') '  State Density (dRho) = ',dRho,' (au-1)'
  write(u6,'(a,g14.3,a)') '  State Density (dRho) = ',dRho/auTocm,' (cm)'
  write(u6,'(a,g17.9,a)') '  1/dRho = ',auTocm/dRho,' (cm-1)'
  write(u6,'(a,f7.3,a)') '  Expansion factor =',dMinWind
  write(u6,'(a,g17.9,a)') '  Window = (+/-)',Half*dMinWind*auTocm/dRho,' (cm-1)'
end if
if (iPrint >= 3) then
  write(u6,*) ' Maximum quantum numbers:',(nMaxQ(i),i=1,nOsc)
  write(u6,*) ' Minimum quantum number: ',minQ
  write(u6,*) ' Suggested n_max (new_n_max)=',new_n_max
  write(u6,*)
end if
call XFlush(u6)

return

end subroutine ISC_Rho

subroutine ISC_Ene(iPrint,nOsc,max_nOrd,nYes,nMat,nTabDim,GE1,GE2,harmfreq1,harmfreq2,x_anharm1,x_anharm2,dMinWind,dRho,lVec)
! Calculate Energy of Levels  GG 30-Dec-08 - 08-Jan-09

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, auTocm
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, nOsc, max_nOrd, nTabDim, nMat(0:nTabDim,nOsc)
integer(kind=iwp), intent(inout) :: nYes, lVec(0:nTabDim)
real(kind=wp), intent(in) :: GE1, GE2, harmfreq1(nOsc), harmfreq2(nOsc), x_anharm1(nOsc,nOsc), x_anharm2(nOsc,nOsc), dRho
real(kind=wp), intent(inout) :: dMinWind
integer(kind=iwp) :: i, iOrd, j, l_harm, loc_n_max, nYes_start
real(kind=wp) :: dEne, dWlow, dWup
logical(kind=iwp) :: lUpdate
integer(kind=iwp), allocatable :: level1(:), level2(:), lTVec(:)
real(kind=wp), allocatable :: EneMat(:)

if (dMinWind == Zero) then
  lUpDate = .true.
  dMinWind = One
else
  lUpDate = .false.
end if
call mma_allocate(lTVec,[0,max_nOrd],label='lTVec')
lTVec(:) = lVec(0:max_nOrd)

! Energy calculation

if (iPrint >= 4) then
  write(u6,*)
  write(u6,*) ' States in the preliminar window :'
  if (nOsc <= 24) then
    write(u6,'(a,108a)') '  ',('=',i=1,108)
    write(u6,*) '     jOrd    ene/au    ene/cm-1 Vibrational quantum numbers'
    write(u6,'(a,108a)') '  ',('-',i=1,108)
  else
    write(u6,'(a,36a)') '  ',('=',i=1,36)
    write(u6,*) '        #    jOrd   ene/au      ene/cm-1'
    write(u6,'(a,36a)') '  ',('-',i=1,36)
  end if
  call XFlush(u6)
end if

call mma_allocate(level1,nOsc,label='level1')
call mma_allocate(level2,nOsc,label='level2')
call mma_allocate(EneMat,[0,max_nOrd],label='EneMat')
level1(:) = 0
do iOrd=0,max_nOrd
  if (lVec(iOrd) == 1) then
    level2(:) = nMat(iOrd,:)
    l_harm = nOsc
    call TransEnergy(GE1,x_anharm1,harmfreq1,level1,GE2,x_anharm2,harmfreq2,level2,dEne,l_harm)
    EneMat(iOrd) = dEne
    if (iPrint >= 4) then
      if (nOsc <= 24) then
        loc_n_max = 0
        do j=1,nOsc
          loc_n_max = loc_n_max+nMat(iOrd,j)
        end do
        write(u6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,dEne*auTocm,loc_n_max,': ',(nMat(iOrd,j),j=1,nOsc)
      else
        write(u6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,dEne*auTocm,loc_n_max
      end if
    end if
  end if
end do

! Energy selection

if (iPrint >= 3) then
  write(u6,*)
  write(u6,*) ' States in the window :'
  if (nOsc <= 24) then
    write(u6,'(a,108a)') '  ',('=',i=1,108)
    write(u6,*) '     jOrd    ene/au    ene/cm-1 Vibrational quantum numbers'
    write(u6,'(a,108a)') '  ',('-',i=1,108)
  else
    write(u6,'(a,36a)') '  ',('=',i=1,36)
    write(u6,*) '        #    jOrd   ene/au      ene/cm-1'
    write(u6,'(a,36a)') '  ',('-',i=1,36)
  end if
  call XFlush(u6)
end if

nYes_start = nYes
do
  dWlow = Half*dMinWind/dRho
  dWup = dWlow
  do iOrd=0,max_nOrd
    lVec(iOrd) = lTVec(iOrd)
    if (lVec(iOrd) == 1) then
      dEne = EneMat(iOrd)
      if ((dEne < -dWlow) .or. (dEne > dWup)) then
        lVec(iOrd) = 0
        nYes = nYes-1
      else
        lVec(iOrd) = 1
        if (iPrint >= 3) then
          if (nOsc <= 24) then
            loc_n_max = 0
            do j=1,nOsc
              loc_n_max = loc_n_max+nMat(iOrd,j)
            end do
            write(u6,'(a2,i8,f11.6,f11.4,i4,a2,24i3)') ' ',iOrd,dEne,dEne*auTocm,loc_n_max,': ',(nMat(iOrd,j),j=1,nOsc)
          else
            write(u6,'(a2,i8,f11.6,f11.4,i4        )') ' ',iOrd,dEne,dEne*auTocm,loc_n_max
          end if
        end if
      end if
    end if
  end do
  if ((nYes > 1) .or. (.not. lUpDate)) exit
  dMinWind = dMinWind+One
  nYes = nYes_start
end do

call mma_deallocate(lTVec)
call mma_deallocate(level1)
call mma_deallocate(level2)
call mma_deallocate(EneMat)

if (iPrint >= 3) then
  if (nOsc <= 30) write(u6,'(a,108a)') '  ',('-',i=1,108)
  if (nOsc > 30) write(u6,'(a,36a)') '  ',('-',i=1,36)
  write(u6,'(a,f12.9,a,f12.9,a)') '  Window: ',-dWlow,' / ',dWup,' (au)'
  write(u6,'(a,f12.6,a,f12.6,a)') '  Window: ',-dWlow*auTocm,' / ',dWup*auTocm,' (cm-1)'
end if
if (iPrint >= 2) then
  write(u6,*) ' Final number of States=',nYes
end if
if ((dMinWind > One) .and. lUpDate .and. (iPrint >= 1)) then
  write(u6,*)
  write(u6,*) ' *** Warning: Expansion factor has been set to ',dMinWind
  write(u6,*)
end if
call XFlush(u6)

return

end subroutine ISC_Ene

subroutine ISC_Rate(iPrint,nOsc,max_nOrd,iMaxYes,nYes,dMinWind,VibWind2,C1,C2,W2,det0,det1,det2,C,W,r01,r02,nTabDim,nMat,nInc, &
                    nDec,m_max,n_max,FC00,dRho)
! Estimate ISC rate  GG 30-Dec-08

use mula_global, only: hbarcm, TranDip
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Pi, auTocm
use Definitions, only: wp, iwp, u6, ItoB

implicit none
integer(kind=iwp), intent(inout) :: nTabDim
integer(kind=iwp), intent(in) :: iPrint, nOsc, iMaxYes, nYes, VibWind2(nYes), nMat(0:nTabDim,nOsc), nInc(0:iMaxYes,nOsc), &
                                 nDec(0:iMaxYes,nOsc), m_max, n_max
integer(kind=iwp), intent(out) :: max_nOrd
real(kind=wp), intent(in) :: dMinWind, C1(nOsc,nOsc), C2(nOsc,nosc), W2(nOsc,nOsc), det0, det1, det2, C(nOsc,nOsc), W(nOsc,nOsc), &
                             r01(nOsc), r02(nOsc), dRho
real(kind=wp), intent(out) :: FC00
integer(kind=iwp) :: ii, m_max_ord, max_nInc, mx_max_ord, n_max_ord, nvTabDim, nx_max_ord
real(kind=wp) :: const, dLT, dRate, dSoc, dSum
real(kind=wp), allocatable :: FCWind2(:)
integer(kind=iwp), parameter :: MB = 1048576

call TabDim(n_max,nOsc,nvTabDim)
max_nOrd = nvTabDim-1
call TabDim(m_max,nOsc,nvTabDim)
m_max_ord = nvTabDim-1
call TabDim(min(n_max,m_max+1),nOsc,nvTabDim)
mx_max_ord = nvTabDim-1
call TabDim(min(m_max,n_max+1),nOsc,nvTabDim)
nx_max_ord = nvTabDim-1
call TabDim(n_max-1,nOsc,nvTabDim)
max_nInc = nvTabDim-1
call TabDim(n_max,nOsc,nvTabDim)
n_max_ord = nvTabDim-1

mx_max_ord = 0 ! CGGn
if (iPrint >= 3) write(u6,*) ' Memory allocated for U matrix:',(n_max_ord+1)*(mx_max_ord+1),' words,  ', &
                             (n_max_ord+1)*(mx_max_ord+1)*ItoB/MB,' MB.'
call XFlush(u6)

call mma_allocate(FCWind2,nYes,label='FCWind2')
call ISC_FCval(iPrint,iMaxYes,nTabDim,C1,det1,r01,C2,W2,det2,r02,m_max_ord,n_max_ord,mx_max_ord,max_nInc,nx_max_ord,nMat,nInc, &
               nDec,C,W,det0,FC00,nOsc,nYes,VibWind2,FCWind2)

const = Two*Pi/hbarcm
if (iPrint >= 4) then
  write(u6,*)
  write(u6,*) '  const =',const
  write(u6,*) '  dRho/cm =',dRho/auTocm
  write(u6,*) '  const*dRho=',const*dRho/auTocm
end if
const = const*dRho/auTocm

dSum = Zero
do ii=1,nYes
  dSum = dSum+FCWind2(ii)**2
end do
call mma_deallocate(FCWind2)

dSoc = Zero
do ii=1,3
  dSOC = dSOC+TranDip(ii)**2
end do

dRate = const*dSum*dSoc/dMinWind
dLT = One/dRate

if (iPrint >= 3) then
  write(u6,*) '  Sum of squares of FC factors =',dSum
  write(u6,*) '  Root-square of the sum =',sqrt(dSum)
  write(u6,*) '  dSOC =',dSOC
end if

if (iPrint >= 1) then
  write(u6,*)
  write(u6,*) ' InterSystem Crossing rate constant:'
  write(u6,*) ' ===================================='
  write(u6,'(a,e10.2,a)') '  ISC Rate Constant  ',dRate,' sec-1'
  write(u6,'(a,e10.2,a)') '  Lifetime           ',dLT,' sec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' msec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' microsec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' nsec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' psec'
  dLT = dLT*1.0e3_wp
  if ((dLT > One) .and. (dLT <= 1.0e3_wp)) write(u6,'(a19,f5.1,a)') ' ',dLT,' fsec'
  write(u6,*) ' ------------------------------------'
  write(u6,*)
  write(u6,*)
  call XFlush(u6)
end if

return

end subroutine ISC_Rate

!  InterSystem Crossing rate evaluation: Multidimensional Franck-Condon
!  Modified copy of FCval by Giovanni Ghigo.
!  Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
!  28-Dec-08 - 06-Jan-09 ; June 2009

subroutine ISC_FCval(iPrint,iMaxYes,nTabDim,C1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,max_nInc,max_nInc2,nMat,nInc, &
                     nDec,C,W,det0,FC00,nOsc,nYes,VibWind2,FCWind2)

use mula_global, only: ndim1, ndim2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iPrint, iMaxYes, max_mOrd, max_nOrd, max_nOrd2, max_nInc2, nOsc, nMat(0:ndim1,ndim2), &
                                 nInc(0:iMaxYes,nOsc), nDec(0:iMaxYes,nOsc), nYes, VibWind2(nYes)
integer(kind=iwp), intent(out) :: nTabDim
real(kind=wp), intent(in) :: C1(nOsc,nOsc), det1, r01(nOsc), C2(nOsc,nOsc), W2(nOsc,nOsc), det2, r02(nOsc), C(nOsc,nOsc), &
                             W(nOsc,nOsc), det0
integer(kind=iwp), intent(inout) :: max_nInc
real(kind=wp), intent(out) :: FC00, FCWind2(nYes)
integer(kind=iwp) :: i, ii, iOrd, j, jOrd, kOsc, kOsc_start, loc_n_max, lOsc, n, nMaxMat
real(kind=wp) :: const, det, dFC, FC00_exp
real(kind=wp), allocatable :: A2(:,:), A2B2T(:,:), Alpha(:,:), Alpha1(:,:), Alpha2(:,:), B2(:,:), Beta(:,:), d2(:), L(:,:), &
                              r_temp1(:), r_temp2(:), sqr(:), temp(:,:), temp1(:,:), temp2(:,:), U(:,:)
real(kind=wp), external :: Ddot_

!write(u6,*) 'CGGt[ISC_FCval] Enter'
!write(u6,*) '                nYes = ',nYes
!write(u6,*) '                VibWind2 :',(VibWind2(i),i=1,nYes)
!write(u6,*) '            L matrix:',max_mOrd,max_nInc2
!write(u6,*) '            U matrix:',max_nOrd,max_nOrd2
!do i=0,iMaxYes
!  write(u6,*) (nInc(i,j),j=1,nOsc)
!end do
!call XFlush(u6)

! Initialize.
nMaxMat = max(max_mOrd+1,max_nOrd+1)
!write(u6,*) '            nMaxMat=',nMaxMat
nTabDim = max(nMaxMat,8)
!write(u6,*) '            nTabDim=',nTabDim
call mma_allocate(temp,nOsc,nOsc,label='temp')
call mma_allocate(temp1,nOsc,nOsc,label='temp1')
call mma_allocate(temp2,nOsc,nOsc,label='temp2')

! Setup sqr table.
n = nTabDim+1
call mma_allocate(sqr,[0,n],label='sqr')
do i=0,nTabDim+1
  sqr(i) = sqrt(real(i,kind=wp))
end do

! Calculate alpha1, alpha2 and alpha.
!write(u6,*) 'CGGt[FCVal] Calculate alpha(s)'
!call XFlush(u6)
call mma_allocate(Alpha,nOsc,nOsc,label='Alpha')
call mma_allocate(Alpha1,nOsc,nOsc,label='Alpha1')
call mma_allocate(Alpha2,nOsc,nOsc,label='Alpha2')
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C1,nOsc,C1,nOsc,Zero,Alpha1,nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOsc,Half,C2,nOsc,C2,nOsc,Zero,Alpha2,nOsc)
temp(:,:) = Alpha1+Alpha2
Alpha(:,:) = Half*temp

!call xxDgemul(C,nOsc,'T',C,nOsc,'N',alpha,nOsc,nOsc,nOsc,nOsc)
!alpha(:,:) = Half*alpha

! Calculate C using a Cholesky factorization of 2*alpha.
!call Cholesky(temp,C)

! Calculate W.
!call unitmat(W,nOsc)
!temp(:,:) = C
!call Dool_MULA(temp,W,det0)

! Calculate r00.
call mma_allocate(r_temp1,nOsc,label='r_temp1')
call mma_allocate(r_temp2,nOsc,label='r_temp2')

! Calculate beta.
!write(u6,*) 'CGGt[FCVal] Calculate beta.'
!call XFlush(u6)
call mma_allocate(Beta,nOsc,nOsc,label='Beta')
do i=1,nOsc
  do j=1,nOsc
    temp1(j,i) = C1(i,j)
  end do
  !write(u6,*) 'CGGt C1(',i,',j)=',(C1(i,jj),jj=1,nOsc)
end do
!call XFlush(u6)
temp1(:,:) = Alpha1
temp(:,:) = Two*Alpha

call Dool_MULA(temp,nOsc,nOsc,temp1,nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,Alpha2,nOsc,temp1,nOsc,Zero,Beta,nOsc)

call mma_deallocate(Alpha1)
call mma_deallocate(Alpha2)

! Calculate FC00.
!r_temp1(:) = r02-r01
r_temp1(:) = r01-r02

call DGEMM_('N','N',nOsc,1,nOsc,One,Beta,nOsc,r_temp1,nOsc,Zero,r_temp2,nOsc)
FC00_exp = Ddot_(nOsc,r_temp1,1,r_temp2,1)
FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)
!write(u6,*) 'CGGt[FCVal] FC00_exp,FC00=',FC00_exp,FC00
!call XFlush(u6)

call mma_deallocate(Beta)

! Calculate A, B and d matrices.
call mma_allocate(A2,nOsc,nOsc,label='A2')
call mma_allocate(B2,nOsc,nOsc,label='B2')
call mma_allocate(d2,nOsc,label='d2')

call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C2,nOsc,W,nOsc,Zero,A2,nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W2,nOsc,C,nOsc,Zero,temp,nOsc)
B2(:,:) = A2-temp

const = -sqr(8)
call DGEMM_('T','N',nOsc,1,nOsc,const,W2,nOsc,r_temp2,nOsc,Zero,d2,nOsc)

! Calculate A1B1T and A2B2T.
call mma_allocate(A2B2T,nOsc,nOsc,label='A2B2T')

call DGEMM_('N','T',nOsc,nOsc,nOsc,One,A2,nOsc,B2,nOsc,Zero,A2B2T,nOsc)

! Initialize L matrix.
call mma_allocate(L,[0,max_mOrd],[0,max_nInc2],label='L')
L(:,:) = Zero
L(0,0) = One

! If max_mOrd > 0 then set up L(m,0).
if (max_mOrd > 0) then
  write(u6,*) '*****************************************'
  write(u6,*) ' Hot initial state not implemented yet !'
  write(u6,*) '*****************************************'
  write(u6,*)
  call Quit_OnUserError()
end if

! Initialize U matrix.
!write(u6,*) 'CGGt[FCVal] Initialize U matrix.'
!call XFlush(u6)
call mma_allocate(U,[0,max_nOrd],[0,max_nOrd2],label='U')
U(:,:) = Zero
U(0,0) = One
!GGt -------------------------------------------------------------------
!do kOsc=1,nOsc
!  write(u6,*) 'CGGt[FCVal] d2(..)=',d2(kOsc)
!  call XFlush(u6)
!end do
!GGt -------------------------------------------------------------------

! If max_nOrd > 0 then set up U(n,0).
!write(u6,*) 'CGGt[FCVal] max_nOrd > 0 then set up U(n,0).'
!call XFlush(u6)
do kOsc=1,nOsc
  !write(u6,*) 'CGGt[FCVal] nInc(0,kOsc)=',nInc(0,kOsc)
  !call XFlush(u6)
  !write(u6,*) 'CGGt[FCVal] d2(..)=',d2(kOsc)
  !call XFlush(u6)
  U(nInc(0,kOsc),0) = d2(kOsc)
end do
!write(u6,*) 'CGGt[FCVal] max_nInc=',max_nInc
!write(u6,*) 'CGGt[FCVal]  iMaxYes=',iMaxYes
!call XFlush(u6)
max_nInc = min(max_nInc,iMaxYes)
if (max_nInc > 0) then
  do iOrd=1,max_nInc
    !write(u6,*) '              iOrd=',iOrd
    !call XFlush(u6)
    kOsc_start = nOsc
    do while ((nMat(iOrd,kOsc_start) == 0) .and. (kOsc_start > 1))
      kOsc_start = kOsc_start-1
    end do
    do kOsc=kOsc_start,nOsc
      do lOsc=1,nOsc
        !write(u6,*) 'iOrd,kOsc_start,kOsc,lOsc==',iOrd,kOsc_start,kOsc,lOsc
        !call XFlush(u6)
        if (nMat(iOrd,lOsc) > 0) then

          U(nInc(iOrd,kOsc),0) = U(nInc(iOrd,kOsc),0)+sqr(nMat(iOrd,lOsc))*A2B2T(kOsc,lOsc)*U(nDec(iOrd,lOsc),0)
        end if
      end do

      U(nInc(iOrd,kOsc),0) = (U(nInc(iOrd,kOsc),0)+d2(kOsc)*U(iOrd,0))/sqr(nMat(nInc(iOrd,kOsc),kOsc))
    end do
  end do
end if

! Use recursion formula to obtain the rest of U.
!write(u6,*) '            Use recursion ... rest of U.'
!write(u6,*) '            max_nOrd2=',max_nOrd2
!call XFlush(u6)
do jOrd=1,max_nOrd2
  lOsc = nOsc
  do while ((nMat(jOrd,lOsc) == 0) .and. (lOsc > 1))
    lOsc = lOsc-1
  end do
  do iOrd=0,max_nOrd
    do kOsc=1,nOsc
      if (nMat(iOrd,kOsc) > 0) then
        !write(u6,*) '              ',iOrd,kOsc,nMat(iOrd,kOsc)
        U(iOrd,jOrd) = U(iOrd,jOrd)+sqr(nMat(iOrd,kOsc))/sqr(nMat(jOrd,lOsc))*A2(kOsc,lOsc)*U(nDec(iOrd,kOsc),nDec(jOrd,lOsc))
      end if
    end do
  end do
end do

call mma_deallocate(A2)
call mma_deallocate(B2)
call mma_deallocate(d2)
call mma_deallocate(A2B2T)

! Calculate Franck-Condon factors.
if (iPrint >= 3) then
  write(u6,*) ' Franck-Condon factors for States in the Window:'
  write(u6,'(a,36a)') '  ',('=',i=1,36)
  write(u6,*) '     #     jOrd   FC factor     jSum'
  write(u6,'(a,36a)') '  ',('-',i=1,36)
end if
do ii=1,nYes
  jOrd = VibWind2(ii)
  dFC = FC00*L(0,0)*U(jOrd,0)
  FCWind2(ii) = dFC
  if (iPrint >= 3) then
    loc_n_max = 0
    do j=1,nOsc
      loc_n_max = loc_n_max+nMat(jOrd,j)
    end do
    write(u6,'(a2,i5,i9,e15.6,a2,i4)') ' ',ii,jOrd,FCWind2(ii),' ',loc_n_max
  end if
end do
if (iPrint >= 3) then
  write(u6,'(a,36a)') '  ',('-',i=1,36)
  write(u6,*) ' FC_00 =',FC00
  write(u6,*)
end if

if (iPrint >= 4) then
  write(u6,*)
  write(u6,*) ' Full Franck-Condon factors (FC_00=',FC00,'):'
  write(u6,*) ' =================================================='
  write(u6,*) '    jOrd   FC            level'
  write(u6,*) ' --------------------------------------------------'
  do jOrd=0,max_nOrd
    loc_n_max = 0
    do j=1,nOsc
      loc_n_max = loc_n_max+nMat(jOrd,j)
    end do
    dFC = FC00*L(0,0)*U(jOrd,0)
    write(u6,'(a,i8,e15.6,a2,i4,a2,24i3)') ' ',jOrd,dFC,' ',loc_n_max,' ',(nMat(jOrd,j),j=1,nOsc)
  end do
  write(u6,*) ' --------------------------------------------------'
end if
!GGt -------------------------------------------------------------------

call mma_deallocate(L)
call mma_deallocate(U)
call mma_deallocate(Alpha)
call mma_deallocate(sqr)
call mma_deallocate(temp)
call mma_deallocate(temp1)
call mma_deallocate(temp2)
call mma_deallocate(r_temp1)
call mma_deallocate(r_temp2)

!write(u6,*) 'CGGt[FCVal] Exit'
!call XFlush(u6)

end subroutine ISC_FCval
