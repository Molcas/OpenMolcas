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

subroutine XT_dMoverdH(exch,nLoc,nCenter,nneq,neqv,neq,nss,nexch,nTempMagn,nT,NM,iopt,mem,Tmin,Tmax,chit_exp,eso,w,T,RROT,zJ, &
                       Xfield,THRS,XT_no_field,dipso,s_so,dipexch,s_exch,tinput,smagn,m_paranoid,m_accurate,doplot)
! this Subroutine computes the XT as M/H as it is observed for most of experiments.
! The M is averaged over the grid for each temperature point.
!       chi*t ----------- the units are cgsemu: [ cm^3*k/mol ]
!
!  WEX0, WEX1 WEX2 : Zeeman exchange energies
!  WL0, WL1, WL2   : Zeeman local energies
! Zeeman local reduced energies, using only NEXCH states;
!  WR0, WR1, WR2
! local statistical sum, Boltzmann distribution
!  ZL0, ZL1, ZL2
! local statistical sum, Boltzmann distribution, using only NEXCH states
!  ZR0, ZR1, ZR2
! spin magnetisation, from the local sites, using ALL states
!  SL0, SL1, SL2
! spin magnetisation, from the local sites, using only NEXCH states
!  SR0, SR1, SR2
! magnetisation, from local sites, using ALL states
!  ML0, ML1, ML2
! magnetisation, from local sites, using only NEXCH states
!  MR0, MR1, MR2
! spin magnetisation, from the exchange block
!  SEX0, SEX1, SEX2
! magnetisation, form the exchange block
!  MEX0, MEX1, MEX2
! exchange statistical sum, Boltzmann distribution
!  ZEX0, ZEX1, ZEX2
! total vectors in general coordinate system:
!  ZRT0, ZLT0, ZRT1, ZLT1, ZRT2, ZLT2, MRT0, MLT0, SRT0, SLT0, MRT1, MLT1, SRT1, SLT1, MRT2, MLT2, SRT2, SLT2
! data for total system:
! total statistical sum, Boltzmann distribution
!  MT0 : total magnetisation
!  MT1 : total magnetisation
!  ST0 : total spin magnetisation,
!  ST1 : total spin magnetisation,
!  ST2 : total spin magnetisation,
!  MT2 : total magnetisation
!  ZT0, ZT1, ZT2
! standard deviation data:
!  dev

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Three, Ten, mBohr, rNAVO
use Definitions, only: wp, iwp, u6, RtoB

implicit none
integer(kind=iwp), intent(in) :: exch, nLoc, nCenter, nneq, neqv, neq(nneq), nss(nneq), nexch(nneq), nTempMagn, nT, NM, iopt, mem
real(kind=wp), intent(in) :: Tmin, Tmax, chit_exp(nT), eso(nneq,nLoc), W(exch), T(nT+nTempMagn), RROT(nneq,neqv,3,3), zJ, Xfield, &
                             THRS, XT_no_field(nT+nTempMagn)
complex(kind=wp), intent(in) :: dipso(nneq,3,nLoc,nLoc), s_so(nneq,3,nLoc,nLoc), dipexch(3,exch,exch), s_exch(3,exch,exch)
logical(kind=iwp), intent(in) :: tinput, smagn, m_paranoid, doplot
integer(kind=iwp) :: i, ibuf, ibuf1, ibuf3, iM, info, isite, iT, j, jT, mem_local, n, nDirX, nTempTotal
real(kind=wp) :: dHX(3), dHY(3), dHZ(3), dltXF, F1, F2, wt(3), Xfield_1, Xfield_2, zt(3,3)
logical(kind=iwp) :: m_accurate
character(len=50) :: label
real(kind=wp), allocatable :: ESO_TMP(:), MEX0(:,:), MEX1(:,:), MEX2(:,:), ML0(:,:,:), ML1(:,:,:), ML2(:,:,:), MLT0(:,:,:), &
                              MLT1(:,:,:), MLT2(:,:,:), MR0(:,:,:), MR1(:,:,:), MR2(:,:,:), MRT0(:,:,:), MRT1(:,:,:), MRT2(:,:,:), &
                              MT0(:,:), MT1(:,:), MT2(:,:), SEX0(:,:), SEX1(:,:), SEX2(:,:), SL0(:,:,:), SL1(:,:,:), SL2(:,:,:), &
                              SLT0(:,:,:), SLT1(:,:,:), SLT2(:,:,:), SR0(:,:,:), SR1(:,:,:), SR2(:,:,:), SRT0(:,:,:), SRT1(:,:,:), &
                              SRT2(:,:,:), ST0(:,:), ST1(:,:), ST2(:,:), WEX0(:), WEX1(:), WEX2(:), WL0(:,:), WL1(:,:), WL2(:,:), &
                              WR0(:,:), WR1(:,:), WR2(:,:), XTM_dMdH(:), XTM_MH(:), XTtens_dMdH(:,:,:), XTtens_MH(:,:,:), ZEX0(:), &
                              ZEX1(:), ZEX2(:), ZL0(:,:), ZL1(:,:), ZL2(:,:), ZLT0(:,:), ZLT1(:,:), ZLT2(:,:), ZR0(:,:), ZR1(:,:), &
                              ZR2(:,:), ZRT0(:,:), ZRT1(:,:), ZRT2(:,:), ZT0(:), ZT1(:), ZT2(:)
complex(kind=wp), allocatable :: dipso_tmp(:,:,:), s_so_tmp(:,:,:)
real(kind=wp), parameter :: cm3tomB = rNAVO*mBohr/Ten ! in cm3 * mol-1 * T
#ifdef _DEBUGPRINT_
#  define _DBG_ .true.
integer(kind=iwp) :: l
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: DBG = _DBG_
real(kind=wp), external :: dev, dnrm2_

#ifndef _DEBUGPRINT_
#include "macros.fh"
unused_var(mem)
#endif

nTempTotal = nT+nTempMagn

!m_paranoid = .true. ! .false.
!ccc-------------------------------------------------------cccc
#ifdef _DEBUGPRINT_
write(u6,'(A,   I5)') '      exch =',exch
write(u6,'(A,   I5)') '      nLoc =',nLoc
write(u6,'(A,   I5)') '   nCenter =',nCenter
write(u6,'(A,   I5)') '      nneq =',nneq
write(u6,'(A,   I5)') '      neqv =',neqv
write(u6,'(A, 10I5)') '    neq(i) =',(neq(i),i=1,nneq)
write(u6,'(A, 10I5)') '  nexch(i) =',(nexch(i),i=1,nneq)
write(u6,'(A, 10I5)') '    nss(i) =',(nss(i),i=1,nneq)
write(u6,'(A,   I5)') ' nTempMagn =',nTempMagn
write(u6,'(A,   I5)') '        nT =',nT
write(u6,'(A,   I5)') '        nM =',nM
write(u6,'(A,   I5)') '      iopt =',iopt
write(u6,'(A,F12.5)') '      Tmin =',Tmin
write(u6,'(A,F12.5)') '      Tmax =',Tmax
write(u6,'(A,F12.5)') '        zJ =',zJ
write(u6,'(A,F12.5)') '    Xfield =',Xfield
write(u6,'(A,ES12.5)') '      THRS =',THRS
write(u6,*) 'm_accurate =',m_accurate
write(u6,*) '     smagn =',smagn
write(u6,*) 'm_paranoid =',m_paranoid
write(u6,*) '    tinput =',tinput
if (tinput) then
  do iT=1,nTempMagn
    write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',iT,')=',T(iT)
  end do
  do iT=1,nT
    jT = iT+nTempMagn
    write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',jT,')=',T(jT),' chiT_exp(',iT,')=',chit_exp(iT)
  end do

else
  do iT=1,nTempTotal
    write(u6,'(2(A,i3,A,F12.6,2x))') 'T(',iT,')=',T(iT)
  end do
end if
write(u6,'(A)') 'local ESO:'
do i=1,nneq
  write(u6,'(A,i3)') 'site:',i
  write(u6,'(10F12.5)') (eso(i,j),j=1,nss(i))
end do
write(u6,'(A)') 'EXCHANGE ENERGY'
write(u6,'(10F12.5)') (W(i),i=1,exch)
write(u6,'(A)') 'Rotation Matrices:'
do i=1,nneq
  write(u6,'(A,i3)') 'site:',i
  do j=1,neq(i)
    do l=1,3
      write(u6,'(10F12.5)') (RROT(i,j,l,n),n=1,3)
    end do
  end do
end do
call xFlush(u6)
#endif
!ccc-------------------------------------------------------cccc
write(u6,*)
write(u6,'(A)') repeat('%',95)
write(u6,'(16X,A)') 'CALCULATION OF THE FIELD-DEPENDENT MAGNETIC SUSCEPTIBILITY'
write(u6,'(18X,A)') 'within true (dM/dH) and "experimentalists" (M/H) models'
write(u6,'(A)') repeat('%',95)
write(u6,*)
write(u6,'(2x,A,F10.6,A)') 'Magnetic field strength:',Xfield,' tesla.'
if (tinput) then
  write(u6,'(2x,a)') 'Temperature dependence of the magnetic susceptibility and'
  write(u6,'(2x,a)') 'high-field magnetization will be calculated according to '
  write(u6,'(2x,a)') 'experimental values provided by the user in file "chitexp.input".'
else
  write(u6,'(2x,a,i3,a)') 'Temperature dependence of the magnetic susceptibility will be calculated in',nT,' points, '
  write(u6,'(2x,a,f4.1,a,f6.1,a)') 'equally distributed in temperature range ',tmin,' ---',tmax,' K.'
end if

!=======================================================================
mem_local = 0
call mma_allocate(WEX0,nM,'WEX0')
call mma_allocate(WEX1,nM,'WEX1')
call mma_allocate(WEX2,nM,'WEX2')
mem_local = mem_local+size(WEX0)*RtoB
mem_local = mem_local+size(WEX1)*RtoB
mem_local = mem_local+size(WEX2)*RtoB

call mma_allocate(WL0,nLoc,nneq,'WL0')
call mma_allocate(WL1,nLoc,nneq,'WL1')
call mma_allocate(WL2,nLoc,nneq,'WL2')
call mma_allocate(WR0,nLoc,nneq,'WR0')
call mma_allocate(WR1,nLoc,nneq,'WR1')
call mma_allocate(WR2,nLoc,nneq,'WR2')
mem_local = mem_local+size(WL0)*RtoB
mem_local = mem_local+size(WL1)*RtoB
mem_local = mem_local+size(WL2)*RtoB
mem_local = mem_local+size(WR0)*RtoB
mem_local = mem_local+size(WR1)*RtoB
mem_local = mem_local+size(WR2)*RtoB

call mma_allocate(ZL0,nTempTotal,nneq,'ZL0')
call mma_allocate(ZR0,nTempTotal,nneq,'ZR0')
call mma_allocate(ZL1,nTempTotal,nneq,'ZL1')
call mma_allocate(ZR1,nTempTotal,nneq,'ZR1')
call mma_allocate(ZL2,nTempTotal,nneq,'ZL2')
call mma_allocate(ZR2,nTempTotal,nneq,'ZR2')
mem_local = mem_local+size(ZL0)*RtoB
mem_local = mem_local+size(ZR0)*RtoB
mem_local = mem_local+size(ZL1)*RtoB
mem_local = mem_local+size(ZR1)*RtoB
mem_local = mem_local+size(ZL2)*RtoB
mem_local = mem_local+size(ZR2)*RtoB

call mma_allocate(SL0,3,nTempTotal,nneq,'SL0')
call mma_allocate(SR0,3,nTempTotal,nneq,'SR0')
call mma_allocate(SL1,3,nTempTotal,nneq,'SL1')
call mma_allocate(SR1,3,nTempTotal,nneq,'SR1')
call mma_allocate(SL2,3,nTempTotal,nneq,'SL2')
call mma_allocate(SR2,3,nTempTotal,nneq,'SR2')
call mma_allocate(ML0,3,nTempTotal,nneq,'ML0')
call mma_allocate(MR0,3,nTempTotal,nneq,'MR0')
call mma_allocate(ML1,3,nTempTotal,nneq,'ML1')
call mma_allocate(MR1,3,nTempTotal,nneq,'MR1')
call mma_allocate(ML2,3,nTempTotal,nneq,'ML2')
call mma_allocate(MR2,3,nTempTotal,nneq,'MR2')
mem_local = mem_local+size(SL0)*RtoB
mem_local = mem_local+size(SR0)*RtoB
mem_local = mem_local+size(SL1)*RtoB
mem_local = mem_local+size(SR1)*RtoB
mem_local = mem_local+size(SL2)*RtoB
mem_local = mem_local+size(SR2)*RtoB
mem_local = mem_local+size(ML0)*RtoB
mem_local = mem_local+size(MR0)*RtoB
mem_local = mem_local+size(ML1)*RtoB
mem_local = mem_local+size(MR1)*RtoB
mem_local = mem_local+size(ML2)*RtoB
mem_local = mem_local+size(MR2)*RtoB

call mma_allocate(ZEX0,nTempTotal,'ZEX0')
call mma_allocate(ZEX1,nTempTotal,'ZEX1')
call mma_allocate(ZEX2,nTempTotal,'ZEX2')
mem_local = mem_local+size(ZEX0)*RtoB
mem_local = mem_local+size(ZEX1)*RtoB
mem_local = mem_local+size(ZEX2)*RtoB

call mma_allocate(ZT0,nTempTotal,'ZT0')
call mma_allocate(ZT1,nTempTotal,'ZT1')
call mma_allocate(ZT2,nTempTotal,'ZT2')
ZT0(:) = Zero
ZT1(:) = Zero
ZT2(:) = Zero
mem_local = mem_local+size(ZT0)*RtoB
mem_local = mem_local+size(ZT1)*RtoB
mem_local = mem_local+size(ZT2)*RtoB

call mma_allocate(SEX0,3,nTempTotal,'SEX0')
call mma_allocate(SEX1,3,nTempTotal,'SEX1')
call mma_allocate(SEX2,3,nTempTotal,'SEX2')
call mma_allocate(MEX0,3,nTempTotal,'MEX0')
call mma_allocate(MEX1,3,nTempTotal,'MEX1')
call mma_allocate(MEX2,3,nTempTotal,'MEX2')
mem_local = mem_local+size(SEX0)*RtoB
mem_local = mem_local+size(SEX1)*RtoB
mem_local = mem_local+size(SEX2)*RtoB
mem_local = mem_local+size(MEX0)*RtoB
mem_local = mem_local+size(MEX1)*RtoB
mem_local = mem_local+size(MEX2)*RtoB

call mma_allocate(ST0,3,nTempTotal,'ST0')
call mma_allocate(ST1,3,nTempTotal,'ST1')
call mma_allocate(ST2,3,nTempTotal,'ST2')
call mma_allocate(MT0,3,nTempTotal,'MT0')
call mma_allocate(MT1,3,nTempTotal,'MT1')
call mma_allocate(MT2,3,nTempTotal,'MT2')
ST0(:,:) = Zero
ST1(:,:) = Zero
ST2(:,:) = Zero
MT0(:,:) = Zero
MT1(:,:) = Zero
ST2(:,:) = Zero
mem_local = mem_local+size(ST0)*RtoB
mem_local = mem_local+size(ST1)*RtoB
mem_local = mem_local+size(ST2)*RtoB
mem_local = mem_local+size(MT0)*RtoB
mem_local = mem_local+size(MT1)*RtoB
mem_local = mem_local+size(MT2)*RtoB

call mma_allocate(XTM_MH,nTempTotal,'XTM_MH')
call mma_allocate(XTM_dMdH,nTempTotal,'XTM_dMdH')
call mma_allocate(XTtens_MH,3,3,nTempTotal,'XTtens_MH')
call mma_allocate(XTtens_dMdH,3,3,nTempTotal,'XTtens_dMdH')
mem_local = mem_local+size(XTM_MH)*RtoB
mem_local = mem_local+size(XTM_dMdH)*RtoB
mem_local = mem_local+size(XTtens_MH)*RtoB
mem_local = mem_local+size(XTtens_dMdH)*RtoB

call mma_allocate(ZRT0,nCenter,nTempTotal,'ZRT0')
call mma_allocate(ZLT0,nCenter,nTempTotal,'ZLT0')
call mma_allocate(ZRT1,nCenter,nTempTotal,'ZRT1')
call mma_allocate(ZLT1,nCenter,nTempTotal,'ZLT1')
call mma_allocate(ZRT2,nCenter,nTempTotal,'ZRT2')
call mma_allocate(ZLT2,nCenter,nTempTotal,'ZLT2')
mem_local = mem_local+size(ZRT0)*RtoB
mem_local = mem_local+size(ZLT0)*RtoB
mem_local = mem_local+size(ZRT1)*RtoB
mem_local = mem_local+size(ZLT1)*RtoB
mem_local = mem_local+size(ZRT2)*RtoB
mem_local = mem_local+size(ZLT2)*RtoB

call mma_allocate(MRT0,nCenter,3,nTempTotal,'MRT0')
call mma_allocate(MLT0,nCenter,3,nTempTotal,'MLT0')
call mma_allocate(SRT0,nCenter,3,nTempTotal,'SRT0')
call mma_allocate(SLT0,nCenter,3,nTempTotal,'SLT0')
call mma_allocate(MRT1,nCenter,3,nTempTotal,'MRT1')
call mma_allocate(MLT1,nCenter,3,nTempTotal,'MLT1')
call mma_allocate(SRT1,nCenter,3,nTempTotal,'SRT1')
call mma_allocate(SLT1,nCenter,3,nTempTotal,'SLT1')
call mma_allocate(MRT2,nCenter,3,nTempTotal,'MRT2')
call mma_allocate(MLT2,nCenter,3,nTempTotal,'MLT2')
call mma_allocate(SRT2,nCenter,3,nTempTotal,'SRT2')
call mma_allocate(SLT2,nCenter,3,nTempTotal,'SLT2')
mem_local = mem_local+size(MRT0)*RtoB
mem_local = mem_local+size(MLT0)*RtoB
mem_local = mem_local+size(SRT0)*RtoB
mem_local = mem_local+size(SLT0)*RtoB
mem_local = mem_local+size(MRT1)*RtoB
mem_local = mem_local+size(MLT1)*RtoB
mem_local = mem_local+size(SRT1)*RtoB
mem_local = mem_local+size(SLT1)*RtoB
mem_local = mem_local+size(MRT2)*RtoB
mem_local = mem_local+size(MLT2)*RtoB
mem_local = mem_local+size(SRT2)*RtoB
mem_local = mem_local+size(SLT2)*RtoB

#ifdef _DEBUGPRINT_
write(u6,*) 'XT_MH:  memory allocated (local):'
write(u6,*) 'mem_local=',mem_local
write(u6,*) 'XT_MH:  memory allocated (total):'
write(u6,*) 'mem_total=',mem+mem_local
#endif

!=======================================================================
nDirX = 3

dHX(1) = One
dHY(1) = Zero
dHZ(1) = Zero

dHX(2) = Zero
dHY(2) = One
dHZ(2) = Zero

dHX(3) = Zero
dHY(3) = Zero
dHZ(3) = One

Xfield_1 = Xfield*0.99999_wp
Xfield_2 = Xfield*1.00001_wp
dltXF = Xfield_2-Xfield_1

! ///  opening the loop over different directions of the magnetic field
do iM=1,nDirX
  WL0(:,:) = Zero
  WL1(:,:) = Zero
  WL2(:,:) = Zero
  WR0(:,:) = Zero
  WR1(:,:) = Zero
  WR2(:,:) = Zero

  MRT0(:,:,:) = Zero
  MLT0(:,:,:) = Zero
  SRT0(:,:,:) = Zero
  SLT0(:,:,:) = Zero
  MRT1(:,:,:) = Zero
  MLT1(:,:,:) = Zero
  SRT1(:,:,:) = Zero
  SLT1(:,:,:) = Zero
  MRT2(:,:,:) = Zero
  MLT2(:,:,:) = Zero
  SRT2(:,:,:) = Zero
  SLT2(:,:,:) = Zero
  ! exchange magnetization:
  call MAGN(EXCH,NM,dHX(iM),dHY(iM),dHZ(iM),XField,W,zJ,THRS,DIPEXCH,S_EXCH,nTempTotal,T,smagn,Wex0,Zex0,Sex0,Mex0,m_paranoid,DBG)
  call MAGN(EXCH,NM,dHX(iM),dHY(iM),dHZ(iM),XField_1,W,zJ,THRS,DIPEXCH,S_EXCH,nTempTotal,T,smagn,Wex1,Zex1,Sex1,Mex1,m_paranoid,DBG)
  call MAGN(EXCH,NM,dHX(iM),dHY(iM),dHZ(iM),XField_2,W,zJ,THRS,DIPEXCH,S_EXCH,nTempTotal,T,smagn,Wex2,Zex2,Sex2,Mex2,m_paranoid,DBG)
# ifdef _DEBUGPRINT_
  do i=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex0(',i,')=',(Mex0(l,i),l=1,3)
  end do
  do i=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex1(',i,')=',(Mex1(l,i),l=1,3)
  end do
  do i=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex2(',i,')=',(Mex2(l,i),l=1,3)
  end do
  do i=1,nTempTotal
    write(u6,'(A,i3,A,9ES22.14)') 'Mex 2-1 (',i,')=',(Mex2(l,i)-Mex1(l,i),l=1,3)
  end do
# endif

  call Add_Info('dM/dH   Mex0',[dnrm2_(3*nTempTotal,Mex0,1)],1,8)
  call Add_Info('dM/dH   Sex0',[dnrm2_(3*nTempTotal,Sex0,1)],1,8)
  call Add_Info('dM/dH   Zex0',[dnrm2_(nTempTotal,Zex0,1)],1,8)
  call Add_Info('dM/dH   Wex0',[dnrm2_(nM,Wex0,1)],1,8)

  call Add_Info('dM/dH   Mex1',[dnrm2_(3*nTempTotal,Mex1,1)],1,8)
  call Add_Info('dM/dH   Sex1',[dnrm2_(3*nTempTotal,Sex1,1)],1,8)
  call Add_Info('dM/dH   Zex1',[dnrm2_(nTempTotal,Zex1,1)],1,8)
  call Add_Info('dM/dH   Wex1',[dnrm2_(nM,Wex1,1)],1,8)

  call Add_Info('dM/dH   Mex2',[dnrm2_(3*nTempTotal,Mex2,1)],1,8)
  call Add_Info('dM/dH   Sex2',[dnrm2_(3*nTempTotal,Sex2,1)],1,8)
  call Add_Info('dM/dH   Zex2',[dnrm2_(nTempTotal,Zex2,1)],1,8)
  call Add_Info('dM/dH   Wex2',[dnrm2_(nM,Wex2,1)],1,8)

  ! compute local magnetizations:
  if (m_accurate) then
    do i=1,nneq
      ! all states:
      if (NSS(i) > NEXCH(i)) then
        call mma_allocate(ESO_TMP,nss(i),label='ESO_TMP')
        call mma_allocate(dipso_tmp,3,nss(i),nss(i),label='dipso_tmp')
        call mma_allocate(s_so_tmp,3,nss(i),nss(i),label='s_so_tmp')
        ! this check is to avoid the unnecessary computation, in cases when no local excited states are present
        ESO_TMP(:) = ESO(i,1:NSS(i))
        dipso_tmp(:,:,:) = DIPSO(i,:,1:NSS(i),1:NSS(i))
        s_so_tmp(:,:,:) = S_SO(i,:,1:NSS(i),1:NSS(i))
        call MAGN(NSS(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField,ESO_TMP,zJ,THRS,dipso_tmp,s_so_tmp,nTempTotal,T,smagn, &
                  WL0(1:NSS(i),i),ZL0(:,i),SL0(:,:,i),ML0(:,:,i),m_paranoid,DBG)

        call MAGN(NSS(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField_1,ESO_TMP,zJ,THRS,dipso_tmp,s_so_tmp,nTempTotal,T,smagn, &
                  WL1(1:NSS(i),i),ZL1(:,i),SL1(:,:,i),ML1(:,:,i),m_paranoid,DBG)

        call MAGN(NSS(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),Xfield_2,ESO_TMP,zJ,THRS,dipso_tmp,s_so_tmp,nTempTotal,T,smagn, &
                  WL2(1:NSS(i),i),ZL2(:,i),SL2(:,:,i),ML2(:,:,i),m_paranoid,DBG)
        call mma_deallocate(ESO_TMP)
        call mma_deallocate(dipso_tmp)
        call mma_deallocate(s_so_tmp)
        call mma_allocate(ESO_TMP,nexch(i),label='ESO_TMP')
        call mma_allocate(dipso_tmp,3,nexch(i),nexch(i),label='dipso_tmp')
        call mma_allocate(s_so_tmp,3,nexch(i),nexch(i),label='s_so_tmp')
        ESO_TMP(:) = ESO(i,1:NEXCH(i))
        dipso_tmp(:,:,:) = DIPSO(i,:,1:NEXCH(i),1:NEXCH(i))
        s_so_tmp(:,:,:) = S_SO(i,:,1:NEXCH(i),1:NEXCH(i))
        ! only local "exchange states":
        call MAGN(NEXCH(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField,ESO_TMP,zJ,THRS,dipso_tmp,s_so_tmp,nTempTotal,T,smagn, &
                  WR0(1:NEXCH(i),i),ZR0(:,i),SR0(:,:,i),MR0(:,:,i),m_paranoid,DBG)

        call MAGN(NEXCH(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField_1,ESO_TMP,zJ,THRS,dipso_tmp,s_so_tmp,nTempTotal,T,smagn, &
                  WR1(1:NEXCH(i),i),ZR1(:,i),SR1(:,:,i),MR1(:,:,i),m_paranoid,DBG)

        call MAGN(NEXCH(i),NEXCH(i),dHX(iM),dHY(iM),dHZ(iM),XField_2,ESO_TMP,zJ,THRS,dipso_tmp,s_so_tmp,nTempTotal,T,smagn, &
                  WR2(1:NEXCH(i),i),ZR2(:,i),SR2(:,:,i),MR2(:,:,i),m_paranoid,DBG)
        call mma_deallocate(ESO_TMP)
        call mma_deallocate(dipso_tmp)
        call mma_deallocate(s_so_tmp)
#       ifdef _DEBUGPRINT_
        do iT=1,nTempTotal
          write(u6,'(A,i1,A,i3,A,3ES22.14)') 'ML(',i,',L,',iT,')=',(ML2(l,iT,i)-ML1(l,iT,i),l=1,3)
          write(u6,'(A,i1,A,i3,A,3ES22.14)') 'MR(',i,',L,',iT,')=',(MR2(l,iT,i)-MR1(l,iT,i),l=1,3)
        end do
#       endif
      else
        ZL0(:,i) = Zero
        ZR0(:,i) = Zero
        ZL1(:,i) = Zero
        ZR1(:,i) = Zero
        ZL2(:,i) = Zero
        ZR2(:,i) = Zero
        SL0(:,:,i) = Zero
        SR0(:,:,i) = Zero
        SL1(:,:,i) = Zero
        SR1(:,:,i) = Zero
        SL2(:,:,i) = Zero
        SR2(:,:,i) = Zero
        ML0(:,:,i) = Zero
        MR0(:,:,i) = Zero
        ML1(:,:,i) = Zero
        MR1(:,:,i) = Zero
        ML2(:,:,i) = Zero
        MR2(:,:,i) = Zero
      end if ! NSS(i) > NEXCH(i)
    end do ! i=1,nneq

    ibuf1 = nTempTotal*nneq
    ibuf3 = 3*ibuf1
    call Add_Info('dM/dH    ML0',[dnrm2_(ibuf3,ML0,1)],1,8)
    call Add_Info('dM/dH    SL0',[dnrm2_(ibuf3,SL0,1)],1,8)
    call Add_Info('dM/dH    ZL0',[dnrm2_(ibuf1,ZL0,1)],1,8)
    call Add_Info('dM/dH    WL0',[dnrm2_(nLoc*nneq,WL0,1)],1,8)

    call Add_Info('dM/dH    MR0',[dnrm2_(ibuf3,MR0,1)],1,8)
    call Add_Info('dM/dH    SR0',[dnrm2_(ibuf3,SR0,1)],1,8)
    call Add_Info('dM/dH    ZR0',[dnrm2_(ibuf1,ZR0,1)],1,8)
    call Add_Info('dM/dH    WR0',[dnrm2_(nLoc*nneq,WR0,1)],1,8)

    call Add_Info('dM/dH    ML1',[dnrm2_(ibuf3,ML1,1)],1,8)
    call Add_Info('dM/dH    SL1',[dnrm2_(ibuf3,SL1,1)],1,8)
    call Add_Info('dM/dH    ZL1',[dnrm2_(ibuf1,ZL1,1)],1,8)
    call Add_Info('dM/dH    WL1',[dnrm2_(nLoc*nneq,WL1,1)],1,8)

    call Add_Info('dM/dH    MR1',[dnrm2_(ibuf3,MR1,1)],1,8)
    call Add_Info('dM/dH    SR1',[dnrm2_(ibuf3,SR1,1)],1,8)
    call Add_Info('dM/dH    ZR1',[dnrm2_(ibuf1,ZR1,1)],1,8)
    call Add_Info('dM/dH    WR1',[dnrm2_(nLoc*nneq,WR1,1)],1,8)

    call Add_Info('dM/dH    ML2',[dnrm2_(ibuf3,ML2,1)],1,8)
    call Add_Info('dM/dH    SL2',[dnrm2_(ibuf3,SL2,1)],1,8)
    call Add_Info('dM/dH    ZL2',[dnrm2_(ibuf1,ZL2,1)],1,8)
    call Add_Info('dM/dH    WL2',[dnrm2_(nLoc*nneq,WL2,1)],1,8)

    call Add_Info('dM/dH    MR2',[dnrm2_(ibuf3,MR2,1)],1,8)
    call Add_Info('dM/dH    SR2',[dnrm2_(ibuf3,SR2,1)],1,8)
    call Add_Info('dM/dH    ZR2',[dnrm2_(ibuf1,ZR2,1)],1,8)
    call Add_Info('dM/dH    WR2',[dnrm2_(nLoc*nneq,WR2,1)],1,8)

    ! expand the basis and rotate local vectors to the general
    ! coordinate system:

    isite = 0
    do i=1,NNEQ
      do j=1,NEQ(i)
        isite = isite+1
        ! statistical distributions
        ZLT0(isite,:) = ZL0(:,i)
        ZRT0(isite,:) = ZR0(:,i)
        ZLT1(isite,:) = ZL1(:,i)
        ZRT1(isite,:) = ZR1(:,i)
        ZLT2(isite,:) = ZL2(:,i)
        ZRT2(isite,:) = ZR2(:,i)
        ! magnetizations:
        !    use R_rot matrices, which have determinant +1.
        !  note that  R_lg matrices have arbitrary determinant.
        do iT=1,nTempTotal
          do n=1,3
            MLT0(isite,:,iT) = MLT0(isite,:,iT)+rrot(i,j,:,n)*ML0(n,iT,i)
            SLT0(isite,:,iT) = SLT0(isite,:,iT)+rrot(i,j,:,n)*SL0(n,iT,i)
            MRT0(isite,:,iT) = MRT0(isite,:,iT)+rrot(i,j,:,n)*MR0(n,iT,i)
            SRT0(isite,:,iT) = SRT0(isite,:,iT)+rrot(i,j,:,n)*SR0(n,iT,i)

            MLT1(isite,:,iT) = MLT1(isite,:,iT)+rrot(i,j,:,n)*ML1(n,iT,i)
            SLT1(isite,:,iT) = SLT1(isite,:,iT)+rrot(i,j,:,n)*SL1(n,iT,i)
            MRT1(isite,:,iT) = MRT1(isite,:,iT)+rrot(i,j,:,n)*MR1(n,iT,i)
            SRT1(isite,:,iT) = SRT1(isite,:,iT)+rrot(i,j,:,n)*SR1(n,iT,i)

            MLT2(isite,:,iT) = MLT2(isite,:,iT)+rrot(i,j,:,n)*ML2(n,iT,i)
            SLT2(isite,:,iT) = SLT2(isite,:,iT)+rrot(i,j,:,n)*SL2(n,iT,i)
            MRT2(isite,:,iT) = MRT2(isite,:,iT)+rrot(i,j,:,n)*MR2(n,iT,i)
            SRT2(isite,:,iT) = SRT2(isite,:,iT)+rrot(i,j,:,n)*SR2(n,iT,i)
          end do
        end do

      end do ! j, neq(i)
    end do ! i, nneq
    ibuf = 0
    ibuf = 3*nTempTotal*nCenter
    call Add_Info('dM/dH    MLT0',[dnrm2_(ibuf,MLT0,1)],1,8)
    call Add_Info('dM/dH    SLT0',[dnrm2_(ibuf,SLT0,1)],1,8)
    call Add_Info('dM/dH    MRT0',[dnrm2_(ibuf,MRT0,1)],1,8)
    call Add_Info('dM/dH    SRT0',[dnrm2_(ibuf,SRT0,1)],1,8)

    call Add_Info('dM/dH    MLT1',[dnrm2_(ibuf,MLT1,1)],1,8)
    call Add_Info('dM/dH    SLT1',[dnrm2_(ibuf,SLT1,1)],1,8)
    call Add_Info('dM/dH    MRT1',[dnrm2_(ibuf,MRT1,1)],1,8)
    call Add_Info('dM/dH    SRT1',[dnrm2_(ibuf,SRT1,1)],1,8)

    call Add_Info('dM/dH    MLT2',[dnrm2_(ibuf,MLT2,1)],1,8)
    call Add_Info('dM/dH    SLT2',[dnrm2_(ibuf,SLT2,1)],1,8)
    call Add_Info('dM/dH    MRT2',[dnrm2_(ibuf,MRT2,1)],1,8)
    call Add_Info('dM/dH    SRT2',[dnrm2_(ibuf,SRT2,1)],1,8)
    ! compute the total magnetizations according to the derived formulas:
#   ifdef _DEBUGPRINT_
    write(u6,*) 'check point MLT0'
#   endif
    do iT=1,nTempTotal
      if (smagn) then
        call MSUM(nCenter,Sex0(:,iT),Zex0(iT),SLT0(:,:,iT),ZLT0(:,iT),SRT0(:,:,iT),ZRT0(:,iT),iopt,ST0(:,iT),ZT0(iT))
        call MSUM(nCenter,Sex1(:,iT),Zex1(iT),SLT1(:,:,iT),ZLT1(:,iT),SRT1(:,:,iT),ZRT1(:,iT),iopt,ST1(:,iT),ZT1(iT))
        call MSUM(nCenter,Sex2(:,iT),Zex2(iT),SLT2(:,:,iT),ZLT2(:,iT),SRT2(:,:,iT),ZRT2(:,iT),iopt,ST2(:,iT),ZT2(iT))
      end if

#     ifdef _DEBUGPRINT_
      write(u6,*) 'check point MSUM'
#     endif
      call MSUM(nCenter,Mex0(:,iT),Zex0(iT),MLT0(:,:,iT),ZLT0(:,iT),MRT0(:,:,iT),ZRT0(:,iT),iopt,MT0(:,iT),ZT0(iT))
      call MSUM(nCenter,Mex1(:,iT),Zex1(iT),MLT1(:,:,iT),ZLT1(:,iT),MRT1(:,:,iT),ZRT1(:,iT),iopt,MT1(:,iT),ZT1(iT))
      call MSUM(nCenter,Mex2(:,iT),Zex2(iT),MLT2(:,:,iT),ZLT2(:,iT),MRT2(:,:,iT),ZRT2(:,iT),iopt,MT2(:,iT),ZT2(iT))
    end do

  end if ! m_accurate

  ibuf = 3*nTempTotal
  call Add_Info('dM/dH    ST0',[dnrm2_(ibuf,ST0,1)],1,6)
  call Add_Info('dM/dH    ST1',[dnrm2_(ibuf,ST1,1)],1,6)
  call Add_Info('dM/dH    ST2',[dnrm2_(ibuf,ST2,1)],1,6)
  call Add_Info('dM/dH    MT0',[dnrm2_(ibuf,MT0,1)],1,6)
  call Add_Info('dM/dH    MT1',[dnrm2_(ibuf,MT1,1)],1,6)
  call Add_Info('dM/dH    MT2',[dnrm2_(ibuf,MT2,1)],1,6)
  ! computing the AVERAGE MOMENTS calculated at different temperatures (T(i))
  do iT=1,nTempTotal
    F1 = T(iT)*cm3tomB/dltXF
    F2 = T(iT)*cm3tomB/Xfield

    ! dM/dH model
    XTtens_dMdH(iM,:,iT) = (MT2(:,iT)-MT1(:,iT))*F1

    ! M/H model
    XTtens_MH(iM,:,iT) = MT0(:,iT)*F2
  end do
  ! ///  closing the loops over field directions
end do ! iM (nDirX)
! computing the XT as tensor's average:
XTM_dMdH(:) = (XTtens_dMdH(1,1,:)+XTtens_dMdH(2,2,:)+XTtens_dMdH(3,3,:))/Three

XTM_MH(:) = (XTtens_MH(1,1,:)+XTtens_MH(2,2,:)+XTtens_MH(3,3,:))/Three

ibuf = 9*nTempTotal
call Add_Info('dM/dH    XTtens_dMdH',[dnrm2_(ibuf,XTtens_dMdH,1)],1,5)
call Add_Info('dM/dH    XTtens_MH',[dnrm2_(ibuf,XTtens_MH,1)],1,6)
ibuf = nTempTotal
call Add_Info('dM/dH    XTM_dMdH',[dnrm2_(ibuf,XTM_dMdH,1)],1,5)
call Add_Info('dM/dH    XTM_MH',[dnrm2_(ibuf,XTM_MH,1)],1,6)
! ----------------------------------------------------------------------
! WRITING SOME OF THE OUTPUT....
! ----------------------------------------------------------------------
write(u6,*)
write(u6,'(A)') '--------------------------------------------------------------------------|'
write(u6,'(A)') '     |     T      | Statistical |   CHI*T     |   CHI*T     |   CHI*T     |'
write(u6,'(A,F6.3,A,F6.3,A)') '     |            |  Sum (Z)    | H = 0.000 T | H =',Xfield,' T | H =',XField,' T |'
write(u6,'(A)') '     |            |             |             | X = dM/dH   | X = M/H     |'
write(u6,'(A)') '-----|--------------------------------------------------------------------|'
write(u6,'(A)') 'Units|   kelvin   |    ---      |  cm3*K/mol  |  cm3*K/mol  |  cm3*K/mol  |'
write(u6,'(A)') '-----|--------------------------------------------------------------------|'

do iT=1,nT
  jT = iT+nTempMagn
  write(u6,'(A,F11.6,A,ES12.5,A,F12.8,A,F12.8,A,F12.7,A)') '     |',T(jT),' |',ZT0(jT),' |',XT_no_field(jT),' |',XTM_dMdH(jT), &
                                                           ' |',XTM_MH(jT),' |'
end do
write(u6,'(A)') '-----|--------------------------------------------------------------------|'

!  calcualtion of the standard deviation:
if (tinput) then
  write(u6,'(a,5x, f20.14)') 'ST.DEV: X= dM/dH:',dev(nT,XTM_dMdH(1+nTempMagn:),chit_exp(1+nTempMagn:))
  write(u6,'(a,5x, f20.14)') 'ST.DEV: X=  M/H :',dev(nT,XTM_MH(1+nTempMagn:),chit_exp(1+nTempMagn:))
  write(u6,'(A)') '-----|--------------------------------------------------------------------|'
end if !tinput

!-------------------------  PLOTs -------------------------------------!
write(label,'(A)') 'with_field_M_over_H'
call xFlush(u6)
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T(1+nTempMagn:),XTM_MH(1+nTempMagn:),chit_exp(1+nTempMagn:))
  else
    call plot_XT_no_Exp(label,nT,T(1+nTempMagn:),XTM_MH(1+nTempMagn:))
  end if
end if

write(label,'(A)') 'with_field_dM_over_dH'
call xFlush(u6)
if (DoPlot) then
  if (tinput) then
    call plot_XT_with_Exp(label,nT,T(1+nTempMagn:),XTM_dMdH(1+nTempMagn:),chit_exp(1+nTempMagn:))
  else
    call plot_XT_no_Exp(label,nT,T(1+nTempMagn:),XTM_dMdH(1+nTempMagn:))
  end if
end if
!------------------------- END PLOTs ----------------------------------!

! print out the main VAN VLECK SUSCEPTIBILITY TENSOR, its main values and main axes:
write(u6,'(/)')
write(u6,'(2A)') repeat('-',110),'|'
write(u6,'(25X,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR the  X=dM/dH  model,  in cm3*K/mol','|'
write(u6,'(2A)') repeat('-',110),'|'
write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
do iT=1,nT
  jT = iT+nTempMagn
  call DIAG_R2(XTtens_dMdH(:,:,jT),3,info,wt,zt)
  write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- z ----|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(XTtens_dMdH(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                  (zt(j,1),j=1,3),'|'
  write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(XTtens_dMdH(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                           (zt(j,2),j=1,3),'|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(XTtens_dMdH(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                  (zt(j,3),j=1,3),'|'
end do
write(u6,'(2A)') repeat('-',110),'|'

write(u6,'(/)')
write(u6,'(2A)') repeat('-',110),'|'

write(u6,'(25X,A,15x,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR FOR the  X = M/H  model,  in cm3*K/mol','|'
write(u6,'(2A)') repeat('-',110),'|'
write(u6,'(A)') '     T(K)   |   |          Susceptibility Tensor      |    Main Values  |               Main Axes             |'
do iT=1,nT
  jT = iT+nTempMagn
  call DIAG_R2(XTtens_MH(:,:,jT),3,info,wt,zt)
  write(u6,'(A)') '------------|---|------- x --------- y --------- z ---|-----------------|------ x --------- y --------- z ----|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | x |',(XTtens_MH(1,j,jT),j=1,3),' |  X:',wt(1),'|', &
                                                  (zt(j,1),j=1,3),'|'
  write(u6,'(F11.6,1x,A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') T(jT),'| y |',(XTtens_MH(2,j,jT),j=1,3),' |  Y:',wt(2),'|', &
                                                           (zt(j,2),j=1,3),'|'
  write(u6,'(A,3F12.6,A,F12.6,1x,A,3F12.8,1x,A)') '            | z |',(XTtens_MH(3,j,jT),j=1,3),' |  Z:',wt(3),'|', &
                                                  (zt(j,3),j=1,3),'|'
end do
write(u6,'(2A)') repeat('-',110),'|'

!=======================================================================
call mma_deallocate(WEX0)
call mma_deallocate(WEX1)
call mma_deallocate(WEX2)
call mma_deallocate(WL0)
call mma_deallocate(WL1)
call mma_deallocate(WL2)
call mma_deallocate(WR0)
call mma_deallocate(WR1)
call mma_deallocate(WR2)

call mma_deallocate(ZL0)
call mma_deallocate(ZR0)
call mma_deallocate(ZL1)
call mma_deallocate(ZR1)
call mma_deallocate(ZL2)
call mma_deallocate(ZR2)
call mma_deallocate(SL0)
call mma_deallocate(SR0)
call mma_deallocate(SL1)
call mma_deallocate(SR1)
call mma_deallocate(SL2)
call mma_deallocate(SR2)
call mma_deallocate(ML0)
call mma_deallocate(MR0)
call mma_deallocate(ML1)
call mma_deallocate(MR1)
call mma_deallocate(ML2)
call mma_deallocate(MR2)

call mma_deallocate(ZEX0)
call mma_deallocate(ZEX1)
call mma_deallocate(ZEX2)
call mma_deallocate(ZT0)
call mma_deallocate(ZT1)
call mma_deallocate(ZT2)
call mma_deallocate(SEX0)
call mma_deallocate(SEX1)
call mma_deallocate(SEX2)
call mma_deallocate(MEX0)
call mma_deallocate(MEX1)
call mma_deallocate(MEX2)
call mma_deallocate(ST0)
call mma_deallocate(ST1)
call mma_deallocate(ST2)
call mma_deallocate(MT0)
call mma_deallocate(MT1)
call mma_deallocate(MT2)
call mma_deallocate(XTM_MH)
call mma_deallocate(XTM_dMdH)
call mma_deallocate(XTtens_MH)
call mma_deallocate(XTtens_dMdH)

call mma_deallocate(ZRT0)
call mma_deallocate(ZLT0)
call mma_deallocate(ZRT1)
call mma_deallocate(ZLT1)
call mma_deallocate(ZRT2)
call mma_deallocate(ZLT2)
call mma_deallocate(MRT0)
call mma_deallocate(MLT0)
call mma_deallocate(SRT0)
call mma_deallocate(SLT0)
call mma_deallocate(MRT1)
call mma_deallocate(MLT1)
call mma_deallocate(SRT1)
call mma_deallocate(SLT1)
call mma_deallocate(MRT2)
call mma_deallocate(MLT2)
call mma_deallocate(SRT2)
call mma_deallocate(SLT2)

return

end subroutine XT_dMoverdH
