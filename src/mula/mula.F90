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
! Copyright (C) Niclas Forsberg                                        *
!               2008,2009, Giovanni Ghigo                              *
!***********************************************************************

subroutine Mula(ireturn)
!  Purpose:
!  Calculate vibronic intensities.
!
!  Written by:
!  Niclas Forsberg,
!  Dept. of Theoretical Chemistry, Chemical Centre, Lund University.
!
!  Modified by:
!  Giovanni Ghigo,
!  Dip. Chimica Generale e Chimica Organica,
!  Universita' di Torino, ITALY.
!  Inter-System Crossing rate constant.

use mula_global, only: AtCoord1, AtCoord2, energy1, energy2, Hess1, Hess2, Huge_Print, inpUnit, ipow, m_plot, maxMax_n, MaxNumAt, &
                       mdim1, mdim2, n_plot, ndata, ndim1, ndim2, ngdim, NormModes, nPolyTerm, nvar, OscStr, t_dipin1, t_dipin2, &
                       TranDip, TranDipGrad, var, WriteVibLevels, yin1, yin2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half, auTocm
use Definitions, only: wp, iwp, u5, u6, ItoB

implicit none
integer(kind=iwp), intent(out) :: ireturn
integer(kind=iwp) :: i, iCode, iMaxYes, iMem, iOrd, iOsc, iPrint, iPrint_Save, iRateMem, iUMem, iv, jOrd, jOsc, k1, k2, l_a, &
                     l_harm, l_Hess_1, l_Hess_2, l_IntensityMat_1, l_IntensityMat_2, l_l, l_m_plot, l_n_plot, l_NormModes, &
                     l_TermMat_1, l_TermMat_2, lBatch, leftBatch, lLeft, lMBatch, lNDEC, lNINC, lNMAT, lNMAT0, m_max, max_dip, &
                     max_mInc, max_mOrd, max_nInc, max_nOrd, max_term, minQ, mTabDim, n_max, nBatch, nBond, nDimTot, new_n_max, &
                     nm, nOsc, nOscOld, nTabDim, nterm, NumOfAt, nvTabDim, nYes
real(kind=wp) :: CPE, CPTF0, CPTF1, CPTF2, CPTF3, CPTF4, CPTF5, CPTF6, det, det0, det1, det2, dh, dMinWind, dRho, FC00, GE1, GE2, &
                 max_err, stand_dev, T0, TIOE, TIOTF0, TIOTF1, TIOTF2, TIOTF3, TIOTF4, TIOTF5, TIOTF6
logical(kind=iwp) :: Cartesian, exists, find_minimum, ForceField, lExpan, lInCore, lISC, lOldCode, MatEl, use_weight
character(len=80) :: Title
character(len=4) :: filnam
integer(kind=iwp), allocatable :: Bond(:), Graph1(:,:), Graph2(:,:,:), InterVec(:), level1(:), level2(:), lVec(:), mDec(:,:), &
                                  mInc(:,:), mMat(:,:), nDec(:,:), nInc(:,:), nIndex(:,:), nMat(:,:), nMat0(:), nMaxQ(:), &
                                  nnTabDim(:), VibWind2(:)
real(kind=wp), allocatable :: alpha1(:,:), alpha2(:,:), anharmfreq1(:), anharmfreq2(:), AtCoord(:,:), Base(:,:), Base2(:,:), &
                              BaseInv(:,:), Bmat(:,:), C(:,:), C1(:,:), C2(:,:), D3_1(:,:,:), D3_2(:,:,:), D4_1(:,:,:,:), &
                              D4_2(:,:,:,:), E1(:), E2(:), eigenVec1(:,:), eigenVec2(:,:), G0(:,:), G1(:,:), G2(:,:), &
                              Gbis0(:,:,:,:), Gbis1(:,:,:,:), Gbis2(:,:,:,:), Gprm0(:,:,:), Gprm1(:,:,:), Gprm2(:,:,:), grad1(:), &
                              grad2(:), H1(:,:), H2(:,:), harmfreq1(:), harmfreq2(:), IntensityMat(:,:), Lambda(:), &
                              Mass(:), OccNumMat1(:,:,:), OccNumMat2(:,:,:), PED(:,:,:), PotCoef(:), qMat(:,:), r0(:), r00(:), &
                              r01(:), r02(:), r1(:), r2(:), rtemp(:), S1(:,:), S2(:,:), Smat(:,:), T4(:), temp(:,:), temp1(:,:), &
                              temp2(:,:), TermMat(:,:), TranDipGradInt(:,:), U1(:,:), U2(:,:), W(:,:), W1(:,:), W2(:,:), &
                              x_anharm1(:,:), x_anharm2(:,:)
character(len=80), allocatable :: trfName1(:), trfName2(:)
character(len=4), allocatable :: AtomLbl(:)
integer(kind=iwp), parameter :: MB = 1048576
integer(kind=iwp), external :: iPrintlevel, isfreeunit
#include "warnings.h"

! Initialize.
iReturn = 20
close(u5)
inpUnit = isfreeunit(12)
call molcas_open(inpUnit,'stdin')
iCode = 00         ! Default: NewCode, InCore
lOldCode = .false. ! X0  CGGn.4
lInCore = .true.   ! 0X  CGGn.4
lNMAT0 = 92        !     CGGn.4
lNMAT = 93         !     CGGn.4
lNINC = 94         !     CGGn.4
lNDEC = 95         !     CGGn.4
iPrint = iPrintLevel(-1)
iPrint_Save = iPrint
!open(inpUnit,'stdin')

! Read input file and Write header to log file.
call mma_allocate(AtomLbl,MaxNumAt,label='AtomLbl')
call mma_allocate(Mass,MaxNumAt,label='Mass')
call mma_allocate(InterVec,MaxNumAt*15,label='InterVec')
call mma_allocate(Bond,MaxNumAt*2,label='Bond')
call mma_allocate(trfName1,MaxNumAt,label='trfName1')
call mma_allocate(trfName2,MaxNumAt,label='trfName2')
call ReadInp(Title,AtomLbl,Mass,InterVec,Bond,nBond,nOsc,NumOfAt,trfName1,trfName2,m_max,n_max,max_dip,max_term,MatEl,ForceField, &
             Cartesian,lExpan,lISC,iCode,dMinWind)
if (lISC) then
  if ((iCode == 01) .or. (iCode == 11)) lOldCode = .true.
  if ((iCode == 10) .or. (iCode == 11)) lInCore = .false.
  if (energy1 < energy2) then
    write(u6,*) ' ******************** ERROR *******************'
    write(u6,*) ' Energy of State #2 must be lower than State #1'
    write(u6,*) ' **********i***********************************'
    write(u6,*)
    call Quit_OnUserError()
  end if
  MatEl = .false.
  m_max = 0
  if (iPrint >= 1) then
    write(u6,*)
    write(u6,*) ' ---------------------------------------------'
    write(u6,*) ' InterSystem Crossing rate constant evaluation'
    write(u6,*) '  - simple harmonic approximation is used.'
    write(u6,*) '  - m_max set to 0.'
    write(u6,*) '  - no plot file.'
    if (lOldCode) then
      write(u6,*) '  - Full index matrices evaluation.'
    else
      write(u6,*) '  - Reduced matrices evaluation.'
    end if
    if (lInCore) then
      write(u6,*) '  - In-core (memory) algorithm.'
    else
      write(u6,*) '  - Out-of-core (disk) algorithm.'
    end if
    write(u6,*) ' ---------------------------------------------'
    write(u6,*)
  end if
  iPrint = 0
end if
if ((iPrint >= 1) .and. (Title(1:1) /= ' ')) call WriteHeader(Title)

close(inpUnit)
call mma_allocate(r01,nOsc,label='r01')
call mma_allocate(r02,nOsc,label='r02')
call mma_allocate(r0,nOsc,label='r0')
call mma_allocate(r1,nOsc,label='r1')
call mma_allocate(r2,nOsc,label='r2')

call mma_allocate(grad1,nOsc,label='grad1')
call mma_allocate(grad2,nOsc,label='grad2')

call mma_allocate(G1,nOsc,nOsc,label='G1')
call mma_allocate(G2,nOsc,nOsc,label='G2')

call mma_allocate(eigenVec1,nOsc,nOsc,label='eigenVec1')
call mma_allocate(eigenVec2,nOsc,nOsc,label='eigenVec2')

call mma_allocate(harmfreq1,nOsc,label='harmfreq1')
call mma_allocate(harmfreq2,nOsc,label='harmfreq2')

call mma_allocate(anharmfreq1,nOsc,label='anharmfreq1')
call mma_allocate(anharmfreq2,nOsc,label='anharmfreq2')

call mma_allocate(qMat,3*NumOfAt,nOsc,label='qMat')
call mma_allocate(PED,nOsc,nOsc,nOsc,label='PED')

call mma_allocate(x_anharm1,nOsc,nOsc,label='x_anharm1')
call mma_allocate(x_anharm2,nOsc,nOsc,label='x_anharm2')

if ((max_term >= 0) .and. (max_term <= 2)) then
  ngdim = 1
else if ((max_term >= 3) .and. (max_term <= 4)) then
  ngdim = nosc
else
  write(u6,*) 'MULA error: max_term=',max_term
  write(u6,*) 'Allowed: 1,2,3, or 4.'
  call Quit_OnUserError()
end if
call mma_allocate(D3_1,ngdim,ngdim,ngdim,label='D3_1')
call mma_allocate(D3_2,ngdim,ngdim,ngdim,label='D3_2')
call mma_allocate(D4_1,ngdim,ngdim,ngdim,ngdim,label='D4_1')
call mma_allocate(D4_2,ngdim,ngdim,ngdim,ngdim,label='D4_2')
call mma_allocate(Gprm1,ngdim,ngdim,ngdim,label='Gprm1')
call mma_allocate(Gprm2,ngdim,ngdim,ngdim,label='Gprm2')
call mma_allocate(Gbis1,ngdim,ngdim,ngdim,ngdim,label='Gbis1')
call mma_allocate(Gbis2,ngdim,ngdim,ngdim,ngdim,label='Gbis2')

!----------------------------------------------------------------------!

! Either fit polynomial to energies and find minimum or use
! forcefield and equilibrium geometry given in input.

!----------------------      First State ------------------------------!

call mma_allocate(AtCoord,3,NumOfAt,label='AtCoord')

l_a = NumOfAt
l_Hess_1 = size(Hess1,1)
l_Hess_2 = size(Hess1,2)

call mma_allocate(PotCoef,nPolyTerm,label='PotCoef')
if (ForceField) then
  AtCoord(:,:) = AtCoord1
else
  find_minimum = .true.
  use_weight = .false.
  call PotFit(nPolyTerm,nvar,ndata,ipow,var,yin1,PotCoef,r01,nOsc,energy1,grad1,Hess1,D3_1,D4_1,trfName1,stand_dev,max_err, &
              find_minimum,max_term,use_weight,l_Hess_1,l_Hess_2,ngdim)
  call Int_To_Cart1(InterVec,r01,AtCoord,l_a,nOsc)
end if
call mma_deallocate(AtCoord1)

! Determine vibrational modes and their frequencies.
!D write(u6,*) ' MULA calling VIBFREQ.'
call VibFreq(AtCoord,r01,InterVec,Mass,Hess1,G1,Gprm1,Gbis1,harmfreq1,eigenVec1,qMat,PED,D3_1,D4_1,x_anharm1,anharmfreq1,max_term, &
             nOsc,NumOfAt)

!call WriteLog(PotCoef,AtomLbl,AtCoord,
if (iPrint >= 1) call WriteLog(PotCoef,AtomLbl,AtCoord,Mass,InterVec,stand_dev,max_err,energy1,Hess1,G1,eigenVec1,harmfreq1,qMat, &
                               Bond,nBond,r01,D3_1,D4_1,PED,x_anharm1,anharmfreq1,max_term,1,ForceField,NumOfAt,nOsc)

!----------------------      Second State -----------------------------!

l_Hess_1 = size(Hess2,1)
l_Hess_2 = size(Hess2,2)

if (ForceField) then
  AtCoord(:,:) = AtCoord2
else
  find_minimum = .true.
  use_weight = .false.
  call PotFit(nPolyTerm,nvar,ndata,ipow,var,yin2,PotCoef,r02,nOsc,energy2,grad2,Hess2,D3_2,D4_2,trfName2,stand_dev,max_err, &
              find_minimum,max_term,use_weight,l_Hess_1,l_Hess_2,ngdim)
  call Int_To_Cart1(InterVec,r02,AtCoord,l_a,nOsc)
end if
call mma_deallocate(AtCoord2)

! Determine vibrational modes and their frequencies.
call VibFreq(AtCoord,r02,InterVec,Mass,Hess2,G2,Gprm2,Gbis2,harmfreq2,eigenVec2,qMat,PED,D3_2,D4_2,x_anharm2,anharmfreq2,max_term, &
             nOsc,NumOfAt)
if (iPrint >= 1) call WriteLog(PotCoef,AtomLbl,AtCoord,Mass,InterVec,stand_dev,max_err,energy2,Hess2,G2,eigenVec2,harmfreq2,qMat, &
                               Bond,nBond,r02,D3_2,D4_2,PED,x_anharm2,anharmfreq2,max_term,2,ForceField,NumOfAt,nOsc)

call mma_deallocate(Bond)
call mma_deallocate(PotCoef)
call mma_deallocate(qMat)
call mma_deallocate(PED)

call mma_deallocate(grad1)
call mma_deallocate(grad2)

call mma_deallocate(anharmfreq1)
call mma_deallocate(anharmfreq2)

call mma_deallocate(D3_1)
call mma_deallocate(D3_2)
call mma_deallocate(D4_1)
call mma_deallocate(D4_2)

!----------------------------------------------------------------------!

! Geometry of intermediate oscillator

!----------------------------------------------------------------------!

T0 = energy2-energy1
call mma_allocate(C,nOsc,nOsc,label='C')
call mma_allocate(C1,nOsc,nOsc,label='C1')
call mma_allocate(C2,nOsc,nOsc,label='C2')
call mma_allocate(W,nOsc,nOsc,label='W')
call mma_allocate(W1,nOsc,nOsc,label='W1')
call mma_allocate(W2,nOsc,nOsc,label='W2')

call mma_allocate(temp,nOsc,nOsc,label='temp')

! Calculate W matrices.
do jOsc=1,nOsc
  W1(:,jOsc) = eigenVec1(:,jOsc)/sqrt(harmfreq1(jOsc))
  W2(:,jOsc) = eigenVec2(:,jOsc)/sqrt(harmfreq2(jOsc))
end do

! Calculate C = W^(-1).
call unitmat(C1,nOsc)
C2(:,:) = C1
temp(:,:) = W1
call Dool_MULA(temp,nOsc,nOsc,C1,nOsc,nOsc,det)
det1 = abs(One/det)
temp(:,:) = W2
call Dool_MULA(temp,nOsc,nOsc,C2,nOsc,nOsc,det)
det2 = abs(One/det)

! Calculate the expansion point geometry and save the full geometries for later use
if (iPrint >= 1) call ExpPointHeader()
call mma_allocate(r00,nOsc,label='r00')
call mma_allocate(alpha1,nOsc,nOsc,label='alpha1')
call mma_allocate(alpha2,nOsc,nOsc,label='alpha2')
call Calc_r00(C1,C2,C,W,alpha1,alpha2,r00,r01,r02,det0,det1,det2,FC00,nOsc)
call mma_deallocate(alpha1)
call mma_deallocate(alpha2)
r0(:) = r00
r1(:) = r01
r2(:) = r02
l_a = NumOfAt
call Int_To_Cart1(InterVec,r00,AtCoord,l_a,nOsc)
if (iPrint >= 1) call WriteCartCoord(AtomLbl,AtCoord,Mass,NumOfAt)
call mma_deallocate(r00)

call mma_deallocate(AtomLbl)

! If we only wanted the expansion point geometry, it's time to quit now.
if (lExpan) then
  write(u6,*) 'Bye bye'
  write(u6,*)
  call Quit(_RC_ALL_IS_WELL_)
end if

if (iPrint >= 1) call IntCalcHeader()

! Calculate term values.
GE1 = Zero
GE2 = GE1+T0
l_NormModes = size(NormModes)
do iOsc=1,nOsc
  jOsc = 1
  exists = .false.
  do while ((.not. exists) .and. (jOsc <= l_NormModes))
    exists = (iOsc == NormModes(jOsc))
    jOsc = jOsc+1
  end do
  if (.not. exists) then
    GE1 = GE1+Half*harmfreq1(iOsc)
    GE2 = GE2+Half*harmfreq2(iOsc)
  end if
end do
T0 = GE2-GE1

! Pick out the chosen modes.
nOscOld = nOsc
nOsc = l_NormModes
call mma_allocate(Lambda,nOsc,label='Lambda')
call mma_allocate(Base,nOscOld,nOsc,label='Base')
call mma_allocate(BaseInv,nOsc,nOscOld,label='BaseInv')
do jOsc=1,nOsc
  nm = NormModes(jOsc)
  do iOsc=1,nOscOld
    Base(iOsc,jOsc) = W2(iOsc,nm)
    BaseInv(jOsc,iOsc) = C2(nm,iOsc)
  end do
end do

call mma_allocate(temp1,nOscOld,nOsc,label='temp1')
call mma_allocate(temp2,nOscOld,nOscOld,label='temp2')

! First state.
! Subroutine SolveRedSec(Hess,Gmat,freq,C,W,det)
call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,Hess1,nOscOld,Base,nOscOld,Zero,temp1,nOscOld)
call mma_deallocate(Hess1)

call mma_allocate(Hess1,nOsc,nOsc,label='Hess1')
call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Base,nOscOld,temp1,nOscOld,Zero,Hess1,nOsc)
call unitmat(temp2,nOscOld)
call Dool_MULA(G1,nOsc,nOsc,temp2,nOscOld,nOscOld,det)
call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,temp2,nOscOld,Base,nOscOld,Zero,temp1,nOscOld)
call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Base,nOscOld,temp1,nOscOld,Zero,temp,nOsc)
call unitmat(G1,nOsc)
call Dool_MULA(temp,nOscOld,nOscOld,G1,nOsc,nOsc,det)
call SolveSecEq(Hess1,nOsc,temp,G1,Lambda)
do iv=1,nOsc
  harmfreq1(iv) = sqrt(abs(Lambda(iv)))
end do
if (iPrint >= 1) call WriteFreq(harmfreq1,NormModes,l_NormModes,'Frequencies of reduced problem, state 1')
do jOsc=1,nOsc
  W1(:,jOsc) = temp(:,jOsc)/sqrt(harmfreq1(jOsc))
end do
call unitmat(C1,nOsc)
temp(:,:) = W1
call Dool_MULA(temp,nOscOld,nOsc,C1,nOsc,nOsc,det)
det1 = abs(One/det)

! Second state.
call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,Hess2,nOscOld,Base,nOscOld,Zero,temp1,nOscOld)
call mma_deallocate(Hess2)

call mma_allocate(Hess2,nOsc,nOsc,label='Hess2')
call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Base,nOscOld,temp1,nOscOld,Zero,Hess2,nOsc)
call unitmat(temp2,nOscOld)
call Dool_MULA(G2,nOsc,nOsc,temp2,nOscOld,nOscOld,det)
call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,temp2,nOscOld,Base,nOscOld,Zero,temp1,nOscOld)
call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Base,nOscOld,temp1,nOscOld,Zero,temp,nOsc)
call unitmat(G2,nOsc)
call Dool_MULA(temp,nOsc,nOsc,G2,nOsc,nOsc,det)
call SolveSecEq(Hess2,nOsc,temp,G2,Lambda)
do iv=1,nOsc
  harmfreq2(iv) = sqrt(abs(Lambda(iv)))
end do
if (iPrint >= 1) call WriteFreq(harmfreq2,NormModes,l_NormModes,'Frequencies of reduced problem, state 2')
do jOsc=1,nOsc
  W2(:,jOsc) = temp(:,jOsc)/sqrt(harmfreq2(jOsc))
end do
call unitmat(C2,nOsc)
temp(:,:) = W2
call Dool_MULA(temp,nOsc,nOsc,C2,nOsc,nOsc,det)
det2 = abs(One/det)

call mma_deallocate(temp)
call mma_deallocate(temp1)
call mma_deallocate(temp2)

call mma_allocate(rtemp,nOscOld,label='rtemp')

rtemp(:) = r01
call DGEMM_('N','N',nOsc,1,nOscOld,One,BaseInv,nOsc,rtemp,nOscOld,Zero,r01,nOsc)

rtemp(:) = r02
call DGEMM_('N','N',nOsc,1,nOscOld,One,BaseInv,nOsc,rtemp,nOscOld,Zero,r02,nOsc)

call mma_deallocate(rtemp)

call mma_allocate(r00,nOsc,label='r00')
call mma_allocate(alpha1,nOsc,nOsc,label='alpha1')
call mma_allocate(alpha2,nOsc,nOsc,label='alpha2')
call Calc_r00(C1,C2,C,W,alpha1,alpha2,r00,r01,r02,det0,det1,det2,FC00,nOsc)
call mma_deallocate(alpha1)
call mma_deallocate(alpha2)

!Base2(:,:) = W2
!BaseInv2(:,:) = C2

!----------------------------------------------------------------------!

!Transform transition dipole gradients from cartesian to internal coordinates

!----------------------------------------------------------------------!

if (Forcefield .and. (max_dip > 0)) then
  call mma_allocate(Bmat,3*NumOfAt,nOscOld,label='Bmat')
  Bmat(:,:) = Zero
  call CalcS(AtCoord,InterVec,Bmat,nOscOld,NumOfAt)

  call mma_allocate(temp,nOscOld,nOscOld,label='temp')
  call mma_allocate(temp2,nOscOld,nOscOld,label='temp2')

  call DGEMM_('T','N',nOscOld,nOscOld,3*NumOfAt,One,Bmat,3*NumOfAt,Bmat,3*NumOfAt,Zero,temp,nOscOld)
  call unitmat(temp2,nOscOld)
  call Dool_MULA(temp,nOscOld,nOscOld,temp2,nOscOld,nOscOld,det)
  call mma_deallocate(temp)
  call mma_allocate(temp1,3*NumOfAt,nOscOld,label='temp1')

  call DGEMM_('N','N',3*NumOfAt,nOscOld,nOscOld,One,Bmat,3*NumOfAt,temp2,nOscOld,Zero,temp1,3*NumOfAt)
  call mma_deallocate(Bmat)
  call mma_deallocate(temp2)
  call mma_allocate(temp2,3*NumOfAt,nOsc,label='temp2')
  call DGEMM_('N','N',3*NumOfAt,nOsc,nOscOld,One,temp1,3*NumOfAt,Base,nOscOld,Zero,temp2,3*NumOfAt)
  call mma_deallocate(temp1)
  call mma_allocate(temp1,3*NumOfAt,nOsc,label='temp1')
  call DGEMM_('N','N',3*NumOfAt,nOsc,nOsc,One,temp2,3*NumOfAt,W,nOsc,Zero,temp1,3*NumOfAt)
  call mma_deallocate(temp2)
  call mma_allocate(TranDipGradInt,3,nOsc,label='TranDipGradInt')
  call DGEMM_('T','N',nOsc,1,3*NumOfAt,One,temp1,3*NumOfAt,TranDipGrad(1,1),3*NumOfAt,Zero,TranDipGradInt(1,1),nOsc)
  call DGEMM_('T','N',nOsc,1,3*NumOfAt,One,temp1,3*NumOfAt,TranDipGrad(2,1),3*NumOfAt,Zero,TranDipGradInt(2,1),nOsc)
  call DGEMM_('T','N',nOsc,1,3*NumOfAt,One,temp1,3*NumOfAt,TranDipGrad(3,1),3*NumOfAt,Zero,TranDipGradInt(3,1),nOsc)
  call mma_deallocate(temp1)

  call WriteDip(TranDipGradInt,NormModes,'Transition dipole gradient',nOsc)
end if

call mma_deallocate(eigenVec1)
call mma_deallocate(eigenVec2)
call mma_deallocate(Lambda)

!----------------------------------------------------------------------!

!  InterSystem Crossing section.
!  Author: Giovanni Ghigo
!          Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
!          Vers. 1  28 Dec-08 - 13 Jan-09

!----------------------------------------------------------------------!

if (iPrint >= 4) write(u6,*) ' FC(0-0)=',FC00
if (lISC) then

  iPrint = iPrint_Save
  if (WriteVibLevels) iPrint = 2
  if (Huge_Print) iPrint = 3
  if (iPrint >= 1) call ISCHeader()

  dRho = Zero
  call mma_allocate(nMaxQ,nOsc,label='nMaxQ')
  call ISC_Rho(iPrint,nOsc,new_n_max,dRho,energy1,energy2,minQ,dMinWind,nMaxQ,harmfreq1,harmfreq2)
  if (new_n_max < n_max) then
    n_max = new_n_max
    write(u6,*) ' n_max reduced to',n_max
  end if
  if (iPrint >= 3) write(u6,*) ' Actual n_max               =',n_max

  ! Set up excitation matrices

  ! Initial State m
  call TabDim(m_max,nOsc,nvTabDim)
  mTabDim = nvTabDim-1
  max_mOrd = mTabDim
  call TabDim(m_max-1,nOsc,nvTabDim)
  max_mInc = nvTabDim-1
  call mma_allocate(mMat,[0,mTabDim],[1,nOsc],label='mMat')
  call mma_allocate(mInc,[0,mTabDim],[1,nOsc],label='mInc')
  call mma_allocate(mDec,[0,mTabDim],[1,nOsc],label='mDec')
  mdim1 = mTabDim
  mdim2 = nOsc
  call MakeTab2(m_max,max_mOrd,max_mInc,mTabDim,mMat,mInc,mDec,nOsc)

  ! Initialize Final State n
  call TabDim(n_max,nOsc,nvTabDim)
  nTabDim = nvTabDim-1
  max_nOrd = nTabDim
  call TabDim(n_max-1,nOsc,nvTabDim)
  max_nInc = nvTabDim-1

  ! Memory estimation and algorithm selection
  call mma_maxINT(lLeft)
  if (lInCore) then
    if (lOldCode) then
      iMem = 7*(nTabDim+1)*nOsc/2
      if (iMem >= lLeft) then
        write(u6,*) ' Too much memory required (',iMem,' words).'
        write(u6,*) ' Switch to reduced matrices evaluation.'
        lOldCode = .false.
      end if
    else
      iMem = 3*(nTabDim+1)*nOsc/2
      if (iMem >= lLeft) then
        write(u6,*) ' Too much memory required (',iMem,' words).'
        write(u6,*) ' Out-of-Core (disk) algorithm will be used.'
        lInCore = .false.
      end if
    end if
  end if

  call Timing(CPTF0,CPE,TIOTF0,TIOE)
  if (lInCore) then

    if (lOldCode) then
      if (iPrint >= 2) then
        write(u6,*)
        write(u6,*) ' Index matrix evaluation.'
        if (iPrint >= 3) write(u6,*) ' Memory allocated for all index matrix:',iMem,' words,  ',iMem*ItoB/MB,' MB.'
        call XFlush(u6)
      end if
      call mma_allocate(nMat,[0,nTabDim],[1,nOsc],label='nMat')
      call mma_allocate(nInc,[0,nTabDim],[1,nOsc],label='nInc')
      call mma_allocate(nDec,[0,nTabDim],[1,nOsc],label='nDec')
      ndim1 = nTabDim
      ndim2 = nOsc
      call MakeTab2(n_max,max_nOrd,max_nInc,nTabDim,nMat,nInc,nDec,nOsc)
    else ! NewCode: nInc & nDec reduced
      iMem = (nTabDim+1)*nOsc
      if (iPrint >= 2) then
        write(u6,*)
        write(u6,*) ' First index matrix evaluation.'
        if (iPrint >= 3) write(u6,*) ' Memory allocated for first index matrix:',iMem,' words,  ',iMem*ItoB/MB,' MB.'
        call XFlush(u6)
      end if
      call mma_allocate(nMat,[0,nTabDim],[1,nOsc],label='nMat')
      call mma_allocate(Graph1,[0,n_max],[1,nOsc+1],label='Graph1')
      call mma_allocate(Graph2,[0,n_max],[0,n_max],[1,nOsc],label='Graph2')
      ndim1 = nTabDim
      ndim2 = nOsc
      call ISC_MakeTab2(n_max,max_nOrd,nTabDim,nMat,Graph1,Graph2,nOsc)
      call mma_deallocate(Graph1)
    end if
    call Timing(CPTF1,CPE,TIOTF1,TIOE)

    ! Definition of Vector with levels in the window

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Energy level screening.'
      call XFlush(u6)
    end if
    call mma_allocate(lVec,[0,nTabDim],label='lVec')
    call LogEVec(iPrint,nOsc,max_nOrd,minQ,nMaxQ,nMat,lVec,nYes)
    if (nYes <= 0) then
      write(u6,*)
      write(u6,*) ' ************ ERROR *************'
      write(u6,*) ' No energy levels in the Window !'
      write(u6,*) ' ********************************'
      write(u6,*)
      call Quit_OnUserError()
    end if
    call Timing(CPTF2,CPE,TIOTF2,TIOE)

    ! Energy levels calculation

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Energy level calculation.'
      if (iPrint >= 3) then
        write(u6,*) ' Memory allocated for energy matrix:',(nTabDim+1),' words,  ',(nTabDim+1)*ItoB/MB,' MB.'
      end if
      call XFlush(u6)
    end if
    call ISC_Ene(iPrint,nOsc,max_nOrd,nYes,nMat,nTabDim,GE1,GE2,harmfreq1,harmfreq2,x_anharm1,x_anharm2,dMinWind,dRho,lVec)
    if (nYes <= 0) then
      write(u6,*)
      write(u6,*) ' ************ ERROR *************'
      write(u6,*) ' No energy levels in the Window !'
      write(u6,*) ' ********************************'
      write(u6,*)
      if (iPrint < 4) call Quit_OnUserError()
    end if

    ! The State Window

    call mma_allocate(VibWind2,nYes,label='VibWind2')
    call MkVibWind2(nYes,iMaxYes,max_nOrd,lVec,VibWind2)
    call mma_deallocate(lVec)
    call Timing(CPTF3,CPE,TIOTF3,TIOE)
    call Timing(CPTF4,CPE,TIOTF4,TIOE)

    ! The remaining excitation matrices nInc & nDec

    if (lOldCode) then
      iMaxYes = nTabDim
    else
      iMem = 2*(iMaxYes+1)*nOsc
      if (iPrint >= 2) then
        write(u6,*)
        write(u6,*) ' nInc and nDec matrices generation.'
        if (iPrint >= 3) write(u6,*) ' Memory allocated for remaining index matrix:',iMem,' words,  ',iMem*ItoB/MB,' MB.'
        call XFlush(u6)
      end if
      call mma_allocate(nInc,[0,iMaxYes],[1,nOsc],label='nInc')
      call mma_allocate(nDec,[0,iMaxYes],[1,nOsc],label='nDec')
      call Mk_nIncDec(n_max,nTabDim,iMaxYes,nInc,nDec,nMat,Graph2,nOsc)
      call mma_deallocate(Graph2)
    end if
    call Timing(CPTF5,CPE,TIOTF5,TIOE)

    ! The ISC rate

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Franck-Condon factors evaluation.'
      call XFlush(u6)
    end if
    call ISC_Rate(iPrint,nOsc,max_nOrd,iMaxYes,nYes,dMinWind,VibWind2,C1,C2,W2,det0,det1,det2,C,W,r01,r02,nTabDim,nMat,nInc,nDec, &
                  m_max,n_max,FC00,dRho)
    call Timing(CPTF6,CPE,TIOTF6,TIOE)

    call mma_deallocate(VibWind2)

    call mma_deallocate(nMat)
    call mma_deallocate(nInc)
    call mma_deallocate(nDec)

  else ! Out-of-core.

    !GGt --- Out-of-core ---------------------------------------------------

    ! Open files

    filnam = 'MAT0'
    call DaName_mf_wa(lNMAT0,filnam)

    call mma_allocate(Graph1,[0,n_max],[1,nOsc+1],label='Graph1')
    call mma_allocate(Graph2,[0,n_max],[0,n_max],[1,nOsc],label='Graph2')
    call ISCD_MakeGraphs(n_max,max_nOrd,Graph1,Graph2,nOsc)
    call mma_deallocate(Graph1)
    call mma_allocate(nnTabDim,[0,nTabDim],label='nnTabDim')
    call mma_allocate(nMat0,nOsc,label='nMat0')
    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Index matrix evaluation.'
      if (iPrint >= 3) then
        write(u6,*) ' Memory allocated for address array:',(nTabDim+1),' words,  ',(nTabDim+1)*ItoB/MB,' MB.'
      end if
      call XFlush(u6)
    end if
    call ISCD_MakenMat(n_max,nOsc,lNMAT0,nTabDim,Graph2,nnTabDim,nMat0)
    call Timing(CPTF1,CPE,TIOTF1,TIOE)

    ! Definition of Vector with levels in the window

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Energy level screening.'
      call XFlush(u6)
    end if
    call mma_allocate(lVec,[0,nTabDim],label='lVec')
    call ISCD_LogEVec(iPrint,nOsc,max_nOrd,minQ,nYes,lNMAT0,nTabDim,nnTabDim,nMaxQ,nMat0,lVec)
    if (nYes <= 0) then
      write(u6,*)
      write(u6,*) ' ************ ERROR *************'
      write(u6,*) ' No energy levels in the Window !'
      write(u6,*) ' ********************************'
      write(u6,*)
      call Quit_OnUserError()
    end if
    call Timing(CPTF2,CPE,TIOTF2,TIOE)

    ! Energy levels calculation

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Energy level calculation.'
      if (iPrint >= 3) then
        write(u6,*) ' Memory allocated for energy matrix:',(nTabDim+1),' words,  ',(nTabDim+1)*ItoB/MB,' MB.'
      end if
      call XFlush(u6)
    end if
    call ISCD_Ene(iPrint,nOsc,max_nOrd,nYes,lNMAT0,nTabDim,GE1,GE2,harmfreq1,harmfreq2,x_anharm1,x_anharm2,dMinWind,dRho,nMat0, &
                  nnTabDim,lVec)
    if (nYes <= 0) then
      write(u6,*)
      write(u6,*) ' ************ ERROR *************'
      write(u6,*) ' No energy levels in the Window !'
      write(u6,*) ' ********************************'
      write(u6,*)
      if (iPrint < 4) call Quit_OnUserError()
    end if

    ! The State Window

    call mma_allocate(VibWind2,nYes,label='VibWind2')
    call MkVibWind2(nYes,iMaxYes,max_nOrd,lVec,VibWind2)
    call mma_deallocate(lVec)
    call Timing(CPTF3,CPE,TIOTF3,TIOE)

    ! Memory estimation & nMat transfer

    if (lOldCode) then
      iMaxYes = nTabDim
    end if
    call mma_maxINT(lLeft)
    if (iPrint >= 3) then
      write(u6,'(A,I10,A,I4,A)') '  Available memory                          :',lLeft,' words,  ',lLeft*iToB/MB,' MB.'
    end if
    iUMem = nTabDim+1  ! Memory for U, 1=> No Hot states
    iRateMem = iUMem+nTabDim+2+13*nOsc*nOsc+5*nOsc
    iRateMem = int(2.1_wp*iRateMem) ! to convert from INTE to REAL
    lLeft = lLeft-iRateMem
    if (iRateMem < 0) then ! .or. (iRateMem > 2048*MB/ItoB)) then
      write(u6,'(A,I10,A,I4,A)') '  Estimated memory for FC factors evaluation:',iRateMem,' words,  ',iRateMem*ItoB/MB,' MB.'
      write(u6,'(A,I10,A,I4,A)') '  Memory left for batches                   :',lLeft,' words,  ',lLeft*ItoB/MB,' MB.'
      write(u6,*)
      write(u6,*) ' ****************** ERROR ********************'
      write(u6,*) ' Not enough memory for FC factors evaluation !'
      write(u6,*) ' *********************************************'
      write(u6,*)
      call Quit_OnUserError()
    end if
    lMBatch = int(lLeft/3)
    lMBatch = min(lMBatch,512*MB/ItoB) ! Max 512 MB
    lMBatch = max(lMBatch,64*MB/ItoB)  ! Min  64 MB
    !GGt lMBatch = MB  ! Test only !!!!
    lBatch = int(lMBatch/nOsc)
    if (lBatch > (iMaxYes+1)) lBatch = (iMaxYes+1)
    lMBatch = lBatch*nOsc
    nBatch = int((iMaxYes+1)/lBatch)
    leftBatch = (iMaxYes+1)-nBatch*lBatch
    if (nBatch+1 > maxMax_n) then
      write(u6,*)
      write(u6,*) ' ******************** ERROR *******************'
      write(u6,*) ' Not enough number of batches for Out-of-core !'
      write(u6,*) ' Increase maxMax_n in src/mula/io_mula.fh'
      write(u6,*) ' **********************************************'
      write(u6,*)
      call Quit_OnUserError()
    end if
    if (iPrint >= 2) then
      write(u6,*)
      if (iPrint >= 3) then
        write(u6,'(A,I10,A,I4,A)') '  Estimated memory for FC factors evaluation:',iRateMem,' words,  ',iRateMem*ItoB/MB,' MB.'
        write(u6,'(A,I10,A,I4,A)') '  Memory left for batches                   :',lLeft,' words,  ',lLeft*ItoB/MB,' MB.'
        write(u6,'(A,I10,A,I4,A)') '  Memory allocated for batch                :',lMBatch,' words,  ',lMBatch*ItoB/MB,' MB.'
        write(u6,'(A,I8)') '  * Elements for batch:',lBatch
        write(u6,'(A,I8,A,I8,A)') '  * Number of batches :',nBatch,' (',lBatch*nBatch,' elements)'
        write(u6,'(A,I8)') '  * Residual elements :',leftBatch
      end if
      write(u6,*) ' Index matrix reloading.'
      call XFlush(u6)
    end if

    call mma_allocate(nIndex,[1,3],[0,maxMax_n],label='nIndex')
    nIndex(:,:) = 0

    filnam = 'NMAT'
    call DaName_mf_wa(lNMAT,filnam)
    call mma_allocate(nMat,nOsc,lBatch,label='nMat')
    call ISCD_ReloadNMAT(nTabDim,iMaxYes,nOsc,lNMAT0,lNMAT,lBatch,nBatch,leftBatch,nIndex,nnTabDim,nMAT0,nMAT)
    call mma_deallocate(nMat0)
    call Timing(CPTF4,CPE,TIOTF4,TIOE)

    ! The remaining excitation matrices nInc & nDec

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' nInc and nDec matrices generation.'
      call XFlush(u6)
    end if
    filnam = 'NINC'
    call DaName_mf_wa(lNINC,filnam)
    filnam = 'NDEC'
    call DaName_mf_wa(lNDEC,filnam)
    call mma_allocate(nInc,nOsc,lBatch,label='nInc')
    call mma_allocate(nDec,nOsc,lBatch,label='nDec')
    call ISCD_MakenIncDec(n_max,iMaxYes,nOsc,lNMAT,lNINC,lNDEC,lBatch,nBatch,leftBatch,nIndex,Graph2,nMAT,nInc,nDec)
    call mma_deallocate(Graph2)
    call Timing(CPTF5,CPE,TIOTF5,TIOE)

    ! The ISC rate

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Franck-Condon factors evaluation.'
      call XFlush(u6)
    end if
    call ISCD_Rate(iPrint,nOsc,max_nOrd,iMaxYes,nYes,dMinWind,lBatch,nIndex,VibWind2,lNMAT0,lNMAT,lNINC,lNDEC,nTabDim,nnTabDim,C1, &
                   C2,W2,det0,det1,det2,C,W,r01,r02,m_max,n_max,FC00,dRho,nMat,nInc,nDec)
    call Timing(CPTF6,CPE,TIOTF6,TIOE)

    call mma_deallocate(nIndex)
    call mma_deallocate(VibWind2)

    call mma_deallocate(nMat)
    call mma_deallocate(nInc)
    call mma_deallocate(nDec)
    call mma_deallocate(nnTabDim)

    call DaClos(lNDEC)
    call DaClos(lNINC)
    call DaClos(lNMAT)
    call DaClos(lNMAT0)

  end if

  if (iPrint >= 3) then
    write(u6,*)
    write(u6,'(A)') ' Timing informations (sec.):'
    write(u6,'(A)') ' ========================================='
    write(u6,'(A)') ' MODULE:                       CPU Elapsed'
    write(u6,'(A,2F8.2)') ' Matrix evaluation        ',CPTF1-CPTF0,TIOTF1-TIOTF0
    write(u6,'(A,2F8.2)') ' Level screening          ',CPTF2-CPTF1,TIOTF2-TIOTF1
    write(u6,'(A,2F8.2)') ' Energy calculation       ',CPTF3-CPTF2,TIOTF3-TIOTF2
    write(u6,'(A,2F8.2)') ' Matrix reloading         ',CPTF4-CPTF3,TIOTF4-TIOTF3
    write(u6,'(A,2F8.2)') ' nInc and nDec generation ',CPTF5-CPTF4,TIOTF5-TIOTF4
    write(u6,'(A,2F8.2)') ' Franck-Condon evaluation ',CPTF6-CPTF5,TIOTF6-TIOTF5
    write(u6,'(A,2F8.2)') ' TOTAL                    ',CPTF6-CPTF0,TIOTF6-TIOTF0
    write(u6,'(A)') ' -----------------------------------------'
    write(u6,*)
  end if

  call mma_deallocate(mMat)
  call mma_deallocate(mInc)
  call mma_deallocate(mDec)

  call mma_deallocate(nMaxQ)

  call mma_deallocate(r00)
  call mma_deallocate(r01)
  call mma_deallocate(r02)

  call mma_deallocate(Hess1)
  call mma_deallocate(Hess2)

  call mma_deallocate(Base)
  call mma_deallocate(BaseInv)

  call mma_deallocate(C)
  call mma_deallocate(C1)
  call mma_deallocate(C2)
  call mma_deallocate(W)
  call mma_deallocate(W1)
  call mma_deallocate(W2)

  call mma_deallocate(AtCoord)

  call mma_deallocate(Gprm1)
  call mma_deallocate(Gprm2)
  call mma_deallocate(Gbis1)
  call mma_deallocate(Gbis2)

  call mma_deallocate(x_anharm1)
  call mma_deallocate(x_anharm2)
  call mma_deallocate(harmfreq1)
  call mma_deallocate(harmfreq2)

  call mma_deallocate(G1)
  call mma_deallocate(G2)

  call mma_deallocate(r0)
  call mma_deallocate(r1)
  call mma_deallocate(r2)

  call mma_deallocate(m_plot)
  call mma_deallocate(n_plot)

  call mma_deallocate(NormModes)
  call mma_deallocate(TranDipGrad)

  lISC = .false.

else
  !GGn -----------------------------------------------------------------

  !--------------------------------------------------------------------!

  ! Set up excitation matrices

  !--------------------------------------------------------------------!

  call mma_allocate(G0,nOsc,nOsc,label='G0')
  call mma_allocate(Gprm0,ngdim,ngdim,ngdim,label='Gprm0')
  call mma_allocate(Gbis0,ngdim,ngdim,ngdim,ngdim,label='Gbis0')

  ! Calculate dimensions given max level of excitation for the different states.
  call TabDim(m_max,nOsc,nvTabDim)
  mTabDim = nvTabDim-1
  call TabDim(n_max,nOsc,nvTabDim)
  nTabDim = nvTabDim-1

  ! Set up mMat for L.
  max_mOrd = mTabDim
  call TabDim(m_max-1,nOsc,nvTabDim)
  max_mInc = nvTabDim-1

  call mma_allocate(mMat,[0,mTabDim],[1,nOsc],label='mMat')
  call mma_allocate(mInc,[0,mTabDim],[1,nOsc],label='mInc')
  call mma_allocate(mDec,[0,mTabDim],[1,nOsc],label='mDec')
  ! Put dimensions into common block:
  mdim1 = mTabDim
  mdim2 = nOsc
  call MakeTab2(m_max,max_mOrd,max_mInc,mTabDim,mMat,mInc,mDec,nOsc)

  ! Set up nMat for U.
  max_nOrd = nTabDim
  call TabDim(n_max-1,nOsc,nvTabDim)
  max_nInc = nvTabDim-1

  call mma_allocate(nMat,[0,nTabDim],[1,nOsc],label='nMat')
  call mma_allocate(nInc,[0,nTabDim],[1,nOsc],label='nInc')
  call mma_allocate(nDec,[0,nTabDim],[1,nOsc],label='nDec')
  ! Put dimensions into common block:
  ndim1 = nTabDim
  ndim2 = nOsc
  call MakeTab2(n_max,max_nOrd,max_nInc,nTabDim,nMat,nInc,nDec,nOsc)

  if (.not. MatEl) then
    l_TermMat_1 = max_mOrd
    l_TermMat_2 = max_nOrd
    call mma_allocate(TermMat,[0,l_TermMat_1],[0,l_TermMat_2],label='TermMat')
    call mma_allocate(level1,nOsc,label='level1')
    call mma_allocate(level2,nOsc,label='level2')
    do jOrd=0,max_nOrd
      level2(:) = nMat(jOrd,:)
      do iOrd=0,max_mOrd
        level1(:) = mMat(iOrd,:)
        l_harm = nOsc
        call TransEnergy(GE1,x_anharm1,harmfreq1,level1,GE2,x_anharm2,harmfreq2,level2,TermMat(iOrd,jOrd),l_harm)
      end do
    end do
    call mma_deallocate(level1)
    call mma_deallocate(level2)

  end if

  !--------------------------------------------------------------------!

  ! Calculate vibrational Hessian if a variational calculation is chosen

  !--------------------------------------------------------------------!

  l_l = 1
  call mma_allocate(U1,l_l,l_l,label='U1')
  call mma_allocate(E1,l_l,label='E1')
  call mma_allocate(U2,l_l,l_l,label='U2')
  call mma_allocate(E2,l_l,label='E2')
  call mma_allocate(OccNumMat1,l_l,l_l,l_l,label='OccNumMat1')
  call mma_allocate(OccNumMat2,l_l,l_l,l_l,label='OccNumMat2')
  if (MatEl) then  ! START: Vibrational Hessian (variational)

    if (Forcefield) then
      nDimTot = max_mOrd+1
    else
      nDimTot = 2*max_mOrd+2
    end if
    l_l = nDimTot

    call mma_allocate(H1,nDimTot,nDimTot,label='H1')
    call mma_allocate(H2,nDimTot,nDimTot,label='H2')
    call mma_allocate(S1,nDimTot,nDimTot,label='S1')
    call mma_allocate(S2,nDimTot,nDimTot,label='S2')

    H1(:,:) = Zero
    S1(:,:) = Zero
    H2(:,:) = Zero
    S2(:,:) = Zero

    ! Calculate inverse mass tensor and its derivatives in r0.
    call mma_allocate(Smat,3*NumOfAt,nOscOld,label='Smat')
    Smat(:,:) = Zero
    !call Int_To_Cart(InterVec,r0,AtCoord,NumOfAt,nOscOld,Mass)
    l_a = NumOfAt
    call Int_To_Cart1(InterVec,r0,AtCoord,l_a,nOsc)
    call CalcS(AtCoord,InterVec,Smat,nOscold,NumOfAt)
    call CalcG(G0,Mass,Smat,nOscold,NumOfAt)
    call mma_deallocate(Smat)
    if (max_term > 2) then
      dh = 1.0e-3_wp
      call CalcGprime(Gprm0,Mass,r0,InterVec,AtCoord,NumOfAt,dh,nOsc)
      dh = 1.0e-2_wp
      call CalcGdblePrime(Gbis0,Mass,r0,InterVec,AtCoord,NumOfAt,dh,nOsc)

      call mma_allocate(T4,nOscOld**4,label='T4')
      call DGEMM_('T','T',nOscOld**2,nOsc,nOscOld,One,Gprm2,nOscOld,BaseInv,nOsc,Zero,T4,nOscOld**2)
      call DGEMM_('T','T',nOsc*nOscOld,nOsc,nOscOld,One,T4,nOscOld,BaseInv,nOsc,Zero,Gprm2,nOsc*nOscOld)
      call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,Gprm2,nOscOld,Base,nOscOld,Zero,T4,nOsc**2)
      call dcopy_(nOsc**3,T4,1,Gprm2,1)

      call DGEMM_('T','T',nOscOld**3,nOsc,nOscOld,One,Gbis2,nOscOld,BaseInv,nOsc,Zero,T4,nOscOld**3)
      call DGEMM_('T','T',nOsc*nOscOld**2,nOsc,nOscOld,One,T4,nOscOld,BaseInv,nOsc,Zero,Gbis2,nOsc*nOscOld**2)
      call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,Gbis2,nOscOld,Base,nOscOld,Zero,T4,nOsc**2*nOscOld)
      call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,T4,nOscOld,Base,nOscOld,Zero,Gbis2,nOsc**3)
      call DGEMM_('T','T',nOscOld**2,nOsc,nOscOld,One,Gprm1,nOscOld,BaseInv,nOsc,Zero,T4,nOscOld**2)
      call DGEMM_('T','T',nOsc*nOscOld,nOsc,nOscOld,One,T4,nOscOld,BaseInv,nOsc,Zero,Gprm1,nOsc*nOscOld)
      call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,Gprm1,nOscOld,Base,nOscOld,Zero,T4,nOsc**2)
      call dcopy_(nOsc**3,T4,1,Gprm1,1)

      call DGEMM_('T','T',nOscOld**3,nOsc,nOscOld,One,Gbis1,nOscOld,BaseInv,nOsc,Zero,T4,nOscOld**3)
      call DGEMM_('T','T',nOsc*nOscOld**2,nOsc,nOscOld,One,T4,nOscOld,BaseInv,nOsc,Zero,Gbis1,nOsc*nOscOld**2)
      call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,Gbis1,nOscOld,Base,nOscOld,Zero,T4,nOsc**2*nOscOld)
      call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,T4,nOscOld,Base,nOscOld,Zero,Gbis1,nOsc**3)
      call DGEMM_('T','T',nOscOld**2,nOsc,nOscOld,One,Gprm0,nOscOld,BaseInv,nOsc,Zero,T4,nOscOld**2)
      call DGEMM_('T','T',nOsc*nOscOld,nOsc,nOscOld,One,T4,nOscOld,BaseInv,nOsc,Zero,Gprm0,nOsc*nOscOld)
      call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,Gprm0,nOscOld,Base,nOscOld,Zero,T4,nOsc**2)
      call dcopy_(nOsc**3,T4,1,Gprm0,1)

      call DGEMM_('T','T',nOscOld**3,nOsc,nOscOld,One,Gbis0,nOscOld,BaseInv,nOsc,Zero,T4,nOscOld**3)
      call DGEMM_('T','T',nOsc*nOscOld**2,nOsc,nOscOld,One,T4,nOscOld,BaseInv,nOsc,Zero,Gbis0,nOsc*nOscOld**2)
      call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,Gbis0,nOscOld,Base,nOscOld,Zero,T4,nOsc**2*nOscOld)
      call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,T4,nOscOld,Base,nOscOld,Zero,Gbis0,nOsc**3)
      call mma_deallocate(T4)
    end if
    call mma_allocate(temp,nOscOld,nOscOld,label='temp')
    call mma_allocate(temp2,nOscOld,nOscOld,label='temp2')

    temp(:,:) = G0
    call unitmat(temp2,nOscOld)
    call Dool_MULA(temp,nOscOld,nOscOld,temp2,nOscOld,nOscOld,det)
    call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,temp2,nOscOld,Base,nOscOld,Zero,temp,nOscOld)
    call mma_deallocate(temp2)
    call mma_allocate(temp1,nOsc,nOsc,label='temp1')
    call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Base,nOscOld,temp,nOscOld,Zero,temp1,nOsc)
    call unitmat(G0,nOsc)
    call Dool_MULA(temp1,nOsc,nOsc,G0,nOsc,nOsc,det)

    call mma_deallocate(temp)
    call mma_deallocate(temp1)

    ! Set up Hamilton matrix for the first state.
    if (Forcefield) then
      call mma_allocate(Base2,nOsc,nOsc,label='Base2')
      call unitmat(Base2,nOsc)

      call SetUpHmat2(energy1,C1,W1,det1,r01,max_mOrd,max_nOrd,max_mInc,max_nInc,mMat,nMat,mInc,nInc,mDec,nDec,H1,S1,Hess1,G1, &
                      Base2,nDimTot,nOsc)
      call SolveSecEq(H1,nDimTot,U1,S1,E1)
      write(u6,'(20f10.1)') ((E1(i)-E1(1))*auTocm,i=2,nOsc)
      call SetUpHmat2(energy1,C2,W2,det2,r01,max_mOrd,max_nOrd,max_nInc,max_nInc,mMat,nMat,mInc,nInc,mDec,nDec,H2,S2,Hess2,G2, &
                      Base2,nDimTot,nOsc)
      call SolveSecEq(H2,nDimTot,U2,S2,E2)
      write(u6,'(20f10.1)') ((E2(i)-E2(1))*auTocm,i=2,nOsc)
    else
      call SetUpHmat(energy1,r1,ipow,var,yin1,trfName1,max_term,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,max_mInc, &
                     max_nInc,max_nInc,mMat,nMat,mInc,nInc,mDec,nDec,H1,S1,G1,G2,G0,Gprm1,Gprm2,Gprm0,Gbis1,Gbis2,Gbis0,det0,Base, &
                     r0,r1,r2,nterm,nvar,ndata,nosc,ndimtot)
      call SolveSort(H1,U1,S1,E1,W1,W2,W1,C1,C2,C1,r01,r02,r01,mInc,mDec,mMat,mdim1,mdim2,OccNumMat1,nOsc,nDimTot)
      call SetUpHmat(energy2,r2,ipow,var,yin2,trfName2,max_term,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd,max_mInc, &
                     max_nInc,max_nInc,mMat,nMat,mInc,nInc,mDec,nDec,H2,S2,G1,G2,G0,Gprm1,Gprm2,Gprm0,Gbis1,Gbis2,Gbis0,det0,Base, &
                     r0,r1,r2,nterm,nvar,ndata,nosc,ndimtot)
      call SolveSort(H2,U2,S2,E2,W1,W2,W2,C1,C2,C2,r01,r02,r02,nInc,nDec,nMat,ndim1,ndim2,OccNumMat2,nOsc,nDimTot)
      call mma_deallocate(yin1)
      call mma_deallocate(yin2)

    end if
    write(u6,*) 'T0',T0*auTocm
    k2 = 1
    l_TermMat_1 = nDimTot-1
    l_TermMat_2 = nDimTot-1

    call mma_allocate(TermMat,[0,l_TermMat_1],[0,l_TermMat_2],label='TermMat')

    do jOrd=0,nDimTot-1
      k1 = 1
      do iOrd=0,nDimTot-1
        TermMat(iOrd,jOrd) = T0+(E2(k2)-E1(k1))
        k1 = k1+1
      end do
      k2 = k2+1
    end do
    call mma_deallocate(H1)
    call mma_deallocate(H2)
    call mma_deallocate(S1)
    call mma_deallocate(S2)
    !call mma_deallocate(G0)
  end if  ! END: Vibrational Hessian (variational)

  call mma_deallocate(AtCoord)
  call mma_deallocate(Hess1)
  call mma_deallocate(Hess2)
  call mma_deallocate(G0)
  call mma_deallocate(G1)
  call mma_deallocate(G2)
  call mma_deallocate(Gprm0)
  call mma_deallocate(Gprm1)
  call mma_deallocate(Gprm2)
  call mma_deallocate(Gbis0)
  call mma_deallocate(Gbis1)
  call mma_deallocate(Gbis2)

  !--------------------------------------------------------------------!

  ! Intensity calculations

  !--------------------------------------------------------------------!

  !PAM: If max_dip=0, there are no transition dipole gradients, but the
  !     array TranDipGradInt must be formally all. anyway:
  if (max_dip == 0) then
    call mma_allocate(TranDipGradInt,3,nOsc,label='TranDipGradInt')
  end if

  if (Matel) then
    l_IntensityMat_1 = ndimtot-1
    l_IntensityMat_2 = ndimtot-1
    l_m_plot = size(m_plot)
    l_n_plot = size(n_plot)
    call mma_allocate(IntensityMat,[0,l_IntensityMat_1],[0,l_IntensityMat_2],label='IntensityMat')
    if (Forcefield) then
      call unitmat(Base2,nOsc)
      call Intensity2(IntensityMat,TermMat,U1,U2,C1,W1,det1,r01,C2,W2,det2,r02,det0,m_max,n_max,max_dip,m_plot,n_plot,TranDip, &
                      TranDipGradInt,Base2,l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nOsc,nDimTot,l_n_plot,l_m_plot)
      call mma_deallocate(TranDipGradInt)
      call mma_deallocate(Base2)
    else
      call Intensity(IntensityMat,TermMat,ipow,var,t_dipin1,t_dipin2,trfName1,U1,U2,C1,W1,det1,r01,C2,W2,det2,r02,det0,m_max, &
                     n_max,max_dip,m_plot,n_plot,r0,r1,r2,Base,l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nOsc, &
                     nDimTot,nPolyTerm,ndata,nvar,MaxNumAt,l_n_plot,l_m_plot)
      call mma_deallocate(t_dipin1)
      call mma_deallocate(t_dipin2)

      call mma_deallocate(var)
      call mma_deallocate(ipow)

    end if
  else
    l_IntensityMat_1 = max_mOrd
    l_IntensityMat_2 = max_nOrd

    call mma_allocate(IntensityMat,[0,l_IntensityMat_1],[0,l_IntensityMat_2],label='IntensityMat')
    call IntForceField(IntensityMat,TermMat,FC00,C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,m_max,n_max,max_dip,Trandip, &
                       TranDipGradInt,harmfreq1,x_anharm1,harmfreq2,x_anharm2,mMat,mInc,mDec,nMat,nInc,nDec,OscStr,max_mOrd, &
                       max_nOrd,nOsc)
    call mma_deallocate(TranDipGradInt)
  end if

  !write results to log.
  write(u6,*) ' Write intensity data to log file.'
  call XFlush(u6)
  call WriteInt(IntensityMat,TermMat,mMat,nMat,OccNumMat2,MatEl,ForceField,E1,E2,T0,harmfreq1,harmfreq2,x_anharm1,x_anharm2, &
                l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nDimTot,nOsc)

  call mma_deallocate(C)
  call mma_deallocate(C1)
  call mma_deallocate(C2)
  call mma_deallocate(W)
  call mma_deallocate(W1)
  call mma_deallocate(W2)

  call mma_deallocate(r00)
  call mma_deallocate(r01)
  call mma_deallocate(r02)

  call mma_deallocate(mMat)
  call mma_deallocate(mInc)
  call mma_deallocate(mDec)
  call mma_deallocate(nMat)
  call mma_deallocate(nInc)
  call mma_deallocate(nDec)

  call mma_deallocate(harmfreq1)
  call mma_deallocate(harmfreq2)

  call mma_deallocate(x_anharm1)
  call mma_deallocate(x_anharm2)
  call mma_deallocate(TermMat)
  call mma_deallocate(IntensityMat)

  call mma_deallocate(U1)
  call mma_deallocate(E1)
  call mma_deallocate(U2)
  call mma_deallocate(E2)

  call mma_deallocate(OccNumMat1)
  call mma_deallocate(OccNumMat2)

  call mma_deallocate(m_plot)
  call mma_deallocate(n_plot)

  call mma_deallocate(NormModes)
  call mma_deallocate(Base)
  call mma_deallocate(BaseInv)

  call mma_deallocate(TranDipGrad)
  call mma_deallocate(r0)
  call mma_deallocate(r1)
  call mma_deallocate(r2)

end if

call mma_deallocate(Mass)
call mma_deallocate(InterVec)
call mma_deallocate(trfName1)
call mma_deallocate(trfName2)

iReturn = 0

end subroutine Mula
