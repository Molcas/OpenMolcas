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
!               2009, Giovanni Ghigo                                   *
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

use Constants, only: Zero, One, Half
use Definitions, only: wp, u5, u6

implicit real*8(a-h,o-z)
#include "inout.fh"
integer InterVec(15*MaxNumAt)
character*80 trfName1(MaxNumAt), trfName2(MaxNumAt)
real*8 Mass(MaxNumAt)
character*4 AtomLbl(MaxNumAt), filnam
integer Bond(2*MaxNumAt)
character*80 Title
real*8 stand_dev, max_err
logical find_minimum
logical lExpan
logical use_weight
logical MatEl, ForceField, lISC, lOldCode, lInCore
logical Cartesian
logical exist
integer nvTabDim
integer iCode, iMaxYes, lNMAT0, lNMAT, lNINC, lNDEC
#include "Constants_mula.fh"
#include "inputdata.fh"
#include "dims.fh"
#include "indims.fh"
#include "WrkSpc.fh"
#include "io_mula.fh"
#include "warnings.h"
integer nIndex(3,0:maxMax_n)

! Initialize.
iReturn = 20
close(u5)
call molcas_open(inpunit,'stdin')
lPotCoef = ip_Dummy
iCode = 00         ! Default: NewCode, InCore
lOldCode = .false. ! X0  CGGn.4
lInCore = .true.   ! 0X  CGGn.4
lNMAT0 = 92        !     CGGn.4
lNMAT = 93         !     CGGn.4
lNINC = 94         !     CGGn.4
lNDEC = 95         !     CGGn.4
iPrint = iPrintLevel(-1)
iPrint_Save = iPrint
!open(unit=inpUnit,file='stdin')

! Read input file and Write header to log file.
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
    write(u6,*) '  - simple harmonic approximation is used.    '
    write(u6,*) '  - m_max set to 0.                           '
    write(u6,*) '  - no plot file.                             '
    if (lOldCode) then
      write(u6,*) '  - Full index matrices evaluation.           '
    else
      write(u6,*) '  - Reduced matrices evaluation.              '
    end if
    if (lInCore) then
      write(u6,*) '  - In-core (memory) algorithm.               '
    else
      write(u6,*) '  - Out-of-core (disk) algorithm.             '
    end if
    write(u6,*) ' ---------------------------------------------'
    write(u6,*)
  end if
  iPrint = 0
  do i=0,maxMax_n
    nIndex(1,i) = 0
    nIndex(2,i) = 0
    nIndex(3,i) = 0
  end do
end if
if ((iPrint >= 1) .and. (Title(1:1) /= ' ')) call WriteHeader(Title)

close(unit=inpUnit)
call GetMem('r01','Allo','Real',ipr01,nOsc)
call GetMem('r02','Allo','Real',ipr02,nOsc)
call GetMem('r1','Allo','Real',ipr1,nOsc)
call GetMem('r2','Allo','Real',ipr2,nOsc)
call GetMem('r0','Allo','Real',ipr0,nOsc)

call GetMem('grad1','Allo','Real',ipgrad1,nOsc)
call GetMem('grad2','Allo','Real',ipgrad2,nOsc)

call GetMem('G1','Allo','Real',ipG1,nOsc*nOsc)
call GetMem('G2','Allo','Real',ipG2,nOsc*nOsc)
call GetMem('G0','Allo','Real',ipG0,nOsc*nOsc)

call GetMem('eigenVec1','Allo','Real',ipeigenVec1,nOsc*nOsc)
call GetMem('eigenVec2','Allo','Real',ipeigenVec2,nOsc*nOsc)

call GetMem('harmfreq1','Allo','Real',ipharmfreq1,nOsc)
call GetMem('harmfreq2','Allo','Real',ipharmfreq2,nOsc)

call GetMem('anharmfreq1','Allo','Real',ipanharmfreq1,nOsc)
call GetMem('anharmfreq2','Allo','Real',ipanharmfreq2,nOsc)

call GetMem('qMat','Allo','Real',ipqMat,3*NumOfAt*nOsc)

call GetMem('PED','Allo','Real',ipPED,nOsc**3)
call GetMem('x_anharm1','Allo','Real',ipx_anharm1,nOsc*nOsc)
call GetMem('x_anharm2','Allo','Real',ipx_anharm2,nOsc*nOsc)

if ((max_term >= 0) .and. (max_term <= 2)) then
  ngdim = 1
else if ((max_term >= 3) .and. (max_term <= 4)) then
  ngdim = nosc
else
  write(u6,*) 'MULA error: max_term=',max_term
  write(u6,*) 'Allowed: 1,2,3, or 4.'
  call Quit_OnUserError()
end if
call getmem('D3_1','Allo','Real',ip_D3_1,ngdim**3)
call getmem('D3_2','Allo','Real',ip_D3_2,ngdim**3)
call getmem('D4_1','Allo','Real',ip_D4_1,ngdim**4)
call getmem('D4_2','Allo','Real',ip_D4_2,ngdim**4)
call getmem('Gprm1','Allo','Real',ip_Gprm1,ngdim**3)
call getmem('Gprm2','Allo','Real',ip_Gprm2,ngdim**3)
call getmem('Gprm0','Allo','Real',ip_Gprm0,ngdim**3)
call getmem('Gbis1','Allo','Real',ip_Gbis1,ngdim**4)
call getmem('Gbis2','Allo','Real',ip_Gbis2,ngdim**4)
call getmem('Gbis0','Allo','Real',ip_Gbis0,ngdim**4)

!----------------------------------------------------------------------!

! Either fit polynomial to energies and find minimum or use
! forcefield and equilibrium geometry given in input.

!----------------------      First State ------------------------------!

call GetMem('AtCoord','Allo','Real',ipAtCoord,3*NumOfAt)

l_a = NumOfAt
l_Hess_1 = l_Hess1
l_Hess_2 = l_Hess1

if (ForceField) then
  call dcopy_(3*NumOfAt,Work(ipAtCoord1),1,Work(ipAtCoord),1)
else
  find_minimum = .true.
  use_weight = .false.
  call GetMem('PotCoef','Allo','Real',lPotCoef,nPolyTerm)
  call PotFit(nPolyTerm,nvar,ndata,iWork(ipipow),Work(ipvar),Work(ipyin1),Work(lPotCoef),Work(ipr01),nOsc,energy1,Work(ipgrad1), &
              Work(ipHess1),work(ip_D3_1),work(ip_D4_1),trfName1,stand_dev,max_err,find_minimum,max_term,use_weight,l_Hess_1, &
              l_Hess_2,ngdim)
  call Int_To_Cart1(InterVec,Work(ipr01),Work(ipAtCoord),l_a,nOsc)
end if
call GetMem('AtCoord1','Free','Real',ipAtCoord1,3*NumOfAt)

! Determine vibrational modes and their frequencies.
!D write(u6,*) ' MULA calling VIBFREQ.'
call VibFreq(Work(ipAtCoord),Work(ipr01),InterVec,Mass,Work(ipHess1),Work(ipG1),work(ip_gprm1),work(ip_gbis1),Work(ipharmfreq1), &
             Work(ipeigenVec1),Work(ipqMat),Work(ipPED),work(ip_D3_1),work(ip_D4_1),Work(ipx_anharm1),Work(ipanharmfreq1), &
             max_term,Cartesian,nOsc,NumOfAt)

!vv use ivv to prevent overoptimization of the code..
ivv = lPotCoef
!call WriteLog(Work(lPotCoef),AtomLbl,Work(ipAtCoord),
if (iPrint >= 1) call WriteLog(Work(ivv),AtomLbl,Work(ipAtCoord),Mass,InterVec,stand_dev,max_err,energy1,Work(ipHess1), &
                               Work(ipG1),Work(ipeigenVec1),Work(ipharmfreq1),Work(ipqMat),Bond,nBond,Work(ipr01),work(ip_D3_1), &
                               work(ip_D4_1),Work(ipPED),Work(ipx_anharm1),Work(ipanharmfreq1),max_term,1,ForceField,NumOfAt,nOsc)
lPotCoef = ivv

!----------------------      Second State -----------------------------!

l_Hess_1 = l_Hess1
l_Hess_2 = l_Hess1

if (ForceField) then
  call dcopy_(3*NumOfAt,Work(ipAtCoord2),1,Work(ipAtCoord),1)
else
  find_minimum = .true.
  use_weight = .false.
  call PotFit(nPolyTerm,nvar,ndata,iWork(ipipow),Work(ipvar),Work(ipyin2),Work(lPotCoef),Work(ipr02),nOsc,energy2,Work(ipgrad2), &
              Work(ipHess2),work(ip_D3_2),work(ip_D4_2),trfName2,stand_dev,max_err,find_minimum,max_term,use_weight,l_Hess_1, &
              l_Hess_2,ngdim)
  call Int_To_Cart1(InterVec,Work(ipr02),Work(ipAtCoord),l_a,nOsc)
end if
call GetMem('AtCoord2','Free','Real',ipAtCoord2,3*NumOfAt)

! Determine vibrational modes and their frequencies.
call VibFreq(Work(ipAtCoord),Work(ipr02),InterVec,Mass,Work(ipHess2),Work(ipG2),work(ip_gprm2),work(ip_gbis2),Work(ipharmfreq2), &
             Work(ipeigenVec2),Work(ipqMat),Work(ipPED),work(ip_D3_2),work(ip_D4_2),Work(ipx_anharm2),Work(ipanharmfreq2), &
             max_term,Cartesian,nOsc,NumOfAt)
if (iPrint >= 1) call WriteLog(Work(lPotCoef),AtomLbl,Work(ipAtCoord),Mass,InterVec,stand_dev,max_err,energy2,Work(ipHess2), &
                               Work(ipG2),Work(ipeigenVec2),Work(ipharmfreq2),Work(ipqMat),Bond,nBond,Work(ipr02),work(ip_D3_2), &
                               work(ip_D4_2),Work(ipPED),Work(ipx_anharm2),Work(ipanharmfreq2),max_term,2,ForceField,NumOfAt,nOsc)

if (.not. Forcefield) call GetMem('PotCoef','Free','Real',lPotCoef,nPolyTerm)
call GetMem('PED','Free','Real',ipPED,nOsc**3)
call GetMem('qMat','Free','Real',ipqMat,3*NumOfAt*nOsc)

call GetMem('grad1','Free','Real',ipgrad1,nOsc)
call GetMem('grad2','Free','Real',ipgrad2,nOsc)

call GetMem('anharmfreq1','Free','Real',ipanharmfreq1,nOsc)
call GetMem('anharmfreq2','Free','Real',ipanharmfreq2,nOsc)

call getmem('D3_1','Free','Real',ip_D3_1,ngdim**3)
call getmem('D3_2','Free','Real',ip_D3_2,ngdim**3)
call getmem('D4_1','Free','Real',ip_D4_1,ngdim**4)
call getmem('D4_2','Free','Real',ip_D4_2,ngdim**4)

!----------------------------------------------------------------------!

! Geometry of intermediate oscillator

!----------------------------------------------------------------------!

T0 = energy2-energy1
call GetMem('C1','Allo','Real',ipC1,nOsc*nOsc)
call GetMem('C2','Allo','Real',ipC2,nOsc*nOsc)
call GetMem('C','Allo','Real',ipC,nOsc*nOsc)
call GetMem('W','Allo','Real',ipW,nOsc*nOsc)
call GetMem('W1','Allo','Real',ipW1,nOsc*nOsc)
call GetMem('W2','Allo','Real',ipW2,nOsc*nOsc)

call GetMem('temp','Allo','Real',iptemp,nOsc*nOsc)

! Calculate W matrices.
do jOsc=1,nOsc
  const1 = One/sqrt(Work(ipharmfreq1+jOsc-1))
  const2 = One/sqrt(Work(ipharmfreq2+jOsc-1))
  do iv=1,nOsc
    Work(ipW1+iv-1+nOsc*(jOsc-1)) = const1*Work(ipeigenVec1+iv-1+nOsc*(jOsc-1))
    Work(ipW2+iv-1+nOsc*(jOsc-1)) = const2*work(ipeigenVec2+iv-1+nOsc*(jOsc-1))
  end do
end do

! Calculate C = W^(-1).
call dcopy_(nOsc**2,[Zero],0,Work(ipC1),1)
call dcopy_(nOsc,[One],0,Work(ipC1),nOsc+1)
call dcopy_(nOsc*nOsc,Work(ipW1),1,Work(iptemp),1)
call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(ipC1),nOsc,nOsc,det)
det1 = abs(One/det)
call dcopy_(nOsc**2,[Zero],0,Work(ipC2),1)
call dcopy_(nOsc,[One],0,Work(ipC2),nOsc+1)
call dcopy_(nOsc*nOsc,Work(ipW2),1,Work(iptemp),1)
call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(ipC2),nOsc,nOsc,det)
det2 = abs(One/det)
call GetMem('temp','Free','Real',iptemp,nOsc*nOsc)

! Calculate the expansion point geometry and save the full geometries for later use
if (iPrint >= 1) call ExpPointHeader()
call GetMem('r00','Allo','Real',ipr00,nOsc)
call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)

call Calc_r00(Work(ipC1),Work(ipC2),Work(ipW1),Work(ipW2),work(ipC),Work(ipW),Work(ipalpha1),Work(ipalpha2),Work(ipr00), &
              Work(ipr01),Work(ipr02),det0,det1,det2,FC00,nOsc)
call dcopy_(nosc,Work(ipr00),1,Work(ipr0),1)
call dcopy_(nosc,Work(ipr01),1,Work(ipr1),1)
call dcopy_(nosc,Work(ipr02),1,Work(ipr2),1)
call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
l_a = NumOfAt
call Int_To_Cart1(InterVec,Work(ipr00),Work(ipAtCoord),l_a,nOsc)
if (iPrint >= 1) call WriteCartCoord(AtomLbl,Work(ipAtCoord),Mass,NumOfAt)
call GetMem('r00','Free','Real',ipr00,nOsc)

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
do iOsc=1,nOsc
  jOsc = 1
  exist = .false.
  do while ((.not. exist) .and. (jOsc <= l_NormModes))
    exist = (iOsc == iWork(ipNormModes+jOsc-1))
    jOsc = jOsc+1
  end do
  if (.not. exist) then
    GE1 = GE1+Half*Work(ipharmfreq1+iOsc-1)
    GE2 = GE2+Half*Work(ipharmfreq2+iOsc-1)
  end if
end do
T0 = GE2-GE1

! Pick out the chosen modes.
nOscOld = nOsc
nOsc = l_NormModes
call GetMem('Lambda','Allo','Real',ipLambda,nOsc)
call GetMem('Base','Allo','Real',ipBase,nOscOld*nOsc)
call GetMem('BaseInv','Allo','Real',ipBaseInv,nOsc*nOscOld)
do jOsc=1,nOsc
  do iOsc=1,nOscOld
    Work(ipBase+iOsc-1+nOscOld*(jOsc-1)) = Work(ipW2+iOsc-1+nOsc*(iWork(ipNormModes+jOsc-1)-1))
    Work(ipBaseInv+jOsc-1+nOsc*(iOsc-1)) = Work(ipC2+iWork(ipNormModes+jOsc-1)-1+nOsc*(iOsc-1))
  end do
end do

! First state.
! Subroutine SolveRedSec(Hess,Gmat,freq,C,W,det)
call GetMem('temp','Allo','Real',iptemp,nOscOld*nOsc)
call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,Work(ipHess1),nOscOld,Work(ipBase),nOscOld,Zero,Work(iptemp),nOscOld)
call GetMem('Hess1','Free','Real',ipHess1,l_Hess1*l_Hess1)

call GetMem('Hess1','Allo','Real',ipHess1,nOsc*nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Work(ipBase),nOscOld,Work(iptemp),nOscOld,Zero,Work(ipHess1),nOsc)
call GetMem('temp2','Allo','Real',iptemp2,nOscOld*nOscOld)
call dcopy_(nOscOld**2,[Zero],0,Work(iptemp2),1)
call dcopy_(nOscOld,[One],0,Work(iptemp2),nOscOld+1)
call Dool_MULA(Work(ipG1),nOsc,nOsc,work(iptemp2),nOscOld,nOscOld,det)
call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,Work(iptemp2),nOscOld,Work(ipBase),nOscOld,Zero,Work(iptemp),nOscOld)
call GetMem('temp2','Free','Real',iptemp2,nOscOld*nOscOld)
call GetMem('temp2','Allo','Real',iptemp2,nOsc*nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Work(ipBase),nOscOld,Work(iptemp),nOscOld,Zero,Work(iptemp2),nOsc)
call dcopy_(nOsc**2,[Zero],0,Work(ipG1),1)
call dcopy_(nOsc,[One],0,Work(ipG1),nOsc+1)
call Dool_MULA(Work(iptemp2),nOscOld,nOscOld,Work(ipG1),nOsc,nOsc,det)
call GetMem('temp2','Free','Real',iptemp2,nOsc*nOsc)
call GetMem('temp','Free','Real',iptemp,nOscOld*nOsc)
call GetMem('temp','Allo','Real',iptemp,nOsc*nOsc)
call SolveSecEq(Work(ipHess1),nOsc,Work(iptemp),Work(ipG1),Work(ipLambda))
do iv=1,nOsc
  Work(ipharmfreq1+iv-1) = sqrt(abs(Work(ipLambda+iv-1)))
end do
if (iPrint >= 1) call WriteFreq(Work(ipharmfreq1),iWork(ipNormModes),l_NormModes,'Frequencies of reduced problem, state 1')
do jOsc=1,nOsc
  const1 = One/sqrt(Work(ipharmfreq1+jOsc-1))
  do iv=1,nOsc
    Work(ipW1+iv-1+nOsc*(jOsc-1)) = const1*Work(iptemp+iv-1+nOsc*(jOsc-1))
  end do
end do
call dcopy_(nOsc**2,[Zero],0,Work(ipC1),1)
call dcopy_(nOsc,[One],0,Work(ipC1),nOsc+1)
call dcopy_(nOsc*nOsc,Work(ipW1),1,Work(iptemp),1)
call Dool_MULA(Work(iptemp),nOscOld,nOsc,Work(ipC1),nOsc,nOsc,det)
det1 = abs(One/det)
call GetMem('temp','Free','Real',iptemp,nOsc*nOsc)

! Second state.
call GetMem('temp','Allo','Real',iptemp,nOscOld*nOsc)
call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,Work(ipHess2),nOscOld,Work(ipBase),nOscOld,Zero,Work(iptemp),nOscOld)

call GetMem('Hess2','Free','Real',ipHess2,l_Hess1*l_Hess1)
call GetMem('Hess2','Allo','Real',ipHess2,nOsc*nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Work(ipBase),nOscOld,Work(iptemp),nOscOld,Zero,Work(ipHess2),nOsc)
call GetMem('temp2','Allo','Real',iptemp2,nOscOld*nOscOld)
call dcopy_(nOscOld**2,[Zero],0,Work(iptemp2),1)
call dcopy_(nOscOld,[One],0,Work(iptemp2),nOscOld+1)
call Dool_MULA(Work(ipG2),nOsc,nOsc,Work(iptemp2),nOscOld,nOscOld,det)
call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,Work(iptemp2),nOscOld,Work(ipBase),nOscOld,Zero,Work(iptemp),nOscOld)
call GetMem('temp2','Free','Real',iptemp2,nOscOld*nOscOld)
call GetMem('temp2','Allo','Real',iptemp2,nOsc*nOsc)
call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Work(ipBase),nOscOld,Work(iptemp),nOscOld,Zero,Work(iptemp2),nOsc)
call dcopy_(nOsc**2,[Zero],0,Work(ipG2),1)
call dcopy_(nOsc,[One],0,Work(ipG2),nOsc+1)
call Dool_MULA(Work(iptemp2),nOsc,nOsc,Work(ipG2),nOsc,nOsc,det)
call GetMem('temp2','Free','Real',iptemp2,nOsc*nOsc)
call GetMem('temp','Free','Real',iptemp,nOscOld*nOsc)
call GetMem('temp','Allo','Real',iptemp,nOsc*nOsc)
call SolveSecEq(Work(ipHess2),nOsc,Work(iptemp),Work(ipG2),Work(ipLambda))
do iv=1,nOsc
  Work(ipharmfreq2+iv-1) = sqrt(abs(Work(ipLambda+iv-1)))
end do
if (iPrint >= 1) call WriteFreq(Work(ipharmfreq2),iWork(ipNormModes),l_NormModes,'Frequencies of reduced problem, state 2')
do jOsc=1,nOsc
  const1 = One/sqrt(Work(ipharmfreq2+jOsc-1))
  do iv=1,nOsc
    Work(ipW2+iv-1+nOsc*(jOsc-1)) = const1*Work(iptemp+iv-1+nOsc*(jOsc-1))
  end do
end do
call dcopy_(nOsc**2,[Zero],0,Work(ipC2),1)
call dcopy_(nOsc,[One],0,Work(ipC2),nOsc+1)
call dcopy_(nOsc*nOsc,Work(ipW2),1,Work(iptemp),1)
call Dool_MULA(Work(iptemp),nOsc,nOsc,Work(ipC2),nOsc,nOsc,det)
det2 = abs(One/det)
call GetMem('temp','Free','Real',iptemp,nOsc*nOsc)

call GetMem('rtemp','Allo','Real',iprtemp,nOscOld)

call dcopy_(nOscOld,Work(ipr01),1,Work(iprtemp),1)
call GetMem('r01','Free','Real',ipr01,nOsc)
call GetMem('r01','Allo','Real',ipr01,nOsc)

call DGEMM_('N','N',nOsc,1,nOscOld,One,Work(ipBaseInv),nOsc,Work(iprtemp),nOscOld,Zero,Work(ipr01),nOsc)
call dcopy_(nOscOld,Work(ipr02),1,Work(iprtemp),1)

call GetMem('r02','Free','Real',ipr02,nOsc)
call GetMem('r02','Allo','Real',ipr02,nOsc)
call DGEMM_('N','N',nOsc,1,nOscOld,One,Work(ipBaseInv),nOsc,Work(iprtemp),nOscOld,Zero,Work(ipr02),nOsc)
call GetMem('rtemp','Free','Real',iprtemp,nOscOld)

call GetMem('r00','Allo','Real',ipr00,nOsc)
call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
call Calc_r00(Work(ipC1),Work(ipC2),Work(ipW1),Work(ipW2),Work(ipC),Work(ipW),Work(ipalpha1),Work(ipalpha2),Work(ipr00), &
              Work(ipr01),Work(ipr02),det0,det1,det2,FC00,nOsc)
call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)

call GetMem('Base2','Allo','Real',ipBase2,nOsc*nOsc)

call dcopy_(nOsc*nOsc,Work(ipW2),1,Work(ipBase2),1)
!call dcopy_(nOsc*nOsc,Work(ipC2),1,Work(ipBaseInv2),1)
!Base2 = W2
!BaseInv2 = C2

!----------------------------------------------------------------------!

!Transform transition dipole gradients from cartesian to internal coordinates

!----------------------------------------------------------------------!

if (Forcefield .and. (max_dip > 0)) then
  call GetMem('Smat','Allo','Real',ipSmat,3*NumOfAt*nOscOld)
  !Smat = Zero
  call dcopy_(3*NumOfAt*nOscOld,[Zero],0,Work(ipSmat),1)
  call CalcS(Work(ipAtCoord),InterVec,Work(ipSmat),nOscOld,NumOfAt)
  call GetMem('Bmat','Allo','Real',ipBmat,3*NumOfAt*nOscOld)

  do j=1,nOscOld
    k = 1
    do i=1,NumOfAt
      Work(ipBmat+k-1+3*NumOfAt*(j-1)) = Work(ipSmat+3*(i-1+NumOfAt*(j-1)))
      Work(ipBmat+k+3*NumOfAt*(j-1)) = Work(ipSmat+1+3*(i-1+NumOfAt*(j-1)))
      Work(ipBmat+k+1+3*NumOfAt*(j-1)) = Work(ipSmat+2+3*(i-1+NumOfAt*(j-1)))
      k = k+3
    end do
  end do
  call GetMem('Smat','Free','Real',ipSmat,3*NumOfAt*nOscOld)
  call GetMem('temp','Allo','Real',iptemp,nOscOld*nOscOld)
  call GetMem('temp3','Allo','Real',iptemp3,nOscOld*nOscOld)

  call DGEMM_('T','N',nOscOld,nOscOld,3*NumOfAt,One,Work(ipBmat),3*NumOfAt,Work(ipBmat),3*NumOfAt,Zero,Work(iptemp),nOscOld)
  call dcopy_(nOscOld**2,[Zero],0,Work(iptemp3),1)
  call dcopy_(nOscOld,[One],0,Work(iptemp3),nOscOld+1)
  call Dool_MULA(Work(iptemp),nOscOld,nOscOld,Work(iptemp3),nOscOld,nOscOld,det)
  call GetMem('temp','Free','Real',iptemp,nOscOld*nOscOld)
  call GetMem('temp1','Allo','Real',iptemp1,3*NumOfAt*nOscOld)

  call DGEMM_('N','N',3*NumOfAt,nOscOld,nOscOld,One,Work(ipBmat),3*NumOfAt,Work(iptemp3),nOscOld,Zero,Work(iptemp1),3*NumOfAt)
  call GetMem('Bmat','Allo','Real',ipBmat,3*NumOfAt*nOscOld)
  call GetMem('temp3','Free','Real',iptemp3,nOscOld*nOscOld)
  call GetMem('temp3','Allo','Real',iptemp3,3*NumOfAt*nOsc)
  call DGEMM_('N','N',3*NumOfAt,nOsc,nOscOld,One,Work(iptemp1),3*NumOfAt,Work(ipBase),nOscOld,Zero,Work(iptemp3),3*NumOfAt)
  call GetMem('temp1','Free','Real',iptemp1,3*NumOfAt*nOscOld)
  call GetMem('temp1','Allo','Real',iptemp1,3*NumOfAt*nOsc)
  call DGEMM_('N','N',3*NumOfAt,nOsc,nOsc,One,Work(iptemp3),3*NumOfAt,Work(ipW),nOsc,Zero,Work(iptemp1),3*NumOfAt)
  call GetMem('temp3','Free','Real',iptemp3,3*NumOfAt*nOsc)
  call GetMem('TranDipGradInt','Allo','Real',ipTranDipGradInt,3*nOsc)
  call DGEMM_('T','N',nOsc,1,3*NumOfAt,One,Work(iptemp1),3*NumOfAt,work(ipTranDipGrad),3*NumOfAt,Zero,Work(ipTranDipGradInt), &
              nOsc)
  call DGEMM_('T','N',nOsc,1,3*NumOfAt,One,Work(iptemp1),3*NumOfAt,work(ipTranDipGrad+1),3*NumOfAt,Zero, &
              Work(ipTranDipGradInt+1),nOsc)
  call DGEMM_('T','N',nOsc,1,3*NumOfAt,One,Work(iptemp1),3*NumOfAt,Work(ipTranDipGrad+2),3*NumOfAt,Zero, &
              Work(ipTranDipGradInt+2),nOsc)
  call GetMem('temp1','Free','Real',iptemp1,3*NumOfAt*nOsc)

  call WriteDip(Work(ipTranDipGradInt),iWork(ipNormModes),'Transition dipole gradient',nOsc)
end if

call GetMem('eigenVec1','Free','Real',ipeigenVec1,nOsc*nOsc)
call GetMem('eigenVec2','Free','Real',ipeigenVec2,nOsc*nOsc)

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
  call GetMem('nMaxQ','Allo','Inte',ipnMaxQ,nOsc)
  call ISC_Rho(iPrint,nOsc,new_n_max,dRho,energy1,energy2,minQ,dMinWind,iWork(ipnMaxQ),Work(ipharmfreq1),Work(ipharmfreq2))
  if (new_n_max < n_max) then
    n_max = new_n_max
    write(u6,*) ' n_max reduced to',n_max
  end if
  if (iPrint >= 3) write(u6,*) ' Actual n_max               =',n_max

  ! Set up excitation matrices

  ! Initial State m
  call TabDim2_drv(m_max,nOsc,nvTabDim)
  mTabDim = nvTabDim-1
  max_mOrd = mTabDim
  call TabDim2_drv(m_max-1,nOsc,nvTabDim)
  max_mInc = nvTabDim-1
  call GetMem('mMat','Allo','Inte',ipmMat,(mTabDim+1)*nOsc)
  call GetMem('mInc','Allo','Inte',ipmInc,(mTabDim+1)*nOsc)
  call GetMem('mDec','Allo','Inte',ipmDec,(mTabDim+1)*nOsc)
  mdim1 = mTabDim
  mdim2 = nOsc
  call MakeTab2(m_max,max_mOrd,max_mInc,mTabDim,iWork(ipmMat),iWork(ipmInc),iWork(ipmDec),nOsc)

  ! Initialize Final State n
  call TabDim2_drv(n_max,nOsc,nvTabDim)
  nTabDim = nvTabDim-1
  max_nOrd = nTabDim
  call TabDim2_drv(n_max-1,nOsc,nvTabDim)
  max_nInc = nvTabDim-1

  ! Memory estimation and algorithm selection
  call GetMem('MULA','MAX','INTE',ipLeft,lLeft)
  if (lInCore) then
    if (lOldCode) then
      iMem = int(7*(nTabDim+1)*nOsc/2)
      if (iMem >= lLeft) then
        write(u6,*) ' Too much memory required (',iMem,' words).'
        write(u6,*) ' Switch to reduced matrices evaluation.'
        lOldCode = .false.
      end if
    else
      iMem = int(3*(nTabDim+1)*nOsc/2)
      if (iMem >= lLeft) then
        write(u6,*) ' Too much memory required (',iMem,' words).'
        write(u6,*) ' Out-of-Core (disk) algorithm will be used.'
        lInCore = .false.
      end if
    end if
  else
    continue
  end if

  call Timing(CPTF0,CPE,TIOTF0,TIOE)
  if (lInCore) then

    if (lOldCode) then
      if (iPrint >= 2) then
        write(u6,*)
        write(u6,*) ' Index matrix evaluation.'
        if (iPrint >= 3) write(u6,*) ' Memory allocated for all index matrix:',iMem,' words,  ',8*iMem/1048576,' MB.'
        call XFlush(u6)
      end if
      call GetMem('nMat','Allo','Inte',ipnMat,(nTabDim+1)*nOsc)
      call GetMem('nInc','Allo','Inte',ipnInc,(nTabDim+1)*nOsc)
      call GetMem('nDec','Allo','Inte',ipnDec,(nTabDim+1)*nOsc)
      ndim1 = nTabDim
      ndim2 = nOsc
      nnsiz = ndim1
      call MakeTab2(n_max,max_nOrd,max_nInc,nTabDim,iWork(ipnMat),iWork(ipnInc),iWork(ipnDec),nOsc)
    else ! NewCode: nInc & nDec reduced
      iMem = (nTabDim+1)*nOsc
      if (iPrint >= 2) then
        write(u6,*)
        write(u6,*) ' First index matrix evaluation.'
        if (iPrint >= 3) write(u6,*) ' Memory allocated for first index matrix:',iMem,' words,  ',(8*iMem+131071)/1048576,' MB.'
        call XFlush(u6)
      end if
      call GetMem('nMat','Allo','Inte',ipnMat,(nTabDim+1)*nOsc)
      call GetMem('Graph2','Allo','Inte',ipGraph2,(n_max+1)*(n_max+1)*nOsc)
      call GetMem('Graph1','Allo','Inte',ipGraph1,(n_max+1)*(nOsc+1))
      ndim1 = nTabDim
      ndim2 = nOsc
      nnsiz = ndim1
      call ISC_MakeTab2(n_max,max_nOrd,max_nInc,nTabDim,iWork(ipnMat),iWork(ipGraph1),iWork(ipGraph2),nOsc)
      call GetMem('Graph1','Free','Inte',ipGraph1,(n_max+1)*(nOsc+1))
    end if
    call Timing(CPTF1,CPE,TIOTF1,TIOE)

    ! Definition of Vector with levels in the window

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Energy level screening.'
      call XFlush(u6)
    end if
    call GetMem('lVec','Allo','Inte',iplVec,(nTabDim+1))
    call LogEVec(iPrint,nOsc,max_nOrd,minQ,iWork(ipnMaxQ),iWork(ipnMat),iWork(iplVec),nYes)
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
        write(u6,*) ' Memory allocated for energy matrix:',(nTabDim+1),' words,  ',(8*(nTabDim+1)+131071)/1048576,' MB.'
      end if
      call XFlush(u6)
    end if
    call GetMem('lTVec','Allo','Inte',iplTVec,(nTabDim+1))
    call GetMem('VibLevel2','Allo','Real',ipVibLevel2,max_nOrd+1)
    call ISC_Ene(iPrint,nOsc,max_nOrd,nYes,iWork(ipnMat),nTabDim,GE1,GE2,Work(ipharmfreq1),Work(ipharmfreq2),Work(ipx_anharm1), &
                 Work(ipx_anharm2),dMinWind,dRho,Work(ipVibLevel2),iWork(iplVec),iWork(iplTVec))
    call GetMem('VibLevel2','Free','Real',ipVibLevel2,max_nOrd+1)
    call GetMem('lTVec','Free','Inte',iplTVec,(nTabDim+1))
    if (nYes <= 0) then
      write(u6,*)
      write(u6,*) ' ************ ERROR *************'
      write(u6,*) ' No energy levels in the Window !'
      write(u6,*) ' ********************************'
      write(u6,*)
      if (iPrint < 4) call Quit_OnUserError()
    end if

    ! The State Window

    call GetMem('VibWind2','Allo','Inte',ipVibWind2,nYes)
    call MkVibWind2(iPrint,nYes,iMaxYes,max_nOrd,iWork(iplVec),iWork(ipVibWind2))
    call GetMem('lVec','Free','INTE',iplVec,(nTabDim+1))
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
        if (iPrint >= 3) write(u6,*) ' Memory allocated for remaining index matrix:',iMem,' words,  ',(8*iMem+131071)/1048576,' MB.'
        call XFlush(u6)
      end if
      call GetMem('nInc','Allo','Inte',ipnInc,(iMaxYes+1)*nOsc)
      call GetMem('nDec','Allo','Inte',ipnDec,(iMaxYes+1)*nOsc)
      call Mk_nIncDec(n_max,nTabDim,iMaxYes,iWork(ipnInc),iWork(ipnDec),iWork(ipnMat),iWork(ipGraph2),nOsc)
      call GetMem('Graph2','Free','Inte',ipGraph2,(n_max+1)*(m_max+1)*nOsc)
    end if
    call Timing(CPTF5,CPE,TIOTF5,TIOE)

    ! The ISC rate

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Franck-Condon factors evaluation.'
      call XFlush(u6)
    end if
    call GetMem('FCWind2','Allo','Real',ipFCWind2,nYes)
    call ISC_Rate(iPrint,nOsc,max_nOrd,iMx_nOrd,iMaxYes,nYes,dMinWind,iWork(ipVibWind2),Work(ipC1),Work(ipC2),Work(ipW1), &
                  Work(ipW2),det0,det1,det2,Work(ipC),Work(ipW),Work(ipr01),Work(ipr02),Work(ipr00),mTabDim,iWork(ipmMat),nTabDim, &
                  iWork(ipnMat),iWork(ipmInc),iWork(ipmDec),iWork(ipnInc),iWork(ipnDec),m_max,n_max,max_dip,nnsiz,FC00, &
                  Work(ipFCWind2),dRho)
    call Timing(CPTF6,CPE,TIOTF6,TIOE)

    call GetMem('FCWind2','Free','Real',ipFCWind2,nYes)
    call GetMem('VibWind2','Free','Inte',ipVibWind2,nYes)

    if (lOldCode) then
      call GetMem('nDec','Free','Inte',ipnDec,(nTabDim+1)*nOsc)
      call GetMem('nInc','Free','Inte',ipnInc,(nTabDim+1)*nOsc)
    else
      call GetMem('nDec','Free','Inte',ipnDec,(nYes+1)*nOsc)
      call GetMem('nInc','Free','Inte',ipnInc,(nYes+1)*nOsc)
    end if
    call GetMem('nMat','Free','Inte',ipnMat,(nTabDim+1)*nOsc)

  else ! Out-of-core.

    !GGt --- Out-of-core ---------------------------------------------------

    ! Open files

    filnam = 'MAT0'
    call DaName_mf_wa(lNMAT0,filnam)

    call GetMem('Graph2','Allo','Inte',ipGraph2,(n_max+1)*(n_max+1)*nOsc)
    call GetMem('Graph1','Allo','Inte',ipGraph1,(n_max+1)*(nOsc+1))
    call ISCD_MakeGraphs(n_max,max_nOrd,max_nInc,iWork(ipGraph1),iWork(ipGraph2),nOsc)
    call GetMem('Graph1','Free','Inte',ipGraph1,(n_max+1)*(nOsc+1))
    call GetMem('nTabDim','Allo','Inte',ipnTabDim,nTabDim+1)
    call GetMem('nMat0','Allo','Inte',ipnMat0,nOsc)
    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Index matrix evaluation.'
      if (iPrint >= 3) then
        write(u6,*) ' Memory allocated for address array:',(nTabDim+1),' words,  ',(8*(nTabDim+1)+131071)/1048576,' MB.'
      end if
      call XFlush(u6)
    end if
    call ISCD_MakenMat(n_max,nOsc,lNMAT0,nTabDim,iWork(ipGraph2),iWork(ipnTabDim),iWork(ipnMat0))
    call Timing(CPTF1,CPE,TIOTF1,TIOE)

    ! Definition of Vector with levels in the window

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Energy level screening.'
      call XFlush(u6)
    end if
    call GetMem('lVec','Allo','Inte',iplVec,(nTabDim+1))
    call ISCD_LogEVec(iPrint,nOsc,max_nOrd,minQ,nYes,lNMAT0,nTabDim,iWork(ipnTabDim),iWork(ipnMaxQ),iWork(ipnMat0),iWork(iplVec))
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
        write(u6,*) ' Memory allocated for energy matrix:',(nTabDim+1),' words,  ',(8*(nTabDim+1)+131071)/1048576,' MB.'
      end if
      call XFlush(u6)
    end if
    call GetMem('lTVec','Allo','Inte',iplTVec,(nTabDim+1))
    call GetMem('VibLevel2','Allo','Real',ipVibLevel2,max_nOrd+1)
    call ISCD_Ene(iPrint,nOsc,max_nOrd,nYes,lNMAT0,nTabDim,GE1,GE2,Work(ipharmfreq1),Work(ipharmfreq2),Work(ipx_anharm1), &
                  Work(ipx_anharm2),dMinWind,dRho,iWork(ipnMat0),iWork(ipnTabDim),iWork(iplVec),iWork(iplTVec),Work(ipVibLevel2))
    call GetMem('VibLevel2','Free','Real',ipVibLevel2,max_nOrd+1)
    call GetMem('lTVec','Free','Inte',iplTVec,(nTabDim+1))
    if (nYes <= 0) then
      write(u6,*)
      write(u6,*) ' ************ ERROR *************'
      write(u6,*) ' No energy levels in the Window !'
      write(u6,*) ' ********************************'
      write(u6,*)
      if (iPrint < 4) call Quit_OnUserError()
    end if

    ! The State Window

    call GetMem('VibWind2','Allo','Inte',ipVibWind2,nYes)
    call MkVibWind2(iPrint,nYes,iMaxYes,max_nOrd,iWork(iplVec),iWork(ipVibWind2))
    call GetMem('lVec','Free','Inte',iplVec,(nTabDim+1))
    call Timing(CPTF3,CPE,TIOTF3,TIOE)

    ! Memory estimation & nMat transfer

    if (lOldCode) then
      iMaxYes = nTabDim
    end if
    call GetMem('MULA','MAX','INTE',ipLeft,lLeft)
    if (iPrint >= 3) then
      write(u6,'(A,I10,A,I4,A)') '  Available memory                          :',lLeft,' words,  ',(lLeft+131071)/131072,' MB.'
    end if
    iUMem = int(1*(nTabDim+1))  ! Memory for U, 1=> No Hot states
    iRateMem = iUMem+nTabDim+2+13*nOsc*nOsc+5*nOsc
    iRateMem = int(2.1_wp*iRateMem) ! to convert from INTE to REAL
    lLeft = lLeft-iRateMem
    if (iRateMem < 0) then ! .or. (iRateMem > 2048*131072)) then
      write(u6,'(A,I10,A,I4,A)') '  Estimated memory for FC factors evaluation:',iRateMem,' words,  ',(iRateMem+131071)/131072, &
                                 ' MB.'
      write(u6,'(A,I10,A,I4,A)') '  Memory left for batches                   :',lLeft,' words,  ',(lLeft+131071)/131072,' MB.'
      write(u6,*)
      write(u6,*) ' ****************** ERROR ********************'
      write(u6,*) ' Not enough memory for FC factors evaluation !'
      write(u6,*) ' *********************************************'
      write(u6,*)
      if (iPrint >= 3) call GetMem('MULA','LIST','INTE',iDum,iDum)
      call Quit_OnUserError()
    end if
    lMBatch = int(lLeft/3)
    lMBatch = min(lMBatch,lMaxMBatch*131072) ! Max 512 MB
    lMBatch = max(lMBatch,lMinMBatch*131072) ! Min  64 MB
    !GGt lMBatch = 131072*8  ! Test only !!!!
    lBatch = int(lMBatch/nOsc)
    if (lBatch > (iMaxYes+1)) lBatch = (iMaxYes+1)
    lMBatch = lBatch*nOsc
    nBatch = int((iMaxYes+1)/lBatch)
    leftBatch = (iMaxYes+1)-nBatch*lBatch
    if (nBatch+1 > maxMax_n) then
      write(u6,*)
      write(u6,*) ' ******************** ERROR *******************'
      write(u6,*) ' Not enough number of batches for Out-of-core !'
      write(u6,*) ' Increase maxMax_n in src/mula/io_mula.fh      '
      write(u6,*) ' **********************************************'
      write(u6,*)
      call Quit_OnUserError()
    end if
    if (iPrint >= 2) then
      write(u6,*)
      if (iPrint >= 3) then
        write(u6,'(A,I10,A,I4,A)') '  Estimated memory for FC factors evaluation:',iRateMem,' words,  ',(iRateMem+131071)/131072, &
                                   ' MB.'
        write(u6,'(A,I10,A,I4,A)') '  Memory left for batches                   :',lLeft,' words,  ',(lLeft+131071)/131072,' MB.'
        write(u6,'(A,I10,A,I4,A)') '  Memory allocated for batch                :',lMBatch,' words,  ',(lMBatch+131071)/131072, &
                                   ' MB.'
        write(u6,'(A,I8)') '  * Elements for batch:',lBatch
        write(u6,'(A,I8,A,I8,A)') '  * Number of batches :',nBatch,' (',lBatch*nBatch,' elements)'
        write(u6,'(A,I8)') '  * Residual elements :',leftBatch
      end if
      write(u6,*) ' Index matrix reloading.'
      call XFlush(u6)
    end if

    filnam = 'NMAT'
    call DaName_mf_wa(lNMAT,filnam)
    call GetMem('nMAT','Allo','Inte',ipnMAT,nOsc*lBatch)
    call ISCD_ReloadNMAT(nTabDim,iMaxYes,nOsc,lNMAT0,lNMAT,lBatch,nBatch,leftBatch,nIndex,iWork(ipnTabDim),iWork(ipnMAT0), &
                         iWork(ipnMAT))
    call GetMem('nMat0','Free','Inte',ipnMat0,nOsc)
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
    call GetMem('nInc','Allo','Inte',ipnInc,nOsc*lBatch)
    call GetMem('nDec','Allo','Inte',ipnDec,nOsc*lBatch)
    call ISCD_MakenIncDec(n_max,iMaxYes,nOsc,lNMAT,lNINC,lNDEC,lBatch,nBatch,leftBatch,nIndex,iWork(ipGraph2),iWork(ipnMAT), &
                          iWork(ipnInc),iWork(ipnDec))
    call GetMem('Graph2','Free','Inte',ipGraph2,(n_max+1)*(m_max+1)*nOsc)
    call Timing(CPTF5,CPE,TIOTF5,TIOE)

    ! The ISC rate

    if (iPrint >= 2) then
      write(u6,*)
      write(u6,*) ' Franck-Condon factors evaluation.'
      call XFlush(u6)
    end if
    call GetMem('FCWind2','Allo','Real',ipFCWind2,nYes)
    call ISCD_Rate(iPrint,nOsc,max_nOrd,iMx_nOrd,iMaxYes,nYes,dMinWind,lBatch,nBatch,leftBatch,nIndex,iWork(ipVibWind2),lNMAT0, &
                   lNMAT,lNINC,lNDEC,nTabDim,iWork(ipnTabDim),Work(ipC1),Work(ipC2),Work(ipW1),Work(ipW2),det0,det1,det2, &
                   Work(ipC),Work(ipW),Work(ipr01),Work(ipr02),Work(ipr00),m_max,n_max,max_dip,nnsiz,FC00,Work(ipFCWind2),dRho, &
                   mTabDim,iWork(ipmMat),iWork(ipmInc),iWork(ipmDec),iWork(ipnMat),iWork(ipnInc),iWork(ipnDec))
    call Timing(CPTF6,CPE,TIOTF6,TIOE)

    call GetMem('FCWind2','Free','Real',ipFCWind2,nYes)
    call GetMem('VibWind2','Free','Inte',ipVibWind2,nYes)
    call GetMem('nDec','Free','Inte',ipnDec,nOsc*lBatch)
    call GetMem('nInc','Free','Inte',ipnInc,nOsc*lBatch)
    call GetMem('nMat','Free','Inte',ipnMAT,nOsc*lBatch)
    call GetMem('nTabDim','Free','Inte',ipnTabDim,nTabDim+1)

    call DaClos(lNDEC)
    call DaClos(lNINC)
    call DaClos(lNMAT)
    call DaClos(lNMAT0)

  end if

  if (iPrint >= 3) then
    write(u6,*)
    write(u6,'(A)') ' Timing informations (sec.):              '
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

  call GetMem('mDec','Free','Inte',ipmDec,(mTabDim+1)*nOsc)
  call GetMem('mInc','Free','Inte',ipmInc,(mTabDim+1)*nOsc)
  call GetMem('mMat','Free','Inte',ipmMat,(mTabDim+1)*nOsc)

  call GetMem('nMaxQ','Free','Inte',ipnMaxQ,nOsc)
  call GetMem('Base2','Free','Real',ipBase2,nOsc*nOsc)

  call GetMem('r00','Free','Real',ipr00,nOsc)
  call GetMem('r02','Free','Real',ipr02,nOsc)
  call GetMem('r01','Free','Real',ipr01,nOsc)

  call GetMem('Hess2','Free','Real',ipHess2,nOsc*nOsc)
  call GetMem('Hess1','Free','Real',ipHess1,nOsc*nOsc)

  call GetMem('BaseInv','Free','Real',ipBaseInv,nOsc*nOscOld)
  call GetMem('Base','Free','Real',ipBase,nOscOld*nOsc)
  call GetMem('Lambda','Free','Real',ipLambda,nOsc)

  call GetMem('W2','Free','Real',ipW2,nOsc*nOsc)
  call GetMem('W1','Free','Real',ipW1,nOsc*nOsc)
  call GetMem('W','Free','Real',ipW,nOsc*nOsc)
  call GetMem('C','Free','Real',ipC,nOsc*nOsc)
  call GetMem('C2','Free','Real',ipC2,nOsc*nOsc)
  call GetMem('C1','Free','Real',ipC1,nOsc*nOsc)

  call GetMem('AtCoord','Free','Real',ipAtCoord,3*NumOfAt)
  call Getmem('Gbis0','Free','Real',ip_Gbis0,ngdim**4)
  call Getmem('Gbis2','Free','Real',ip_Gbis2,ngdim**4)
  call Getmem('Gbis1','Free','Real',ip_Gbis1,ngdim**4)

  call Getmem('Gprm0','Free','Real',ip_Gprm0,ngdim**3)
  call Getmem('Gprm2','Free','Real',ip_Gprm2,ngdim**3)
  call Getmem('Gprm1','Free','Real',ip_Gprm1,ngdim**3)

  call GetMem('x_anharm2','Free','Real',ipx_anharm2,nOsc*nOsc)
  call GetMem('x_anharm1','Free','Real',ipx_anharm1,nOsc*nOsc)
  call GetMem('harmfreq2','Free','Real',ipharmfreq2,nOsc)
  call GetMem('harmfreq1','Free','Real',ipharmfreq1,nOsc)

  call GetMem('G0','Free','Real',ipG0,nOsc*nOsc)
  call GetMem('G2','Free','Real',ipG2,nOsc*nOsc)
  call GetMem('G1','Free','Real',ipG1,nOsc*nOsc)

  call GetMem('r0','Free','Real',ipr0,nOsc)
  call GetMem('r2','Free','Real',ipr2,nOsc)
  call GetMem('r1','Free','Real',ipr1,nOsc)

  call GetMem('n_plot','Free','Inte',ipn_plot,l_n_plot)
  call GetMem('m_plot','Free','Inte',ipm_plot,l_m_plot)

  call GetMem('NormModes','Free','Inte',ipNormModes,l_NormModes)
  call GetMem('TranDipGrad','Free','Real',ipTranDipGrad,3*NumOfAt)

  lISC = .false.
  !call GetMem('Test_F','LIST','INTE',iDum,iDum)
  goto 999

end if
!GGn -------------------------------------------------------------------

!----------------------------------------------------------------------!

! Set up excitation matrices

!----------------------------------------------------------------------!

! Calculate dimensions given max level of excitation for the different states.
call TabDim2_drv(m_max,nOsc,nvTabDim)
mTabDim = nvTabDim-1
call TabDim2_drv(n_max,nOsc,nvTabDim)
nTabDim = nvTabDim-1

! Set up mMat for L.
max_mOrd = mTabDim
call TabDim2_drv(m_max-1,nOsc,nvTabDim)
max_mInc = nvTabDim-1
call GetMem('mMat','Allo','Inte',ipmMat,(mTabDim+1)*nOsc)

call GetMem('mInc','Allo','Inte',ipmInc,(mTabDim+1)*nOsc)
call GetMem('mDec','Allo','Inte',ipmDec,(mTabDim+1)*nOsc)
! Put dimensions into common block:
mdim1 = mTabDim
mdim2 = nOsc
call MakeTab2(m_max,max_mOrd,max_mInc,mTabDim,iWork(ipmMat),iWork(ipmInc),iWork(ipmDec),nOsc)

! Set up nMat for U.
max_nOrd = nTabDim
call TabDim2_drv(n_max-1,nOsc,nvTabDim)
max_nInc = nvTabDim-1
call GetMem('nMat','Allo','Inte',ipnMat,(nTabDim+1)*nOsc)

call GetMem('nInc','Allo','Inte',ipnInc,(nTabDim+1)*nOsc)
call GetMem('nDec','Allo','Inte',ipnDec,(nTabDim+1)*nOsc)
! Put dimensions into common block:
ndim1 = nTabDim
ndim2 = nOsc
nnsiz = ndim1
call MakeTab2(n_max,max_nOrd,max_nInc,nTabDim,iWork(ipnMat),iWork(ipnInc),iWork(ipnDec),nOsc)

if (.not. MatEl) then
  l_TermMat_1 = max_mOrd
  l_TermMat_2 = max_nOrd
  call GetMem('TermMat','Allo','Real',ipTermMat,(l_TermMat_1+1)*(l_TermMat_2+1))
  call GetMem('level1','Allo','Inte',iplevel1,nOsc)
  call GetMem('level2','Allo','Inte',iplevel2,nOsc)
  do jOrd=0,max_nOrd
    do iv=1,nOsc
      iWork(iplevel2+iv-1) = iWork(ipnMat+jOrd+(nTabDim+1)*(iv-1))
    end do
    do iOrd=0,max_mOrd
      do iv=1,nOsc
        iWork(iplevel1+iv-1) = iWork(ipmMat+iOrd+(mTabDim+1)*(iv-1))
      end do
      l_harm = nOsc
      call TransEnergy(GE1,Work(ipx_anharm1),Work(ipharmfreq1),iWork(iplevel1),GE2,Work(ipx_anharm2),Work(ipharmfreq2), &
                       iWork(iplevel2),Work(ipTermMat+iOrd+(l_TermMat_1+1)*jOrd),l_harm)
    end do
  end do
  call GetMem('level1','Free','Inte',iplevel1,nOsc)
  call GetMem('level2','Free','Inte',iplevel2,nOsc)

end if

!----------------------------------------------------------------------!

! Calculate vibrational Hessian if a variational calculation is chosen

!----------------------------------------------------------------------!

l_l = 1
call GetMem('U1','Allo','Real',ipU1,l_l*l_l)
call GetMem('E1','Allo','Real',ipE1,l_l)
call GetMem('U2','Allo','Real',ipU2,l_l*l_l)
call GetMem('E2','Allo','Real',ipE2,l_l)
!l_l = 1
call GetMem('OccNumMat1','Allo','Real',ipOccNumMat1,l_l*3)
call GetMem('OccNumMat2','Allo','Real',ipOccNumMat2,l_l*3)
if (MatEl) then  ! START: Vibrational Hessian (variational)

  call GetMem('U2','Free','Real',ipU2,l_l*l_l)
  call GetMem('E2','Free','Real',ipE2,l_l)
  call GetMem('U1','Free','Real',ipU1,l_l*l_l)
  call GetMem('E1','Free','Real',ipE1,l_l)

  call GetMem('OccNumMat1','Free','Real',ipOccNumMat1,l_l*3)
  call GetMem('OccNumMat2','Free','Real',ipOccNumMat2,l_l*3)

  if (Forcefield) then
    nDimTot = max_mOrd+1
  else
    nDimTot = 2*max_mOrd+2
  end if
  l_l = nDimTot
  call GetMem('U1','Allo','Real',ipU1,l_l*l_l)
  call GetMem('E1','Allo','Real',ipE1,l_l)
  call GetMem('U2','Allo','Real',ipU2,l_l*l_l)
  call GetMem('E2','Allo','Real',ipE2,l_l)

  call GetMem('H1','Allo','Real',ipH1,nDimTot*nDimTot)
  call GetMem('H2','Allo','Real',ipH2,nDimTot*nDimTot)
  call GetMem('S1','Allo','Real',ipS1,nDimTot*nDimTot)
  call GetMem('S2','Allo','Real',ipS2,nDimTot*nDimTot)
  call GetMem('OccNumMat1','Allo','Real',ipOccNumMat1,l_l*3)
  call GetMem('OccNumMat2','Allo','Real',ipOccNumMat2,l_l*3)

  !H1 = Zero
  call dcopy_(nDimTot*nDimTot,[Zero],0,Work(ipH1),1)
  call dcopy_(nDimTot*nDimTot,[Zero],0,Work(ipS1),1)
  call dcopy_(nDimTot*nDimTot,[Zero],0,Work(ipH2),1)
  call dcopy_(nDimTot*nDimTot,[Zero],0,Work(ipS2),1)
  !S1 = Zero
  !H2 = Zero
  !S2 = Zero

  ! Calculate inverse mass tensor and its derivatives in r0.
  call GetMem('Smat','Allo','Real',ipSmat,3*NumOfAt*nOscOld)
  call dcopy_(3*NumOfAt*nOscOld,[Zero],0,Work(ipSmat),1)
  !call Int_To_Cart(InterVec,r0,AtCoord,NumOfAt,nOscOld,Mass)
  l_a = NumOfAt
  call Int_To_Cart1(InterVec,Work(ipr0),Work(ipAtCoord),l_a,nOsc)
  call CalcS(Work(ipAtCoord),InterVec,Work(ipSmat),nOscold,NumOfAt)
  call CalcG(Work(ipG0),Mass,Work(ipSmat),nOscold,NumOfAt)
  if (max_term > 2) then
    dh = 1.0e-3_wp
    call CalcGprime(work(ip_gprm0),Mass,Work(ipr0),InterVec,Work(ipAtCoord),NumOfAt,dh,nOsc)
    dh = 1.0e-2_wp
    call CalcGdblePrime(work(ip_gbis0),Mass,Work(ipr0),InterVec,Work(ipAtCoord),NumOfAt,dh,nOsc)

    call GetMem('T4','Allo','Real',ipT4,nOscOld**4)
    call DGEMM_('T','T',nOscOld**2,nOsc,nOscOld,One,work(ip_gprm2),nOscOld,Work(ipBaseInv),nOsc,Zero,Work(ipT4),nOscOld**2)
    call DGEMM_('T','T',nOsc*nOscOld,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBaseInv),nOsc,Zero,work(ip_gprm2),nOsc*nOscOld)
    call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,work(ip_gprm2),nOscOld,Work(ipBase),nOscOld,Zero,Work(ipT4),nOsc**2)
    call getmem('Gprm2','Free','Real',ip_Gprm2,ngdim**3)

    call getmem('Gprm2','Allo','Real',ip_Gprm2,ngdim**3)
    call dcopy_(nOsc**3,Work(ipT4),1,work(ip_gprm2),1)

    call DGEMM_('T','T',nOscOld**3,nOsc,nOscOld,One,work(ip_gbis2),nOscOld,Work(ipBaseInv),nOsc,Zero,Work(ipT4),nOscOld**3)
    call DGEMM_('T','T',nOsc*nOscOld**2,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBaseInv),nOsc,Zero,work(ip_gbis2), &
                nOsc*nOscOld**2)
    call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,work(ip_gbis2),nOscOld,Work(ipBase),nOscOld,Zero,Work(ipT4), &
                nOsc**2*nOscOld)
    call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBase),nOscOld,Zero,work(ip_gbis2),nOsc**3)
    call DGEMM_('T','T',nOscOld**2,nOsc,nOscOld,One,work(ip_gprm1),nOscOld,Work(ipBaseInv),nOsc,Zero,Work(ipT4),nOscOld**2)
    call DGEMM_('T','T',nOsc*nOscOld,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBaseInv),nOsc,Zero,work(ip_gprm1),nOsc*nOscOld)
    call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,work(ip_gprm1),nOscOld,Work(ipBase),nOscOld,Zero,Work(ipT4),nOsc**2)
    call getmem('Gprm1','Free','Real',ip_Gprm1,ngdim**3)
    call getmem('Gprm1','Allo','Real',ip_Gprm1,ngdim**3)
    call dcopy_(nOsc**3,Work(ipT4),1,work(ip_gprm1),1)

    call DGEMM_('T','T',nOscOld**3,nOsc,nOscOld,One,work(ip_gbis1),nOscOld,Work(ipBaseInv),nOsc,Zero,Work(ipT4),nOscOld**3)
    call DGEMM_('T','T',nOsc*nOscOld**2,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBaseInv),nOsc,Zero,work(ip_gbis1), &
                nOsc*nOscOld**2)
    call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,work(ip_gbis1),nOscOld,Work(ipBase),nOscOld,Zero,Work(ipT4), &
                nOsc**2*nOscOld)
    call getmem('Gbis1','Free','Real',ip_Gbis1,ngdim**4)
    call getmem('Gbis1','Allo','Real',ip_Gbis1,ngdim**4)
    call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBase),nOscOld,Zero,work(ip_gbis1),nOsc**3)
    call DGEMM_('T','T',nOscOld**2,nOsc,nOscOld,One,work(ip_gprm0),nOscOld,Work(ipBaseInv),nOsc,Zero,Work(ipT4),nOscOld**2)
    call DGEMM_('T','T',nOsc*nOscOld,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBaseInv),nOsc,Zero,work(ip_gprm0),nOsc*nOscOld)
    call DGEMM_('T','N',nOsc**2,nOsc,nOscOld,One,work(ip_gprm0),nOscOld,Work(ipBase),nOscOld,Zero,Work(ipT4),nOsc**2)
    call getmem('Gprm0','Free','Real',ip_Gprm0,ngdim**3)
    ngdim = nosc
    call getmem('Gprm0','Allo','Real',ip_Gprm0,ngdim**3)
    call dcopy_(nOsc**3,Work(ipT4),1,work(ip_gprm0),1)

    call DGEMM_('T','T',nOscOld**3,nOsc,nOscOld,One,work(ip_gbis0),nOscOld,Work(ipBaseInv),nOsc,Zero,Work(ipT4),nOscOld**3)
    call DGEMM_('T','T',nOsc*nOscOld**2,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBaseInv),nOsc,Zero,work(ip_gbis0), &
                nOsc*nOscOld**2)
    call DGEMM_('T','N',nOsc**2*nOscOld,nOsc,nOscOld,One,work(ip_gbis0),nOscOld,Work(ipBase),nOscOld,Zero,Work(ipT4), &
                nOsc**2*nOscOld)
    call getmem('Gbis0','Free','Real',ip_Gbis0,ngdim**4)
    call getmem('Gbis0','Allo','Real',ip_Gbis0,ngdim**4)
    call DGEMM_('T','N',nOsc**3,nOsc,nOscOld,One,Work(ipT4),nOscOld,Work(ipBase),nOscOld,Zero,work(ip_gbis0),nOsc**3)
    call GetMem('T4','Free','Real',ipT4,nOscOld**4)
  end if
  call GetMem('temp','Allo','Real',iptemp,nOscOld*nOscOld)
  call GetMem('temp2','Allo','Real',iptemp2,nOscOld*nOscOld)

  call dcopy_(nOscOld*nOscOld,Work(ipG0),1,Work(iptemp),1)
  call dcopy_(nOscOld**2,[Zero],0,Work(iptemp2),1)
  call dcopy_(nOscOld,[One],0,Work(iptemp2),nOscOld+1)
  call Dool_MULA(Work(iptemp),nOscOld,nOscOld,Work(iptemp2),nOscOld,nOscOld,det)
  call DGEMM_('N','N',nOscOld,nOsc,nOscOld,One,Work(iptemp2),nOscOld,Work(ipBase),nOscOld,Zero,Work(iptemp),nOscOld)
  call GetMem('temp2','Free','Real',iptemp2,nOscOld*nOscOld)
  call GetMem('temp2','Allo','Real',iptemp2,nOsc*nOsc)
  call DGEMM_('T','N',nOsc,nOsc,nOscOld,One,Work(ipBase),nOscOld,Work(iptemp),nOscOld,Zero,Work(iptemp2),nOsc)
  call dcopy_(nOsc**2,[Zero],0,Work(ipG0),1)
  call dcopy_(nOsc,[One],0,Work(ipG0),nOsc+1)
  call Dool_MULA(Work(iptemp2),nOsc,nOsc,Work(ipG0),nOsc,nOsc,det)

  call GetMem('temp','Free','Real',iptemp,nOscOld*nOscOld)
  call GetMem('temp2','Free','Real',iptemp2,nOsc*nOsc)

  call GetMem('Smat','Free','Real',ipSmat,3*NumOfAt*nOscOld)

  ! Set up Hamilton matrix for the first state.
  if (Forcefield) then
    !Base2 = Zero
    call dcopy_(nOsc*nOsc,[Zero],0,Work(ipBase2),1)
    do i=1,nOsc
      Work(ipBase2+i-1+nOsc*(i-1)) = One
    end do

    call SetUpHmat2(energy1,energy2,Work(ipC1),Work(ipW1),det1,Work(ipr01),Work(ipr02),max_mOrd,max_nOrd,max_nOrd,max_mInc, &
                    max_nInc,max_nInc,iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),iWork(ipnInc),iWork(ipmDec),iWork(ipnDec), &
                    Work(ipH1),Work(ipS1),Work(ipHess1),Work(ipG1),Work(ipBase2),Work(ipr01),nnsiz,nDimTot,nOsc)
    call SolveSecEq(Work(ipH1),nDimTot,Work(ipU1),Work(ipS1),Work(ipE1))
    write(u6,'(20f10.1)') ((Work(ipE1+i-1)-Work(ipE1))*auTocm,i=2,nOsc)
    call SetUpHmat2(energy1,energy2,Work(ipC2),Work(ipW2),det2,Work(ipr01),Work(ipr02),max_mOrd,max_nOrd,max_nOrd,max_mInc, &
                    max_nInc,max_nInc,iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),iWork(ipnInc),iWork(ipmDec),iWork(ipnDec), &
                    Work(ipH2),Work(ipS2),Work(ipHess2),Work(ipG2),Work(ipBase2),Work(ipr02),nnsiz,nDimTot,nOsc)
    call SolveSecEq(Work(ipH2),nDimTot,Work(ipU2),Work(ipS2),Work(ipE2))
    write(u6,'(20f10.1)') ((Work(ipE2+i-1)-Work(ipE2))*auTocm,i=2,nOsc)
  else
    call SetUpHmat(energy1,Work(ipr1),iWork(ipipow),Work(ipvar),Work(ipyin1),Work(ipr00),trfName1,max_term,Work(ipC1),Work(ipW1), &
                   det1,Work(ipr01),Work(ipC2),Work(ipW2),det2,Work(ipr02),max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,max_nInc, &
                   iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),Work(ipH1),Work(ipS1), &
                   Work(ipG1),Work(ipG2),Work(ipG0),work(ip_gprm1),work(ip_gprm2),work(ip_gprm0),work(ip_gbis1),work(ip_gbis2), &
                   work(ip_gbis0),Work(ipC),Work(ipW),det0,Mass,Work(ipr00),Work(ipBase),Work(ipr0),Work(ipr1),Work(ipr2),nnsiz, &
                   nterm,nvar,ndata,nosc,ndimtot,numofat)
    call SolveSort(Work(ipH1),Work(ipU1),Work(ipS1),Work(ipE1),Work(ipW1),Work(ipW2),Work(ipW1),Work(ipC1),Work(ipC2),Work(ipC1), &
                   Work(ipr01),Work(ipr02),Work(ipr01),iWork(ipmInc),iWork(ipmDec),iWork(ipmMat),mdim1,mdim2,Work(ipOccNumMat1), &
                   nOsc,nDimTot)
    call SetUpHmat(energy2,Work(ipr2),iWork(ipipow),Work(ipvar),Work(ipyin2),Work(ipr00),trfName2,max_term,Work(ipC1),Work(ipW1), &
                   det1,Work(ipr01),Work(ipC2),Work(ipW2),det2,Work(ipr02),max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,max_nInc, &
                   iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),Work(ipH2),Work(ipS2), &
                   Work(ipG1),Work(ipG2),Work(ipG0),work(ip_gprm1),work(ip_gprm2),work(ip_gprm0),work(ip_gbis1),work(ip_gbis2), &
                   work(ip_gbis0),Work(ipC),Work(ipW),det0,Mass,Work(ipr00),Work(ipBase),Work(ipr0),Work(ipr1),Work(ipr2),nnsiz, &
                   nterm,nvar,ndata,nosc,ndimtot,numofat)
    call SolveSort(Work(ipH2),Work(ipU2),Work(ipS2),Work(ipE2),Work(ipW1),Work(ipW2),Work(ipW2),Work(ipC1),Work(ipC2),Work(ipC2), &
                   Work(ipr01),Work(ipr02),Work(ipr02),iWork(ipnInc),iWork(ipnDec),iWork(ipnMat),ndim1,ndim2,Work(ipOccNumMat2), &
                   nOsc,nDimTot)
    call GetMem('yin1','Free','Real',ipyin1,ndata)
    call GetMem('yin2','Free','Real',ipyin2,ndata)

  end if
  write(u6,*) 'T0',T0*HarToRcm
  k2 = 1
  l_TermMat_1 = nDimTot-1
  l_TermMat_2 = nDimTot-1

  call GetMem('TermMat','Allo','Real',ipTermMat,(l_TermMat_1+1)*(l_TermMat_2+1))

  do jOrd=0,nDimTot-1
    k1 = 1
    do iOrd=0,nDimTot-1
      Work(ipTermMat+iOrd+(l_TermMat_1+1)*jOrd) = T0+(Work(ipE2+k2-1)-Work(ipE1+k1-1))
      k1 = k1+1
    end do
    k2 = k2+1
  end do
  call GetMem('H1','Free','Real',ipH1,nDimTot*nDimTot)
  call GetMem('H2','Free','Real',ipH2,nDimTot*nDimTot)
  call GetMem('S1','Free','Real',ipS1,nDimTot*nDimTot)
  call GetMem('S2','Free','Real',ipS2,nDimTot*nDimTot)
  !call GetMem('G0','Free','Real',ipG0,nOsc*nOsc)
end if  ! END: Vibrational Hessian (variational)
call getmem('Gprm0','Free','Real',ip_Gprm0,ngdim**3)
call getmem('Gbis0','Free','Real',ip_Gbis0,ngdim**4)

call GetMem('Hess1','Free','Real',ipHess1,nOsc*nOsc)
call GetMem('Hess2','Free','Real',ipHess2,nOsc*nOsc)
call GetMem('AtCoord','Free','Real',ipAtCoord,3*NumOfAt)
call GetMem('G0','Free','Real',ipG0,nOsc*nOsc)
call GetMem('G1','Free','Real',ipG1,nOsc*nOsc)
call getmem('Gprm1','Free','Real',ip_Gprm1,ngdim**3)
call getmem('Gbis1','Free','Real',ip_Gbis1,ngdim**4)
call GetMem('G2','Free','Real',ipG2,nOsc*nOsc)
call getmem('Gprm2','Free','Real',ip_Gprm2,ngdim**3)
call getmem('Gbis2','Free','Real',ip_Gbis2,ngdim**4)

!----------------------------------------------------------------------!

! Intensity calculations

!----------------------------------------------------------------------!

!PAM: If max_dip=0, there are no transition dipole gradients, but the
!     array TranDipGradInt must be formally all. anyway:
if (max_dip == 0) then
  call GetMem('TranDipGradInt','Allo','Real',ipTranDipGradInt,3*nOsc)

end if

if (Matel) then
  l_IntensityMat_1 = ndimtot-1
  l_IntensityMat_2 = ndimtot-1
  call GetMem('IntensityMat','Allo','Real',ipIntensityMat,(l_IntensityMat_1+1)*(l_IntensityMat_2+1))
  !call GetMem('IntensityMat','Allo','Real',ipIntensityMat,ndimtot*ndimtot)
  if (Forcefield) then
    call dcopy_(nOsc*nOsc,[Zero],0,Work(ipBase2),1)
    !Base2 = Zero
    do i=1,nOsc
      Work(ipBase2+i-1+nOsc*(i-1)) = One
    end do
    call Intensity2(Work(ipIntensityMat),Work(ipTermMat),T0,max_term,Work(ipU1),Work(ipU2),Work(ipE1),Work(ipE2),Work(ipC1), &
                    Work(ipW1),det1,Work(ipr01),Work(ipC2),Work(ipW2),det2,Work(ipr02),Work(ipC),Work(ipW),det0,Work(ipr00),m_max, &
                    n_max,max_dip,iWork(ipm_plot),iWork(ipn_plot),TranDip,Work(ipTranDipGradInt),Work(ipharmfreq1), &
                    Work(ipx_anharm1),Work(ipharmfreq2),Work(ipx_anharm2),Work(ipr0),Work(ipr1),Work(ipr2),Work(ipBase2), &
                    l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nOsc,nDimTot,l_n_plot,l_m_plot)
    call GetMem('TranDipGradInt','Free','real',ipTranDipGradInt,3*nOsc)
  else
    call Intensity(Work(ipIntensityMat),Work(ipTermMat),T0,max_term,iWork(ipipow),Work(ipvar),Work(ipt_dipin1),Work(ipt_dipin2), &
                   Work(ipt_dipin3),trfName1,Work(ipU1),Work(ipU2),Work(ipE1),Work(ipE2),Work(ipC1),Work(ipW1),det1,Work(ipr01), &
                   Work(ipC2),Work(ipW2),det2,Work(ipr02),work(ipC),Work(ipW),det0,Work(ipr00),m_max,n_max,max_dip, &
                   iWork(ipm_plot),iWork(ipn_plot),Work(ipharmfreq1),Work(ipharmfreq2),Work(ipr0),Work(ipr1),Work(ipr2), &
                   Work(ipBase),l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nOsc,nDimTot,nPolyTerm,ndata,nvar, &
                   MaxNumAt,l_n_plot,l_m_plot)
    call GetMem('t_dipin1','Free','Real',ipt_dipin1,ndata)
    call GetMem('t_dipin2','Free','Real',ipt_dipin2,ndata)
    call GetMem('t_dipin3','Free','Real',ipt_dipin3,ndata)

    call GetMem('var','Free','Real',ipvar,ndata*nvar)
    call GetMem('ipow','Free','Inte',ipipow,nPolyTerm*nvar)

  end if
else
  l_IntensityMat_1 = max_mOrd
  l_IntensityMat_2 = max_nOrd

  call GetMem('IntensityMat','Allo','Real',ipIntensityMat,(l_IntensityMat_1+1)*(l_IntensityMat_2+1))
  call IntForceField(Work(ipIntensityMat),Work(ipTermMat),T0,max_term,FC00,Work(ipC1),Work(ipW1),det1,Work(ipr01),Work(ipC2), &
                     Work(ipW2),det2,Work(ipr02),Work(ipC),Work(ipW),det0,Work(ipr00),m_max,n_max,max_dip,Trandip, &
                     Work(ipTranDipGradInt),Work(ipharmfreq1),Work(ipx_anharm1),Work(ipharmfreq2),Work(ipx_anharm2),iWork(ipmMat), &
                     iWork(ipmInc),iWork(ipmDec),iWork(ipnMat),iWork(ipnInc),iWork(ipnDec),OscStr,nsize,max_mOrd,max_nOrd,nDimTot, &
                     nOsc)
  call GetMem('TranDipGradInt','Free','Real',ipTranDipGradInt,3*nOsc)
end if

!write results to log.
write(u6,*) ' Write intensity data to log file.'
call XFlush(u6)
call WriteInt(Work(ipIntensityMat),Work(ipTermMat),iWork(ipmMat),iWork(ipnMat),Work(ipOccNumMat1),Work(ipOccNumMat2),MatEl, &
              ForceField,Work(ipE1),Work(ipE2),T0,Work(ipharmfreq1),Work(ipharmfreq2),Work(ipx_anharm1),Work(ipx_anharm2), &
              l_IntensityMat_1,l_IntensityMat_2,l_TermMat_1,l_TermMat_2,nDimTot,nOsc)

call GetMem('C','Free','Real',ipC,nOsc*nOsc)
call GetMem('W','Free','Real',ipW,nOsc*nOsc)
call GetMem('W1','Free','Real',ipW1,nOsc*nOsc)
call GetMem('W2','Free','Real',ipW2,nOsc*nOsc)

call GetMem('C1','Free','Real',ipC1,nOsc*nOsc)
call GetMem('C2','Free','Real',ipC2,nOsc*nOsc)

call GetMem('r00','Free','Real',ipr00,nOsc)
call GetMem('r01','Free','Real',ipr01,nOsc)
call GetMem('r02','Free','Real',ipr02,nOsc)

call GetMem('nInc','Free','Inte',ipnInc,(nTabDim+1)*nOsc)
call GetMem('nDec','Free','Inte',ipnDec,(nTabDim+1)*nOsc)
call GetMem('mInc','Free','Inte',ipmInc,(mTabDim+1)*nOsc)
call GetMem('mDec','Free','Inte',ipmDec,(mTabDim+1)*nOsc)

call GetMem('harmfreq1','Free','Real',ipharmfreq1,nOsc)
call GetMem('harmfreq2','Free','Real',ipharmfreq2,nOsc)

call GetMem('x_anharm1','Free','Real',ipx_anharm1,nOsc*nOsc)
call GetMem('x_anharm2','Free','Real',ipx_anharm2,nOsc*nOsc)
call GetMem('TermMat','Free','Real',ipTermMat,(l_TermMat_1+1)*(l_TermMat_2+1))
call GetMem('IntensityMat','Free','Real',ipIntensityMat,(l_IntensityMat_1+1)*(l_IntensityMat_2+1))

call GetMem('U1','Free','Real',ipU1,l_l*l_l)
call GetMem('E1','Free','Real',ipE1,l_l)
call GetMem('U2','Free','Real',ipU2,l_l*l_l)
call GetMem('E2','Free','Real',ipE2,l_l)
call GetMem('mMat','Free','Inte',ipmMat,(mTabDim+1)*nOsc)
call GetMem('nMat','Free','Inte',ipnMat,(nTabDim+1)*nOsc)

call GetMem('OccNumMat1','Free','Real',ipOccNumMat1,l_l*3)
call GetMem('OccNumMat2','Free','Real',ipOccNumMat2,l_l*3)

call GetMem('m_plot','Free','Inte',ipm_plot,l_m_plot)
call GetMem('n_plot','Free','Inte',ipn_plot,l_n_plot)

call GetMem('NormModes','Free','Inte',ipNormModes,l_NormModes)
call GetMem('Base','Free','Real',ipBase,nOscOld*nOsc)
call GetMem('BaseInv','Free','Real',ipBaseInv,nOsc*nOscOld)

! I don't understand why is it here??
call GetMem('TranDipGrad','Free','Real',ipTranDipGrad,3*NumOfAt)
call GetMem('Base2','Free','Real',ipBase2,nOsc*nOsc)
call GetMem('Lambda','Free','Real',ipLambda,nOsc)
call GetMem('r1','Free','Real',ipr1,nOsc)
call GetMem('r2','Free','Real',ipr2,nOsc)
call GetMem('r0','Free','Real',ipr0,nOsc)

999 continue

iReturn = 0

end subroutine Mula
