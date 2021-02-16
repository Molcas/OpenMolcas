************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Niclas Forsberg                                        *
*               2009, Giovanni Ghigo                                   *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine Mula(ireturn)
C!
C!      Purpose:
C!      Calculate vibronic intensities.
C!
C!
C!      Written by:
C!      Niclas Forsberg,
C!      Dept. of Theoretical Chemistry, Chemical Centre, Lund University.
C!
C!      Modified by:
C!      Giovanni Ghigo,
C!      Dip. Chimica Generale e Chimica Organica,
C!      Universita' di Torino, ITALY.
C!      Inter-System Crossing rate constant.
C!
C!-----------------------------------------------------------------------!
C!
C!
      Implicit Real*8 ( a-h,o-z )
#include "inout.fh"
      Integer InterVec(15*MaxNumAt )
      Character*80 trfName1(MaxNumAt) ,trfName2(MaxNumAt)
      Real*8 Mass( MaxNumAt )
      Character*4 AtomLbl(MaxNumAt), filnam
      Integer Bond(2*MaxNumAt)
      Character*80 Title
      Real*8       stand_dev,max_err
      Logical      find_minimum
      Logical      lExpan
      Logical      use_weight
      Logical      MatEl,ForceField, lISC,lOldCode,lInCore
      Logical      Cartesian
      Logical      exist
      Integer  nvTabDim
      Integer iCode, iMaxYes, lNMAT0, lNMAT, lNINC, lNDEC
#include "Constants_mula.fh"
#include "inputdata.fh"
#include "dims.fh"
#include "indims.fh"
#include "WrkSpc.fh"
#include "io_mula.fh"
#include "warnings.fh"
C!
C!
C!----      Initialize.
      iReturn=20
      close(5)
      call molcas_open(inpunit,'stdin')
      lPotCoef=ip_Dummy
      iCode = 00       ! Default: NewCode, InCore
      lOldCode=.False. ! X0  CGGn.4
      lInCore=.True.   ! 0X  CGGn.4
      lNMAT0= 92       !     CGGn.4
      lNMAT = 93       !     CGGn.4
      lNINC = 94       !     CGGn.4
      lNDEC = 95       !     CGGn.4
      iPrint = iPrintLevel(-1)
      iPrint_Save = iPrint
c      open(unit=inpUnit,file='stdin')
C!
C!----      Read input file and Write header to log file.
      Call ReadInp(Title,AtomLbl,Mass,InterVec,Bond,nBond,nOsc,
     &             NumOfAt,trfName1,trfName2,m_max,n_max,max_dip,
     &  max_term,MatEl,ForceField,Cartesian,lExpan,lISC,iCode,dMinWind)
      If (lISC) then
        If (iCode.EQ.01 .or. iCode.EQ. 11) lOldCode=.True.
        If (iCode.EQ.10 .or. iCode.EQ. 11) lInCore=.False.
        If (energy1.LT.energy2) then
           Write(6,*) ' ******************** ERROR *******************'
           write(6,*) ' Energy of State #2 must be lower than State #1'
           Write(6,*) ' **********i***********************************'
           Write(6,*)
           Call Quit_OnUserError()
        EndIf
        MatEl = .False.
        m_max = 0
        If (iPrint.GE.1) then
          Write(6,*)
          Write(6,*) ' ---------------------------------------------'
          Write(6,*) ' InterSystem Crossing rate constant evaluation'
          Write(6,*) '  - simple harmonic approximation is used.    '
          Write(6,*) '  - m_max set to 0.                           '
          Write(6,*) '  - no plot file.                             '
          If (lOldCode) then
            Write(6,*) '  - Full index matrices evaluation.           '
          else
            Write(6,*) '  - Reduced matrices evaluation.              '
          EndIf
          If (lInCore) then
            Write(6,*) '  - In-core (memory) algorithm.               '
          else
            Write(6,*) '  - Out-of-core (disk) algorithm.             '
          EndIf
          Write(6,*) ' ---------------------------------------------'
          Write(6,*)
        EndIf
        iPrint = 0
        Do i = 0, maxMax_n
          nIndex(1,i) = 0
          nIndex(2,i) = 0
          nIndex(3,i) = 0
        EndDo
      EndIf
      If (iPrint.GE.1 .and. Title(1:1).NE.' ') Call WriteHeader(Title)

C!
      close(unit=inpUnit)
      Call GetMem('r01','Allo','Real',ipr01,nOsc)
      Call GetMem('r02','Allo','Real',ipr02,nOsc)
      Call GetMem('r1','Allo','Real',ipr1,nOsc)
      Call GetMem('r2','Allo','Real',ipr2,nOsc)
      Call GetMem('r0','Allo','Real',ipr0,nOsc)

      Call GetMem('grad1','Allo','Real',ipgrad1,nOsc)
      Call GetMem('grad2','Allo','Real',ipgrad2,nOsc)

      Call GetMem('G1','Allo','Real',ipG1,nOsc*nOsc)
      Call GetMem('G2','Allo','Real',ipG2,nOsc*nOsc)
      Call GetMem('G0','Allo','Real',ipG0,nOsc*nOsc)

      Call GetMem('eigenVec1','Allo','Real',ipeigenVec1,nOsc*nOsc)
      Call GetMem('eigenVec2','Allo','Real',ipeigenVec2,nOsc*nOsc)

      Call GetMem('harmfreq1','Allo','Real',ipharmfreq1,nOsc)
      Call GetMem('harmfreq2','Allo','Real',ipharmfreq2,nOsc)

      Call GetMem('anharmfreq1','Allo','Real',ipanharmfreq1,nOsc)
      Call GetMem('anharmfreq2','Allo','Real',ipanharmfreq2,nOsc)

      Call GetMem('qMat','Allo','Real',ipqMat,3*NumOfAt*nOsc)

      Call GetMem('PED','Allo','Real',ipPED,nOsc**3)
      Call GetMem('x_anharm1','Allo','Real',ipx_anharm1,nOsc*nOsc)
      Call GetMem('x_anharm2','Allo','Real',ipx_anharm2,nOsc*nOsc)


      If(max_term.ge.0.and.max_term.le.2) then
      ngdim=1
      Else If(max_term.ge.3.and.max_term.le.4) then
      ngdim=nosc
      Else
      Write(6,*) 'MULA error: max_term=',max_term
      Write(6,*) 'Allowed: 1,2,3, or 4.'
      Call Quit_OnUserError()
      Endif
      Call getmem('D3_1','Allo','Real',ip_D3_1,ngdim**3)
      Call getmem('D3_2','Allo','Real',ip_D3_2,ngdim**3)
      Call getmem('D4_1','Allo','Real',ip_D4_1,ngdim**4)
      Call getmem('D4_2','Allo','Real',ip_D4_2,ngdim**4)
      Call getmem('Gprm1','Allo','Real',ip_Gprm1,ngdim**3)
      Call getmem('Gprm2','Allo','Real',ip_Gprm2,ngdim**3)
      Call getmem('Gprm0','Allo','Real',ip_Gprm0,ngdim**3)
      Call getmem('Gbis1','Allo','Real',ip_Gbis1,ngdim**4)
      Call getmem('Gbis2','Allo','Real',ip_Gbis2,ngdim**4)
      Call getmem('Gbis0','Allo','Real',ip_Gbis0,ngdim**4)
C!
C!-----------------------------------------------------------------------!
C!
C!      Either fit polynomial to energies and find minimum or use
C!      forcefield and equilibrium geometry given in input.
C!
C!----------------------------      First State ------------------------------!
C!
      Call GetMem('AtCoord','Allo','Real',ipAtCoord,3*NumOfAt)

      l_a=NumOfAt
      l_Hess_1=l_Hess1
      l_Hess_2=l_Hess1

      If ( ForceField ) Then
      call dcopy_(3*NumOfAt,Work(ipAtCoord1),1,
     &      Work(ipAtCoord),1)
      Else
      find_minimum = .true.
      use_weight = .false.
      Call GetMem('PotCoef','Allo','Real',lPotCoef,nPolyTerm)
      Call PotFit(nPolyTerm,nvar,ndata,iWork(ipipow),
     &      Work(ipvar),Work(ipyin1),
     &      Work(lPotCoef),Work(ipr01),nOsc,energy1,Work(ipgrad1),
     &      Work(ipHess1),
     &              work(ip_D3_1),work(ip_D4_1),
     &              trfName1,stand_dev,max_err,find_minimum,max_term,
     &              use_weight,l_Hess_1,l_Hess_2,ngdim)
      Call Int_To_Cart1(InterVec,Work(ipr01),Work(ipAtCoord),
     &      l_a,nOsc)
      End If
      Call GetMem('AtCoord1','Free','Real',ipAtCoord1,3*NumOfAt)
C!
C!----      Determine vibrational modes and their frequencies.
C!D      Write(6,*)' MULA calling VIBFREQ.'
      Call VibFreq(Work(ipAtCoord),Work(ipr01),InterVec,
     &      Mass,Work(ipHess1),Work(ipG1),
     &      work(ip_gprm1),work(ip_gbis1),
     &      Work(ipharmfreq1),Work(ipeigenVec1),Work(ipqMat),
     &      Work(ipPED),work(ip_D3_1),work(ip_D4_1),
     &      Work(ipx_anharm1),Work(ipanharmfreq1),max_term,Cartesian,
     &      nOsc,NumOfAt)
C!

cvv use ivv to prevent overoptimization of the code..
         ivv=lPotCoef
c      Call WriteLog(Work(lPotCoef),AtomLbl,Work(ipAtCoord),
      If (iPrint.GE.1) Call WriteLog(Work(ivv),AtomLbl,Work(ipAtCoord),
     &      Mass,InterVec,
     &       stand_dev,max_err, energy1,Work(ipHess1),
     &     Work(ipG1),Work(ipeigenVec1),
     &                 Work(ipharmfreq1),Work(ipqMat),Bond,
     &       nBond,Work(ipr01),work(ip_D3_1),
     &                 work(ip_D4_1),Work(ipPED),
     &          Work(ipx_anharm1),Work(ipanharmfreq1),
     &                 max_term,1,ForceField,NumOfAt,nOsc)
         lPotCoef=ivv

C!
C!----------------------------      Second State -----------------------------!
C!
      l_Hess_1=l_Hess1
      l_Hess_2=l_Hess1

      If ( ForceField ) Then
      call dcopy_(3*NumOfAt,Work(ipAtCoord2),1,Work(ipAtCoord),1)
      Else
      find_minimum = .true.
      use_weight = .false.
      Call PotFit(nPolyTerm,nvar,ndata,iWork(ipipow),
     &      Work(ipvar),Work(ipyin2),
     &              Work(lPotCoef),Work(ipr02),nOsc,energy2,
     &     Work(ipgrad2),Work(ipHess2),
     &              work(ip_D3_2),work(ip_D4_2),
     &              trfName2,stand_dev,max_err,find_minimum,max_term,
     &              use_weight,l_Hess_1,l_Hess_2,ngdim)
      Call Int_To_Cart1(InterVec,Work(ipr02),Work(ipAtCoord),
     &      l_a,nOsc)
      End If
      Call GetMem('AtCoord2','Free','Real',ipAtCoord2,3*NumOfAt)
C!
C!----      Determine vibrational modes and their frequencies.
      Call VibFreq(Work(ipAtCoord),Work(ipr02),InterVec,
     &      Mass,Work(ipHess2),Work(ipG2),
     &              work(ip_gprm2),work(ip_gbis2),
     &              Work(ipharmfreq2),Work(ipeigenVec2),
     &     Work(ipqMat),Work(ipPED),work(ip_D3_2),
     &              work(ip_D4_2),Work(ipx_anharm2),
     &              Work(ipanharmfreq2),max_term,Cartesian,
     &              nOsc,NumOfAt)
      If (iPrint.GE.1)
     &   Call WriteLog(Work(lPotCoef),AtomLbl,Work(ipAtCoord),
     &      Mass,InterVec,
     &                 stand_dev,max_err, energy2,
     &   Work(ipHess2),Work(ipG2),Work(ipeigenVec2),
     &                 Work(ipharmfreq2),Work(ipqMat),Bond,nBond,
     &    Work(ipr02),work(ip_D3_2),
     &                 work(ip_D4_2),Work(ipPED),Work(ipx_anharm2),
     &    Work(ipanharmfreq2),
     &                 max_term,2,ForceField,NumOfAt,nOsc)
C!
      If ( .not.Forcefield )
     &   Call GetMem('PotCoef','Free','Real',lPotCoef,nPolyTerm)
      Call GetMem('PED','Free','Real',ipPED,nOsc**3)
      Call GetMem('qMat','Free','Real',ipqMat,3*NumOfAt*nOsc)

      Call GetMem('grad1','Free','Real',ipgrad1,nOsc)
      Call GetMem('grad2','Free','Real',ipgrad2,nOsc)

      Call GetMem('anharmfreq1','Free','Real',ipanharmfreq1,nOsc)
      Call GetMem('anharmfreq2','Free','Real',ipanharmfreq2,nOsc)

      Call getmem('D3_1','Free','Real',ip_D3_1,ngdim**3)
      Call getmem('D3_2','Free','Real',ip_D3_2,ngdim**3)
      Call getmem('D4_1','Free','Real',ip_D4_1,ngdim**4)
      Call getmem('D4_2','Free','Real',ip_D4_2,ngdim**4)
C!
C!-----------------------------------------------------------------------!
C!
C!      Geometry of intermediate oscillator
C!
C!-----------------------------------------------------------------------!
C!
      T0 = energy2-energy1
      Call GetMem('C1','Allo','Real',ipC1,nOsc*nOsc)
      Call GetMem('C2','Allo','Real',ipC2,nOsc*nOsc)
      Call GetMem('C','Allo','Real',ipC,nOsc*nOsc)
      Call GetMem('W','Allo','Real',ipW,nOsc*nOsc)
      Call GetMem('W1','Allo','Real',ipW1,nOsc*nOsc)
      Call GetMem('W2','Allo','Real',ipW2,nOsc*nOsc)

      Call GetMem('temp','Allo','Real',iptemp,nOsc*nOsc)
C!
C!----      Calculate W matrices.
      Do jOsc = 1,nOsc
      const1 = 1.0d0/sqrt(Work(ipharmfreq1+jOsc-1))
      const2 = 1.0d0/sqrt(Work(ipharmfreq2+jOsc-1))
      do iv=1,nOsc
      Work(ipW1+iv-1+nOsc*(jOsc-1)) = const1*
     &        Work(ipeigenVec1+iv-1+nOsc*(jOsc-1))
      Work(ipW2+iv-1+nOsc*(jOsc-1)) = const2*
     &         work(ipeigenVec2+iv-1+nOsc*(jOsc-1))
      enddo
      End Do
C!
C!----      Calculate C = W^(-1).
      call dcopy_(nOsc**2,0.0d0,0,Work(ipC1),1)
      call dcopy_(nOsc,1.0d0,0,Work(ipC1),nOsc+1)
      call dcopy_(nOsc*nOsc,Work(ipW1),1,Work(iptemp),1)
      Call Dool(Work(iptemp),nOsc,nOsc,Work(ipC1),
     &        nOsc,nOsc,det)
      det1 = abs(1.0d0/det)
      call dcopy_(nOsc**2,0.0d0,0,Work(ipC2),1)
      call dcopy_(nOsc,1.0d0,0,Work(ipC2),nOsc+1)
      call dcopy_(nOsc*nOsc,Work(ipW2),1,Work(iptemp),1)
      Call Dool(Work(iptemp),nOsc,nOsc,Work(ipC2),
     &        nOsc,nOsc,det)
      det2 = abs(1.0d0/det)
      Call GetMem('temp','Free','Real',iptemp,nOsc*nOsc)
C!
C!----      Calculate the expansion point geometry and save the full
C!      geometries for later use
      If (iPrint.GE.1) Call ExpPointHeader
      call GetMem('r00','Allo','Real',ipr00,nOsc)
      Call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)

      Call Calc_r00(Work(ipC1),Work(ipC2),Work(ipW1),
     &      Work(ipW2),work(ipC),
     &      Work(ipW),Work(ipalpha1),Work(ipalpha2),
     &      Work(ipr00),Work(ipr01),
     &      Work(ipr02),det0,det1,det2,FC00,nOsc)
      call dcopy_(nosc,Work(ipr00),1,Work(ipr0),1)
      call dcopy_(nosc,Work(ipr01),1,Work(ipr1),1)
      call dcopy_(nosc,Work(ipr02),1,Work(ipr2),1)
      Call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
      l_a=NumOfAt
      Call Int_To_Cart1(InterVec,Work(ipr00),
     &      Work(ipAtCoord),l_a,nOsc)
      If (iPrint.GE.1)
     &  Call WriteCartCoord(AtomLbl,Work(ipAtCoord),Mass,NumOfAt)
      call GetMem('r00','Free','Real',ipr00,nOsc)
C!
C!----      If we only wanted the expansion point geometry, its time to quit now.
      If ( lExpan ) Then
      Write(6,*) 'Bye bye'
      Write(6,*)
      Call Quit(_RC_ALL_IS_WELL_)
      End If
C!
      If (iPrint.GE.1) Call IntCalcHeader
C!
C!----      Calculate term values.
      GE1 = 0.0d0
      GE2 = GE1+T0
      Do iOsc = 1,nOsc
      jOsc = 1
      exist = .false.
      Do While (( .not.exist ).and.( jOsc.le.l_NormModes ))
      exist = ( iOsc.eq.iWork(ipNormModes+jOsc-1) )
      jOsc = jOsc+1
      End Do
      If ( .not.exist ) Then
      GE1 = GE1+0.5d0*Work(ipharmfreq1+iOsc-1)
      GE2 = GE2+0.5d0*Work(ipharmfreq2+iOsc-1)
      End If
      End Do
      T0 = GE2-GE1
C!
C!----      Pick out the chosen modes.
      nOscOld = nOsc
      nOsc = l_NormModes
      Call GetMem('Lambda','Allo','Real',ipLambda,nOsc)
      Call GetMem('Base','Allo','Real',ipBase,nOscOld*nOsc)
      Call GetMem('BaseInv','Allo','Real',
     &      ipBaseInv,nOsc*nOscOld)
      Do jOsc = 1,nOsc
      Do iOsc = 1,nOscOld
      Work(ipBase+iOsc-1+nOscOld*(jOsc-1))    =
     &      Work(ipW2+iOsc-1+nOsc*(iWork(ipNormModes+jOsc-1)-1))
      Work(ipBaseInv+jOsc-1+nOsc*(iOsc-1)) =
     &      Work(ipC2+iWork(ipNormModes+jOsc-1)-1+nOsc*(iOsc-1))
      End Do
      End Do
C!
C!----      First state.
C!      Subroutine SolveRedSec(Hess,Gmat,freq,C,W,det)
      Call GetMem('temp','Allo','Real',iptemp,nOscOld*nOsc)
      Call DGEMM_('N','N',
     &            nOscOld,nOsc,nOscOld,
     &            1.0d0,Work(ipHess1),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(iptemp),nOscOld)
      Call GetMem('Hess1','Free','Real',ipHess1,l_Hess1*l_Hess1)

      Call GetMem('Hess1','Allo','Real',ipHess1,nOsc*nOsc)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOscOld,
     &            1.0d0,Work(ipBase),nOscOld,
     &            Work(iptemp),nOscOld,
     &            0.0d0,Work(ipHess1),nOsc)
      Call GetMem('temp2','Allo','Real',iptemp2,nOscOld*nOscOld)
      call dcopy_(nOscOld**2,0.0d0,0,Work(iptemp2),1)
      call dcopy_(nOscOld,1.0d0,0,Work(iptemp2),nOscOld+1)
      Call Dool(Work(ipG1),nOsc,nOsc,work(iptemp2),
     &      nOscOld,nOscOld,det)
      Call DGEMM_('N','N',
     &            nOscOld,nOsc,nOscOld,
     &            1.0d0,Work(iptemp2),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(iptemp),nOscOld)
      Call GetMem('temp2','Free','Real',iptemp2,nOscOld*nOscOld)
      Call GetMem('temp2','Allo','Real',iptemp2,nOsc*nOsc)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOscOld,
     &            1.0d0,Work(ipBase),nOscOld,
     &            Work(iptemp),nOscOld,
     &            0.0d0,Work(iptemp2),nOsc)
      call dcopy_(nOsc**2,0.0d0,0,Work(ipG1),1)
      call dcopy_(nOsc,1.0d0,0,Work(ipG1),nOsc+1)
      Call Dool(Work(iptemp2),nOscOld,nOscOld,
     &      Work(ipG1),nOsc,nOsc,det)
      Call GetMem('temp2','Free','Real',iptemp2,nOsc*nOsc)
      Call GetMem('temp','Free','Real',iptemp,nOscOld*nOsc)
      Call GetMem('temp','Allo','Real',iptemp,nOsc*nOsc)
      Call SolveSecEq(Work(ipHess1),nOsc,Work(iptemp),
     &      Work(ipG1),Work(ipLambda))
      do iv=1,nOsc
      Work(ipharmfreq1+iv-1) = sqrt(abs(Work(ipLambda+iv-1)))
      enddo
      If (iPrint.GE.1)
     &  Call WriteFreq(Work(ipharmfreq1),iWork(ipNormModes),
     &      l_NormModes,
     &      'Frequencies of reduced problem, state 1')
      Do jOsc = 1,nOsc
      const1 = 1.0/sqrt(Work(ipharmfreq1+jOsc-1))
      do iv=1,nOsc
      Work(ipW1+iv-1+nOsc*(jOsc-1)) = const1*
     &      Work(iptemp+iv-1+nOsc*(jOsc-1))
      enddo
      End Do
      call dcopy_(nOsc**2,0.0d0,0,Work(ipC1),1)
      call dcopy_(nOsc,1.0d0,0,Work(ipC1),nOsc+1)
      call dcopy_(nOsc*nOsc,Work(ipW1),1,Work(iptemp),1)
      Call Dool(Work(iptemp),nOscOld,nOsc,Work(ipC1),
     &        nOsc,nOsc,det)
      det1 = abs(1.0d0/det)
      Call GetMem('temp','Free','Real',iptemp,nOsc*nOsc)
C!
C!----      Second state.
      Call GetMem('temp','Allo','Real',iptemp,nOscOld*nOsc)
      Call DGEMM_('N','N',
     &            nOscOld,nOsc,nOscOld,
     &            1.0d0,Work(ipHess2),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(iptemp),nOscOld)

      Call GetMem('Hess2','Free','Real',ipHess2,l_Hess1*l_Hess1)
      Call GetMem('Hess2','Allo','Real',ipHess2,nOsc*nOsc)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOscOld,
     &            1.0d0,Work(ipBase),nOscOld,
     &            Work(iptemp),nOscOld,
     &            0.0d0,Work(ipHess2),nOsc)
      Call GetMem('temp2','Allo','Real',iptemp2,nOscOld*nOscOld)
      call dcopy_(nOscOld**2,0.0d0,0,Work(iptemp2),1)
      call dcopy_(nOscOld,1.0d0,0,Work(iptemp2),nOscOld+1)
      Call Dool(Work(ipG2),nOsc,nOsc,Work(iptemp2),
     &      nOscOld,nOscOld,det)
      Call DGEMM_('N','N',
     &            nOscOld,nOsc,nOscOld,
     &            1.0d0,Work(iptemp2),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(iptemp),nOscOld)
      Call GetMem('temp2','Free','Real',iptemp2,nOscOld*nOscOld)
      Call GetMem('temp2','Allo','Real',iptemp2,nOsc*nOsc)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOscOld,
     &            1.0d0,Work(ipBase),nOscOld,
     &            Work(iptemp),nOscOld,
     &            0.0d0,Work(iptemp2),nOsc)
      call dcopy_(nOsc**2,0.0d0,0,Work(ipG2),1)
      call dcopy_(nOsc,1.0d0,0,Work(ipG2),nOsc+1)
      Call Dool(Work(iptemp2),nOsc,nOsc,Work(ipG2),
     &      nOsc,nOsc,det)
      Call GetMem('temp2','Free','Real',iptemp2,nOsc*nOsc)
      Call GetMem('temp','Free','Real',iptemp,nOscOld*nOsc)
      Call GetMem('temp','Allo','Real',iptemp,nOsc*nOsc)
      Call SolveSecEq(Work(ipHess2),nOsc,Work(iptemp),
     &      Work(ipG2),Work(ipLambda))
      do iv=1,nOsc
      Work(ipharmfreq2+iv-1) = sqrt(abs(Work(ipLambda+iv-1)))
      enddo
      If (iPrint.GE.1)
     &  Call WriteFreq(Work(ipharmfreq2),iWork(ipNormModes),
     &      l_NormModes,
     &      'Frequencies of reduced problem, state 2')
      Do jOsc = 1,nOsc
      const1 = 1.0d0/sqrt(Work(ipharmfreq2+jOsc-1))
      do iv=1,nOsc
      Work(ipW2+iv-1+nOsc*(jOsc-1)) = const1*
     &      Work(iptemp+iv-1+nOsc*(jOsc-1))
      enddo
      End Do
      call dcopy_(nOsc**2,0.0d0,0,Work(ipC2),1)
      call dcopy_(nOsc,1.0d0,0,Work(ipC2),nOsc+1)
      call dcopy_(nOsc*nOsc,Work(ipW2),1,Work(iptemp),1)
      Call Dool(Work(iptemp),nOsc,nOsc,Work(ipC2),
     &        nOsc,nOsc,det)
      det2 = abs(1.0d0/det)
      Call GetMem('temp','Free','Real',iptemp,nOsc*nOsc)
C!
      Call GetMem('rtemp','Allo','Real',iprtemp,nOscOld)

      call dcopy_(nOscOld,Work(ipr01),1,Work(iprtemp),1)
      Call GetMem('r01','Free','Real',ipr01,nOsc)
      Call GetMem('r01','Allo','Real',ipr01,nOsc)

      Call DGEMM_('N','N',
     &            nOsc,1,nOscOld,
     &            1.0d0,Work(ipBaseInv),nOsc,
     &            Work(iprtemp),nOscOld,
     &            0.0d0,Work(ipr01),nOsc)
      call dcopy_(nOscOld,Work(ipr02),1,Work(iprtemp),1)

      Call GetMem('r02','Free','Real',ipr02,nOsc)
      Call GetMem('r02','Allo','Real',ipr02,nOsc)
      Call DGEMM_('N','N',
     &            nOsc,1,nOscOld,
     &            1.0d0,Work(ipBaseInv),nOsc,
     &            Work(iprtemp),nOscOld,
     &            0.0d0,Work(ipr02),nOsc)
      Call GetMem('rtemp','Free','Real',iprtemp,nOscOld)
C!
C!----
      call GetMem('r00','Allo','Real',ipr00,nOsc)
      Call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
      Call Calc_r00(Work(ipC1),Work(ipC2),Work(ipW1),
     &      Work(ipW2),Work(ipC),
     &      Work(ipW),Work(ipalpha1),Work(ipalpha2),Work(ipr00),
     &      Work(ipr01),
     &      Work(ipr02),det0,det1,det2,FC00,nOsc)
      Call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
C!
      Call GetMem('Base2','Allo','Real',ipBase2,nOsc*nOsc)

      call dcopy_(nOsc*nOsc,Work(ipW2),1,Work(ipBase2),1)
c            call dcopy_(nOsc*nOsc,Work(ipC2),1,Work(ipBaseInv2),1)
c            Base2 = W2
c            BaseInv2 = C2
C!
C!-----------------------------------------------------------------------!
C!
C!      Transform transition dipole gradients from cartesian to
C!      internal coordinates
C!
C!-----------------------------------------------------------------------!
C!
      If ( Forcefield.and.(max_dip.gt.0) ) Then
      call GetMem('Smat','Allo','Real',ipSmat,3*NumOfAt*nOscOld)
c            Smat = 0.0d0
      call dcopy_(3*NumOfAt*nOscOld,0.0d0,0,Work(ipSmat),1)
      Call CalcS(Work(ipAtCoord),InterVec,Work(ipSmat),
     &      nOscOld,NumOfAt)
      Call GetMem('Bmat','Allo','Real',ipBmat,3*NumOfAt*nOscOld)

      Do j = 1,nOscOld
      k = 1
      Do i = 1,NumOfAt
      Work(ipBmat+k-1+3*NumOfAt*(j-1)) =
     &      Work(ipSmat+  3*(i-1+NumOfAt*(j-1)))
      Work(ipBmat+k  +3*NumOfAt*(j-1)) =
     &      Work(ipSmat+1+3*(i-1+NumOfAt*(j-1)))
      Work(ipBmat+k+1+3*NumOfAt*(j-1)) =
     &      Work(ipSmat+2+3*(i-1+NumOfAt*(j-1)))
      k = k+3
      End Do
      End Do
      call GetMem('Smat','Free','Real',ipSmat,3*NumOfAt*nOscOld)
      Call GetMem('temp','Allo','Real',iptemp,nOscOld*nOscOld)
      Call GetMem('temp3','Allo','Real',iptemp3,nOscOld*nOscOld)

      Call DGEMM_('T','N',
     &            nOscOld,nOscOld,3*NumOfAt,
     &            1.0d0,Work(ipBmat),3*NumOfAt,
     &            Work(ipBmat),3*NumOfAt,
     &            0.0d0,Work(iptemp),nOscOld)
      call dcopy_(nOscOld**2,0.0d0,0,Work(iptemp3),1)
      call dcopy_(nOscOld,1.0d0,0,Work(iptemp3),nOscOld+1)
      Call Dool(Work(iptemp),nOscOld,nOscOld,Work(iptemp3),
     &      nOscOld,nOscOld,det)
      Call GetMem('temp','Free','Real',iptemp,nOscOld*nOscOld)
      Call GetMem('temp1','Allo','Real',
     &      iptemp1,3*NumOfAt*nOscOld)

      Call DGEMM_('N','N',
     &            3*NumOfAt,nOscOld,nOscOld,
     &            1.0d0,Work(ipBmat),3*NumOfAt,
     &            Work(iptemp3),nOscOld,
     &            0.0d0,Work(iptemp1),3*NumOfAt)
      Call GetMem('Bmat','Allo','Real',ipBmat,3*NumOfAt*nOscOld)
      Call GetMem('temp3','Free','Real',iptemp3,nOscOld*nOscOld)
      Call GetMem('temp3','Allo','Real',iptemp3,3*NumOfAt*nOsc)
      Call DGEMM_('N','N',
     &            3*NumOfAt,nOsc,nOscOld,
     &            1.0d0,Work(iptemp1),3*NumOfAt,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(iptemp3),3*NumOfAt)
      Call GetMem('temp1','Free','Real',
     &      iptemp1,3*NumOfAt*nOscOld)
      Call GetMem('temp1','Allo','Real',iptemp1,3*NumOfAt*nOsc)
      Call DGEMM_('N','N',
     &            3*NumOfAt,nOsc,nOsc,
     &            1.0d0,Work(iptemp3),3*NumOfAt,
     &            Work(ipW),nOsc,
     &            0.0d0,Work(iptemp1),3*NumOfAt)
      Call GetMem('temp3','Free','Real',iptemp3,3*NumOfAt*nOsc)
      call GetMem('TranDipGradInt','Allo','Real',
     &      ipTranDipGradInt,3*nOsc)
      Call DGEMM_('T','N',
     &            nOsc,1,3*NumOfAt,
     &            1.0d0,Work(iptemp1),3*NumOfAt,
     &            work(ipTranDipGrad),3*NumOfAt,
     &            0.0d0,Work(ipTranDipGradInt),nOsc)
      Call DGEMM_('T','N',
     &            nOsc,1,3*NumOfAt,
     &            1.0d0,Work(iptemp1),3*NumOfAt,
     &            work(ipTranDipGrad+1),3*NumOfAt,
     &            0.0d0,Work(ipTranDipGradInt+1),nOsc)
      Call DGEMM_('T','N',
     &            nOsc,1,3*NumOfAt,
     &            1.0d0,Work(iptemp1),3*NumOfAt,
     &            Work(ipTranDipGrad+2),3*NumOfAt,
     &            0.0d0,Work(ipTranDipGradInt+2),nOsc)
      Call GetMem('temp1','Free','Real',iptemp1,3*NumOfAt*nOsc)
C!
      Call WriteDip(Work(ipTranDipGradInt),iWork(ipNormModes),
     &      'Transition dipole gradient',nOsc)
      End If
C!
      Call GetMem('eigenVec1','Free','Real',ipeigenVec1,nOsc*nOsc)
      Call GetMem('eigenVec2','Free','Real',ipeigenVec2,nOsc*nOsc)
C!
C!-----------------------------------------------------------------------!
C!
C!    InterSystem Crossing section.
C!    Author: Giovanni Ghigo
C!            Dip. Chimica Generale e Chimica Organica, Torino (ITALY)
C!            Vers. 1  28 Dec-08 - 13 Jan-09
C!
C!-----------------------------------------------------------------------!
C!
      If (iPrint.GE.4) Write(6,*) ' FC(0-0)=',FC00
      If (lISC) then
C!
        iPrint = iPrint_Save
        If (WriteVibLevels) iPrint = 2
        If (Huge_Print)     iPrint = 3
        If (iPrint.GE.1) Call ISCHeader
C!
        dRho = 0.0d0
        Call GetMem('nMaxQ','Allo','Inte',ipnMaxQ,nOsc)
        Call ISC_Rho(iPrint,nOsc,new_n_max,dRho,energy1,energy2,minQ,
     &    dMinWind,iWork(ipnMaxQ),Work(ipharmfreq1),Work(ipharmfreq2))
        If (new_n_max.LT.n_max) then
          n_max = new_n_max
          Write(6,*) ' n_max reduced to',n_max
        EndIf
        If (iPrint.GE.3) Write(6,*)' Actual n_max               =',n_max
C!
C!      Set up excitation matrices
C!
C!      Initial State m
        Call TabDim2_drv(m_max,nOsc,nvTabDim)
        mTabDim  = nvTabDim-1
        max_mOrd = mTabDim
        Call TabDim2_drv(m_max-1,nOsc,nvTabDim)
        max_mInc = nvTabDim-1
        Call GetMem('mMat','Allo','Inte',ipmMat,(mTabDim+1)*nOsc)
        Call GetMem('mInc','Allo','Inte',ipmInc,(mTabDim+1)*nOsc)
        Call GetMem('mDec','Allo','Inte',ipmDec,(mTabDim+1)*nOsc)
        mdim1=mTabDim
        mdim2=nOsc
        Call MakeTab2(m_max,max_mOrd,max_mInc,mTabDim,iWork(ipmMat),
     &  iWork(ipmInc),iWork(ipmDec),nOsc)
C!
C!      Initialize Final State n
        n_max_orig = n_max
        Call TabDim2_drv(n_max,nOsc,nvTabDim)
        nTabDim  = nvTabDim-1
        max_nOrd = nTabDim
        Call TabDim2_drv(n_max-1,nOsc,nvTabDim)
        max_nInc = nvTabDim-1
C!
C!      Memory estimation and algorithm selection
        Call GetMem('MULA','MAX','INTE',ipLeft,lLeft)
        If (lInCore) then
          If (lOldCode) then
            iMem = INT(7*(nTabDim+1)*nOsc/2)
            If (iMem.ge.lLeft) then
              Write(6,*) ' Too much memory required (',iMem,' words).'
              Write(6,*) ' Switch to reduced matrices evaluation.'
              lOldCode=.False.
            EndIf
          else
            iMem = INT(3*(nTabDim+1)*nOsc/2)
            If (iMem.ge.lLeft) then
              Write(6,*) ' Too much memory required (',iMem,' words).'
              write(6,*) ' Out-of-Core (disk) algorithm will be used.'
              lInCore=.False.
            EndIf
          EndIf
        else
          Continue
        EndIf
C!
        Call Timing(CPTF0,CPE,TIOTF0,TIOE)                        ! CGGn
        IF (lInCore) THEN

          If (lOldCode) then
            If (iPrint.GE.2) then
              Write(6,*)
              Write(6,*) ' Index matrix evaluation.'
              If (iPrint.GE.3)
     &        Write(6,*) ' Memory allocated for all index matrix:',
     &        iMem,' words,  ',8*iMem/1048576,' MB.'
              Call XFlush(6)
            EndIf
            Call GetMem('nMat','Allo','Inte',ipnMat,(nTabDim+1)*nOsc)
            Call GetMem('nInc','Allo','Inte',ipnInc,(nTabDim+1)*nOsc)
            Call GetMem('nDec','Allo','Inte',ipnDec,(nTabDim+1)*nOsc)
            ndim1=nTabDim
            ndim2=nOsc
            nnsiz=ndim1
            Call MakeTab2(n_max,max_nOrd,max_nInc,nTabDim,
     &      iWork(ipnMat),iWork(ipnInc),iWork(ipnDec),nOsc)
          else ! NewCode: nInc & nDec reduced
            iMem = (nTabDim+1)*nOsc
            If (iPrint.GE.2) then
              Write(6,*)
              Write(6,*) ' First index matrix evaluation.'
              If (iPrint.GE.3)
     &        Write(6,*) ' Memory allocated for first index matrix:',
     &        iMem,' words,  ',(8*iMem+131071)/1048576,' MB.'
              Call XFlush(6)
            EndIf
            Call GetMem('nMat','Allo','Inte',ipnMat,(nTabDim+1)*nOsc)
            Call GetMem('Graph2','Allo','Inte',
     &                  ipGraph2,(n_max+1)*(n_max+1)*nOsc)
            Call GetMem('Graph1','Allo','Inte',
     &                  ipGraph1,(n_max+1)*(nOsc+1))
            ndim1=nTabDim
            ndim2=nOsc
            nnsiz=ndim1
            Call ISC_MakeTab2(n_max,max_nOrd,max_nInc,nTabDim,
     &           iWork(ipnMat),iWork(ipGraph1),iWork(ipGraph2),nOsc)
            Call GetMem('Graph1','Free','Inte',
     &                    ipGraph1,(n_max+1)*(nOsc+1))
          EndIf
          Call Timing(CPTF1,CPE,TIOTF1,TIOE)                      ! CGGn
C!
C!      Definition of Vector with levels in the window
C!
          If (iPrint.GE.2) then
            Write(6,*)
            Write(6,*) ' Energy level screening.'
            Call XFlush(6)
          EndIf
          Call GetMem('lVec','Allo','Inte',iplVec,(nTabDim+1))
          Call LogEVec(iPrint,nOsc,max_nOrd,minQ,iWork(ipnMaxQ),
     &                        iWork(ipnMat),iWork(iplVec),nYes)
          If (nYes.LE.0) then
            Write(6,*)
            Write(6,*) ' ************ ERROR *************'
            write(6,*) ' No energy levels in the Window !'
            Write(6,*) ' ********************************'
            Write(6,*)
            Call Quit_OnUserError()
          EndIf
          Call Timing(CPTF2,CPE,TIOTF2,TIOE)                      ! CGGn
C!
C!      Energy levels calculation
C!
          If (iPrint.GE.2) then
            Write(6,*)
            Write(6,*) ' Energy level calculation.'
            If (iPrint.GE.3) then
              Write(6,*) ' Memory allocated for energy matrix:',
     &          (nTabDim+1),' words,  ',
     &       (8*(nTabDim+1)+131071)/1048576,' MB.'
            EndIf
            Call XFlush(6)
          EndIf
          Call GetMem('lTVec','Allo','Inte',iplTVec,(nTabDim+1))
          Call GetMem('VibLevel2','Allo','Real',ipVibLevel2,max_nOrd+1)
          Call ISC_Ene(iPrint,nOsc,max_nOrd,nYes,iWork(ipnMat),nTabDim,
     &    GE1,GE2,Work(ipharmfreq1),Work(ipharmfreq2),
     &    Work(ipx_anharm1),Work(ipx_anharm2),
     &    dMinWind,dRho,Work(ipVibLevel2),iWork(iplVec),iWork(iplTVec))
          Call GetMem('VibLevel2','Free','Real',ipVibLevel2,max_nOrd+1)
          Call GetMem('lTVec','Free','Inte',iplTVec,(nTabDim+1))
          If (nYes.LE.0) then
            Write(6,*)
            Write(6,*) ' ************ ERROR *************'
            write(6,*) ' No energy levels in the Window !'
            Write(6,*) ' ********************************'
            Write(6,*)
            If (iPrint.LT.4) Call Quit_OnUserError()
          EndIf
C!
C!      The State Window
C!
          Call GetMem('VibWind2','Allo','Inte',ipVibWind2,nYes)
          Call MkVibWind2(iPrint,nYes,iMaxYes,max_nOrd,
     &                    iWork(iplVec),iWork(ipVibWind2))
          Call GetMem('lVec','Free','INTE',iplVec,(nTabDim+1))
          Call Timing(CPTF3,CPE,TIOTF3,TIOE)                      ! CGGn
          Call Timing(CPTF4,CPE,TIOTF4,TIOE)                      ! CGGn
C!
C!      The remaining excitation matrices nInc & nDec
C!
          If (lOldCode) then
            iMaxYes = nTabDim
          else
            iMem = 2*(iMaxYes+1)*nOsc
            If (iPrint.GE.2) then
              Write(6,*)
              Write(6,*) ' nInc and nDec matrices generation.'
              If (iPrint.GE.3) Write(6,*)
     &        ' Memory allocated for remaining index matrix:',
     &        iMem,' words,  ',(8*iMem+131071)/1048576,' MB.'
              Call XFlush(6)
            EndIf
            Call GetMem('nInc','Allo','Inte',ipnInc,(iMaxYes+1)*nOsc)
            Call GetMem('nDec','Allo','Inte',ipnDec,(iMaxYes+1)*nOsc)
            Call Mk_nIncDec(n_max,nTabDim,iMaxYes,
     &          iWork(ipnInc),iWork(ipnDec),iWork(ipnMat),
     &          iWork(ipGraph2),nOsc)
            Call GetMem('Graph2','Free','Inte',
     &                  ipGraph2,(n_max+1)*(m_max+1)*nOsc)
          EndIf
          Call Timing(CPTF5,CPE,TIOTF5,TIOE)                      ! CGGn
C!
C!      The ISC rate
C!
          If (iPrint.GE.2) then
            Write(6,*)
            Write(6,*) ' Franck-Condon factors evaluation.'
            Call XFlush(6)
          EndIf
          Call GetMem('FCWind2','Allo','Real',ipFCWind2,nYes)
          Call ISC_Rate(iPrint,nOsc,max_nOrd,iMx_nOrd,iMaxYes,
     &    nYes,dMinWind,              iWork(ipVibWind2),
     &    Work(ipC1),Work(ipC2),Work(ipW1),Work(ipW2),det0,det1,det2,
     &    Work(ipC),Work(ipW),Work(ipr01),Work(ipr02),Work(ipr00),
     &    mTabDim,iWork(ipmMat),nTabDim,iWork(ipnMat),
     &    iWork(ipmInc),iWork(ipmDec),iWork(ipnInc),iWork(ipnDec),
     &    m_max,n_max,max_dip,nnsiz,FC00,Work(ipFCWind2),dRho)
          Call Timing(CPTF6,CPE,TIOTF6,TIOE)                      ! CGGn

          Call GetMem('FCWind2','Free','Real',ipFCWind2,nYes)
          Call GetMem('VibWind2','Free','Inte',ipVibWind2,nYes)

          If (lOldCode) then
            Call GetMem('nDec','Free','Inte',ipnDec,(nTabDim+1)*nOsc)
            Call GetMem('nInc','Free','Inte',ipnInc,(nTabDim+1)*nOsc)
          else
            Call GetMem('nDec','Free','Inte',ipnDec,(nYes+1)*nOsc)
            Call GetMem('nInc','Free','Inte',ipnInc,(nYes+1)*nOsc)
          EndIf
          Call GetMem('nMat','Free','Inte',ipnMat,(nTabDim+1)*nOsc)

        ELSE ! Out-of-core.

CGGt --- Out-of-core ---------------------------------------------------
C!
C!      Open files
C!
          filnam='MAT0'
          Call DaName_mf_wa(lNMAT0,filnam)
C!
          Call GetMem('Graph2','Allo','Inte',
     &                  ipGraph2,(n_max+1)*(n_max+1)*nOsc)
          Call GetMem('Graph1','Allo','Inte',
     &                  ipGraph1,(n_max+1)*(nOsc+1))
          Call ISCD_MakeGraphs(n_max,max_nOrd,max_nInc,      ! Mk_Graphs
     &                  iWork(ipGraph1),iWork(ipGraph2),nOsc)
          Call GetMem('Graph1','Free','Inte',
     &                  ipGraph1,(n_max+1)*(nOsc+1))
          Call GetMem('nTabDim','Allo','Inte',ipnTabDim,nTabDim+1)
          Call GetMem('nMat0','Allo','Inte',ipnMat0,nOsc)
          If (iPrint.GE.2) then
            Write(6,*)
            Write(6,*) ' Index matrix evaluation.'
            If (iPrint.GE.3) then
              Write(6,*) ' Memory allocated for address array:',
     &          (nTabDim+1),' words,  ',
     &       (8*(nTabDim+1)+131071)/1048576,' MB.'
            EndIf
            Call XFlush(6)
          EndIf
          Call ISCD_MakenMat(n_max,nOsc,lNMAT0,nTabDim,iWork(ipGraph2),
     &                  iWork(ipnTabDim),iWork(ipnMat0))
          Call Timing(CPTF1,CPE,TIOTF1,TIOE)                      ! CGGn
C!
C!      Definition of Vector with levels in the window
C!
          If (iPrint.GE.2) then
            Write(6,*)
            Write(6,*) ' Energy level screening.'
            Call XFlush(6)
          EndIf
          Call GetMem('lVec','Allo','Inte',iplVec,(nTabDim+1))
          Call ISCD_LogEVec(iPrint,nOsc,max_nOrd,minQ,nYes,
     &         lNMAT0,nTabDim, iWork(ipnTabDim),iWork(ipnMaxQ),
     &         iWork(ipnMat0),iWork(iplVec))
          If (nYes.LE.0) then
            Write(6,*)
            Write(6,*) ' ************ ERROR *************'
            write(6,*) ' No energy levels in the Window !'
            Write(6,*) ' ********************************'
            Write(6,*)
            Call Quit_OnUserError()
          EndIf
          Call Timing(CPTF2,CPE,TIOTF2,TIOE)                      ! CGGn
C!
C!      Energy levels calculation
C!
          If (iPrint.GE.2) then
            Write(6,*)
            Write(6,*) ' Energy level calculation.'
            If (iPrint.GE.3) then
              Write(6,*) ' Memory allocated for energy matrix:',
     &        (nTabDim+1),' words,  ',
     &        (8*(nTabDim+1)+131071)/1048576,' MB.'
            EndIf
            Call XFlush(6)
          EndIf
          Call GetMem('lTVec','Allo','Inte',iplTVec,(nTabDim+1))
          Call GetMem('VibLevel2','Allo','Real',ipVibLevel2,max_nOrd+1)
          Call ISCD_Ene(iPrint,nOsc,max_nOrd,nYes,lNMAT0,nTabDim,
     &    GE1,GE2,Work(ipharmfreq1),Work(ipharmfreq2),
     &    Work(ipx_anharm1),Work(ipx_anharm2), dMinWind,dRho,
     &    iWork(ipnMat0),iWork(ipnTabDim),iWork(iplVec),iWork(iplTVec),
     &    Work(ipVibLevel2))
          Call GetMem('VibLevel2','Free','Real',ipVibLevel2,max_nOrd+1)
          Call GetMem('lTVec','Free','Inte',iplTVec,(nTabDim+1))
          If (nYes.LE.0) then
            Write(6,*)
            Write(6,*) ' ************ ERROR *************'
            write(6,*) ' No energy levels in the Window !'
            Write(6,*) ' ********************************'
            Write(6,*)
            If (iPrint.LT.4) Call Quit_OnUserError()
          EndIf
C!
C!      The State Window
C!
          Call GetMem('VibWind2','Allo','Inte',ipVibWind2,nYes)
          Call MkVibWind2(iPrint,nYes,iMaxYes,max_nOrd,
     &                    iWork(iplVec),iWork(ipVibWind2))
          Call GetMem('lVec','Free','Inte',iplVec,(nTabDim+1))
          Call Timing(CPTF3,CPE,TIOTF3,TIOE)                      ! CGGn
C!
C!      Memory estimation & nMat transfer
C!
          If (lOldCode) then
            iMaxYes = nTabDim
          EndIf
          Call GetMem('MULA','MAX','INTE',ipLeft,lLeft)
          If (iPrint.GE.3) then
            Write(6,'(A,I10,A,I4,A)')
     &            '  Available memory                          :',
     &            lLeft,' words,  ',(lLeft+131071)/131072,' MB.'
          EndIf
          iUMem = INT(1*(nTabDim+1))  ! Memory for U, 1=> No Hot states
          iRateMem = iUMem + nTabDim+2 + 13*nOsc*nOsc + 5*nOsc
          iRateMem = INT(2.1D0*iRateMem) ! to convert from INTE to REAL
          lLeft = lLeft - iRateMem
          If (iRateMem.LT.0) then ! .or. iRateMem.GT.2048*131072) then
            Write(6,'(A,I10,A,I4,A)')
     &            '  Estimated memory for FC factors evaluation:',
     &            iRateMem,' words,  ',(iRateMem+131071)/131072,' MB.'
            Write(6,'(A,I10,A,I4,A)')
     &            '  Memory left for batches                   :',
     &            lLeft,' words,  ',(lLeft+131071)/131072,' MB.'
            Write(6,*)
            Write(6,*) ' ****************** ERROR ********************'
            Write(6,*) ' Not enough memory for FC factors evaluation !'
            Write(6,*) ' *********************************************'
            Write(6,*)
            If (iPrint.GE.3) Call GetMem('MULA','LIST','INTE',iDum,iDum)
            Call Quit_OnUserError()
          EndIf
          lMBatch = INT(lLeft/3)
          lMBatch = Min(lMBatch,lMaxMBatch*131072) ! Max 512 MB
          lMBatch = Max(lMBatch,lMinMBatch*131072) ! Min  64 MB
CGGt          lMBatch = 131072*8                       ! CGGt Test only !!!!
          lBatch = INT(lMBatch/nOsc)
          If (lBatch.GT.(iMaxYes+1)) lBatch=(iMaxYes+1)
          lMBatch = lBatch*nOsc
          nBatch = INT( (iMaxYes+1) / lBatch )
          leftBatch = (iMaxYes+1) - nBatch * lBatch
          If (nBatch+1.GT.maxMax_n) then
            Write(6,*)
            Write(6,*) ' ******************** ERROR *******************'
            write(6,*) ' Not enough number of batches for Out-of-core !'
            write(6,*) ' Increase maxMax_n in src/mula/io_mula.fh      '
            Write(6,*) ' **********************************************'
            Write(6,*)
            Call Quit_OnUserError()
          EndIf
          If (iPrint.GE.2) then
            Write(6,*)
            If (iPrint.GE.3) then
              Write(6,'(A,I10,A,I4,A)')
     &              '  Estimated memory for FC factors evaluation:',
     &              iRateMem,' words,  ',(iRateMem+131071)/131072,' MB.'
              Write(6,'(A,I10,A,I4,A)')
     &              '  Memory left for batches                   :',
     &              lLeft,   ' words,  ',(lLeft+131071)/131072,' MB.'
              Write(6,'(A,I10,A,I4,A)')
     &              '  Memory allocated for batch                :',
     &              lMBatch, ' words,  ',(lMBatch+131071)/131072,' MB.'
              Write(6,'(A,I8)')        '  * Elements for batch:',lBatch
              Write(6,'(A,I8,A,I8,A)') '  * Number of batches :',nBatch,
     &                                 ' (',lBatch*nBatch,' elements)'
              Write(6,'(A,I8)')      '  * Residual elements :',leftBatch
            EndIf
            Write(6,*) ' Index matrix reloading.'
            Call XFlush(6)
          EndIf
C!
          filnam='NMAT'
          Call DaName_mf_wa(lNMAT,filnam)
          Call GetMem('nMAT','Allo','Inte',ipnMAT,nOsc*lBatch)
          Call ISCD_ReloadNMAT(nTabDim,iMaxYes,nOsc,lNMAT0,lNMAT,
     &         lBatch,nBatch,leftBatch,nIndex,
     &         iWork(ipnTabDim),iWork(ipnMAT0),iWork(ipnMAT))
          Call GetMem('nMat0','Free','Inte',ipnMat0,nOsc)
          Call Timing(CPTF4,CPE,TIOTF4,TIOE)                      ! CGGn
C!
C!      The remaining excitation matrices nInc & nDec
C!
          If (iPrint.GE.2) then
            Write(6,*)
            Write(6,*) ' nInc and nDec matrices generation.'
            Call XFlush(6)
          EndIf
          filnam='NINC'
          Call DaName_mf_wa(lNINC,filnam)
          filnam='NDEC'
          Call DaName_mf_wa(lNDEC,filnam)
          Call GetMem('nInc','Allo','Inte',ipnInc,nOsc*lBatch)
          Call GetMem('nDec','Allo','Inte',ipnDec,nOsc*lBatch)
          Call ISCD_MakenIncDec(n_max,iMaxYes,nOsc,lNMAT,lNINC,lNDEC,
     &                   lBatch,nBatch,leftBatch,nIndex,iWork(ipGraph2),
     &                   iWork(ipnMAT),iWork(ipnInc),iWork(ipnDec))
          Call GetMem('Graph2','Free','Inte',
     &                ipGraph2,(n_max+1)*(m_max+1)*nOsc)
          Call Timing(CPTF5,CPE,TIOTF5,TIOE)                      ! CGGn
C!
C!      The ISC rate
C!
          If (iPrint.GE.2) then
            Write(6,*)
            Write(6,*) ' Franck-Condon factors evaluation.'
            Call XFlush(6)
          EndIf
          Call GetMem('FCWind2','Allo','Real',ipFCWind2,nYes)
          Call ISCD_Rate(iPrint,nOsc,max_nOrd,iMx_nOrd,iMaxYes,nYes,
     &    dMinWind,lBatch,nBatch,leftBatch,nIndex,iWork(ipVibWind2),
     &    lNMAT0,lNMAT,lNINC,lNDEC,nTabDim,iWork(ipnTabDim),
     &    Work(ipC1),Work(ipC2),Work(ipW1),Work(ipW2),det0,det1,det2,
     &    Work(ipC),Work(ipW),Work(ipr01),Work(ipr02),Work(ipr00),
     &    m_max,n_max,max_dip,nnsiz,FC00,Work(ipFCWind2),dRho,
     &    mTabDim,iWork(ipmMat),iWork(ipmInc),iWork(ipmDec),
     &            iWork(ipnMat),iWork(ipnInc),iWork(ipnDec))
          Call Timing(CPTF6,CPE,TIOTF6,TIOE)                      ! CGGn
C!
          Call GetMem('FCWind2','Free','Real',ipFCWind2,nYes)
          Call GetMem('VibWind2','Free','Inte',ipVibWind2,nYes)
          Call GetMem('nDec','Free','Inte',ipnDec,nOsc*lBatch)
          Call GetMem('nInc','Free','Inte',ipnInc,nOsc*lBatch)
          Call GetMem('nMat','Free','Inte',ipnMAT,nOsc*lBatch)
          Call GetMem('nTabDim','Free','Inte',ipnTabDim,nTabDim+1)
C!
          Call DaClos(lNDEC)
          Call DaClos(lNINC)
          Call DaClos(lNMAT)
          Call DaClos(lNMAT0)
C!
        ENDIF
C!
        If (iPrint.GE.3) then
          Write(6,*)
          Write(6,'(A)') ' Timing informations (sec.):              '
          Write(6,'(A)') ' ========================================='
          Write(6,'(A)') ' MODULE:                       CPU Elapsed'
          Write(6,'(A,2F8.2)') ' Matrix evaluation        ',
     &                           CPTF1-CPTF0, TIOTF1-TIOTF0
          Write(6,'(A,2F8.2)') ' Level screening          ',
     &                           CPTF2-CPTF1, TIOTF2-TIOTF1
          Write(6,'(A,2F8.2)') ' Energy calculation       ',
     &                           CPTF3-CPTF2, TIOTF3-TIOTF2
          Write(6,'(A,2F8.2)') ' Matrix reloading         ',
     &                           CPTF4-CPTF3, TIOTF4-TIOTF3
          Write(6,'(A,2F8.2)') ' nInc and nDec generation ',
     &                           CPTF5-CPTF4, TIOTF5-TIOTF4
          Write(6,'(A,2F8.2)') ' Franck-Condon evaluation ',
     &                           CPTF6-CPTF5, TIOTF6-TIOTF5
          Write(6,'(A,2F8.2)') ' TOTAL                    ',
     &                           CPTF6-CPTF0, TIOTF6-TIOTF0
          Write(6,'(A)') ' -----------------------------------------'
          Write(6,*)
        EndIf

        Call GetMem('mDec','Free','Inte',ipmDec,(mTabDim+1)*nOsc)
        Call GetMem('mInc','Free','Inte',ipmInc,(mTabDim+1)*nOsc)
        Call GetMem('mMat','Free','Inte',ipmMat,(mTabDim+1)*nOsc)

        Call GetMem('nMaxQ','Free','Inte',ipnMaxQ,nOsc)
        Call GetMem('Base2','Free','Real',ipBase2,nOsc*nOsc)

        Call GetMem('r00','Free','Real',ipr00,nOsc)
        Call GetMem('r02','Free','Real',ipr02,nOsc)
        Call GetMem('r01','Free','Real',ipr01,nOsc)

        Call GetMem('Hess2','Free','Real',ipHess2,nOsc*nOsc)
        Call GetMem('Hess1','Free','Real',ipHess1,nOsc*nOsc)

        Call GetMem('BaseInv','Free','Real',ipBaseInv,nOsc*nOscOld)
        Call GetMem('Base','Free','Real',ipBase,nOscOld*nOsc)
        Call GetMem('Lambda','Free','Real',ipLambda,nOsc)

        Call GetMem('W2','Free','Real',ipW2,nOsc*nOsc)
        Call GetMem('W1','Free','Real',ipW1,nOsc*nOsc)
        Call GetMem('W','Free','Real',ipW,nOsc*nOsc)
        Call GetMem('C','Free','Real',ipC,nOsc*nOsc)
        Call GetMem('C2','Free','Real',ipC2,nOsc*nOsc)
        Call GetMem('C1','Free','Real',ipC1,nOsc*nOsc)

        Call GetMem('AtCoord','Free','Real',ipAtCoord,3*NumOfAt)
        Call Getmem('Gbis0','Free','Real',ip_Gbis0,ngdim**4)
        Call Getmem('Gbis2','Free','Real',ip_Gbis2,ngdim**4)
        Call Getmem('Gbis1','Free','Real',ip_Gbis1,ngdim**4)

        Call Getmem('Gprm0','Free','Real',ip_Gprm0,ngdim**3)
        Call Getmem('Gprm2','Free','Real',ip_Gprm2,ngdim**3)
        Call Getmem('Gprm1','Free','Real',ip_Gprm1,ngdim**3)

        Call GetMem('x_anharm2','Free','Real',ipx_anharm2,nOsc*nOsc)
        Call GetMem('x_anharm1','Free','Real',ipx_anharm1,nOsc*nOsc)
        Call GetMem('harmfreq2','Free','Real',ipharmfreq2,nOsc)
        Call GetMem('harmfreq1','Free','Real',ipharmfreq1,nOsc)

        Call GetMem('G0','Free','Real',ipG0,nOsc*nOsc)
        Call GetMem('G2','Free','Real',ipG2,nOsc*nOsc)
        Call GetMem('G1','Free','Real',ipG1,nOsc*nOsc)

        Call GetMem('r0','Free','Real',ipr0,nOsc)
        Call GetMem('r2','Free','Real',ipr2,nOsc)
        Call GetMem('r1','Free','Real',ipr1,nOsc)

        Call GetMem('n_plot','Free','Inte',ipn_plot,l_n_plot)
        Call GetMem('m_plot','Free','Inte',ipm_plot,l_m_plot)

        Call GetMem('NormModes','Free','Inte',ipNormModes,l_NormModes)
        Call GetMem('TranDipGrad','Free','Real',ipTranDipGrad,3*NumOfAt)

        lISC = .False.
c      Call GetMem('Test_F','LIST','INTE',iDum,iDum)               ! CGGt
        GoTo 999

      EndIf
CGGn -------------------------------------------------------------------

C!-----------------------------------------------------------------------!
C!
C!      Set up excitation matrices
C!
C!-----------------------------------------------------------------------!
C!
C!----      Calculate dimensions given max level of excitation for the
C!      different states.
      Call TabDim2_drv(m_max,nOsc,nvTabDim)
      mTabDim = nvTabDim-1
      Call TabDim2_drv(n_max,nOsc,nvTabDim)
      nTabDim = nvTabDim-1
C!
C!----      Set up mMat for L.
      max_mOrd = mTabDim
      Call TabDim2_drv(m_max-1,nOsc,nvTabDim)
      max_mInc = nvTabDim-1
      call GetMem('mMat','Allo','Inte',ipmMat,(mTabDim+1)*nOsc)

      call GetMem('mInc','Allo','Inte',ipmInc,(mTabDim+1)*nOsc)
      call GetMem('mDec','Allo','Inte',ipmDec,(mTabDim+1)*nOsc)
c Put dimensions into common block:
      mdim1=mTabDim
      mdim2=nOsc
      Call MakeTab2(
     &          m_max,max_mOrd,max_mInc,mTabDim,iWork(ipmMat),
     &    iWork(ipmInc),iWork(ipmDec),nOsc)
C!
C!----      Set up nMat for U.
      max_nOrd = nTabDim
      Call TabDim2_drv(n_max-1,nOsc,nvTabDim)
      max_nInc = nvTabDim-1
      call GetMem('nMat','Allo','Inte',ipnMat,(nTabDim+1)*nOsc)

      call GetMem('nInc','Allo','Inte',ipnInc,(nTabDim+1)*nOsc)
      call GetMem('nDec','Allo','Inte',ipnDec,(nTabDim+1)*nOsc)
c Put dimensions into common block:
      ndim1=nTabDim
      ndim2=nOsc
      nnsiz=ndim1
      Call MakeTab2(
     &          n_max,max_nOrd,max_nInc,nTabDim,iWork(ipnMat),
     &      iWork(ipnInc),iWork(ipnDec),nOsc)
C!
      If ( .not.MatEl ) Then
      l_TermMat_1=max_mOrd
      l_TermMat_2=max_nOrd
      Call GetMem('TermMat','Allo','Real',
     &      ipTermMat,(l_TermMat_1+1)*(l_TermMat_2+1))
      call GetMem('level1','Allo','Inte',iplevel1,nOsc)
      call GetMem('level2','Allo','Inte',iplevel2,nOsc)
      Do jOrd = 0,max_nOrd
      do iv=1,nOsc
      iWork(iplevel2+iv-1) =
     &      iWork(ipnMat+jOrd+(nTabDim+1)*(iv-1))
      enddo
      Do iOrd = 0,max_mOrd
      do iv=1,nOsc
      iWork(iplevel1+iv-1) =
     &      iWork(ipmMat+iOrd+(mTabDim+1)*(iv-1))
      enddo
      l_harm=nOsc
      Call TransEnergy(GE1,Work(ipx_anharm1),
     &      Work(ipharmfreq1),iWork(iplevel1),
     &              GE2,Work(ipx_anharm2),Work(ipharmfreq2),
     &      iWork(iplevel2),
     &      Work(ipTermMat+iOrd+(l_TermMat_1+1)*jOrd),l_harm)
      End Do
      End Do
      call GetMem('level1','Free','Inte',iplevel1,nOsc)
      call GetMem('level2','Free','Inte',iplevel2,nOsc)

      End If
C!
C!-----------------------------------------------------------------------!
C!
C!      Calculate vibrational Hessian if a variational calculation
C!      is chosen
C!
C!-----------------------------------------------------------------------!
C!
      l_l=1
      call GetMem('U1','Allo','Real',ipU1,l_l*l_l)
      call GetMem('E1','Allo','Real',ipE1,l_l)
      call GetMem('U2','Allo','Real',ipU2,l_l*l_l)
      call GetMem('E2','Allo','Real',ipE2,l_l)
c           l_l=1
      Call GetMem('OccNumMat1','Allo','Real',ipOccNumMat1,l_l*3)
      Call GetMem('OccNumMat2','Allo','Real',ipOccNumMat2,l_l*3)
      If ( MatEl ) Then       ! START: Vibrational Hessian (variational)

      call GetMem('U2','Free','Real',ipU2,l_l*l_l)
      call GetMem('E2','Free','Real',ipE2,l_l)
      call GetMem('U1','Free','Real',ipU1,l_l*l_l)
      call GetMem('E1','Free','Real',ipE1,l_l)

      Call GetMem('OccNumMat1','Free','Real',ipOccNumMat1,l_l*3)
      Call GetMem('OccNumMat2','Free','Real',ipOccNumMat2,l_l*3)

      If ( Forcefield ) Then
      nDimTot = max_mOrd+1
      Else
      nDimTot = 2*max_mOrd+2
      End If
      l_l=nDimTot
      call GetMem('U1','Allo','Real',ipU1,l_l*l_l)
      call GetMem('E1','Allo','Real',ipE1,l_l)
      call GetMem('U2','Allo','Real',ipU2,l_l*l_l)
      call GetMem('E2','Allo','Real',ipE2,l_l)

      Call GetMem('H1','Allo','Real',ipH1,nDimTot*nDimTot)
      Call GetMem('H2','Allo','Real',ipH2,nDimTot*nDimTot)
      Call GetMem('S1','Allo','Real',ipS1,nDimTot*nDimTot)
      Call GetMem('S2','Allo','Real',ipS2,nDimTot*nDimTot)
      Call GetMem('OccNumMat1','Allo','Real',ipOccNumMat1,l_l*3)
      Call GetMem('OccNumMat2','Allo','Real',ipOccNumMat2,l_l*3)

c            H1 = 0.0d0
      call dcopy_(nDimTot*nDimTot,0.0d0,0,Work(ipH1),1)
      call dcopy_(nDimTot*nDimTot,0.0d0,0,Work(ipS1),1)
      call dcopy_(nDimTot*nDimTot,0.0d0,0,Work(ipH2),1)
      call dcopy_(nDimTot*nDimTot,0.0d0,0,Work(ipS2),1)
c           S1 = 0.0d0
c            H2 = 0.0d0
c            S2 = 0.0d0
C!
C!----      Calculate inverse mass tensor and its derivatives in r0.
      call GetMem('Smat','Allo','Real',ipSmat,3*NumOfAt*nOscOld)
      call dcopy_(3*NumOfAt*nOscOld,0.0d0,0,Work(ipSmat),1)
C!      Call Int_To_Cart(InterVec,r0,AtCoord,NumOfAt,nOscOld,Mass)
      l_a=NumOfAt
      Call Int_To_Cart1(InterVec,Work(ipr0),
     &      Work(ipAtCoord),l_a,nOsc)
      Call CalcS(Work(ipAtCoord),InterVec,Work(ipSmat),
     &      nOscold,NumOfAt)
      Call CalcG(Work(ipG0),Mass,Work(ipSmat),nOscold,NumOfAt)
      If ( max_term.gt.2 ) Then
      dh = 1.0d-3
      Call CalcGprime(work(ip_gprm0),Mass,Work(ipr0),
     &      InterVec,Work(ipAtCoord),
     &      NumOfAt,dh,nOsc)
      dh = 1.0d-2
      Call CalcGdblePrime(work(ip_gbis0),Mass,Work(ipr0),
     &      InterVec,
     &      Work(ipAtCoord),NumOfAt,dh,nOsc)
C!
      call GetMem('T4','Allo','Real',ipT4,nOscOld**4)
      Call DGEMM_('T','T',
     &            nOscOld**2,nOsc,nOscOld,
     &            1.0d0,work(ip_gprm2),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,Work(ipT4),nOscOld**2)
      Call DGEMM_('T','T',
     &            nOsc*nOscOld,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,work(ip_gprm2),nOsc*nOscOld)
      Call DGEMM_('T','N',
     &            nOsc**2,nOsc,nOscOld,
     &            1.0d0,work(ip_gprm2),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(ipT4),nOsc**2)
      Call getmem('Gprm2','Free','Real',ip_Gprm2,ngdim**3)


      Call getmem('Gprm2','Allo','Real',ip_Gprm2,ngdim**3)
      call dcopy_(nOsc**3,Work(ipT4),1,work(ip_gprm2),1)
C!
      Call DGEMM_('T','T',
     &            nOscOld**3,nOsc,nOscOld,
     &            1.0d0,work(ip_gbis2),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,Work(ipT4),nOscOld**3)
      Call DGEMM_('T','T',
     &            nOsc*nOscOld**2,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,work(ip_gbis2),nOsc*nOscOld**2)
      Call DGEMM_('T','N',
     &            nOsc**2*nOscOld,nOsc,nOscOld,
     &            1.0d0,work(ip_gbis2),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(ipT4),nOsc**2*nOscOld)
      Call DGEMM_('T','N',
     &            nOsc**3,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,work(ip_gbis2),nOsc**3)
      Call DGEMM_('T','T',
     &            nOscOld**2,nOsc,nOscOld,
     &            1.0d0,work(ip_gprm1),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,Work(ipT4),nOscOld**2)
      Call DGEMM_('T','T',
     &            nOsc*nOscOld,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,work(ip_gprm1),nOsc*nOscOld)
      Call DGEMM_('T','N',
     &            nOsc**2,nOsc,nOscOld,
     &            1.0d0,work(ip_gprm1),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(ipT4),nOsc**2)
      Call getmem('Gprm1','Free','Real',ip_Gprm1,ngdim**3)
      Call getmem('Gprm1','Allo','Real',ip_Gprm1,ngdim**3)
      call dcopy_(nOsc**3,Work(ipT4),1,work(ip_gprm1),1)
C!
      Call DGEMM_('T','T',
     &            nOscOld**3,nOsc,nOscOld,
     &            1.0d0,work(ip_gbis1),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,Work(ipT4),nOscOld**3)
      Call DGEMM_('T','T',
     &            nOsc*nOscOld**2,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,work(ip_gbis1),nOsc*nOscOld**2)
      Call DGEMM_('T','N',
     &            nOsc**2*nOscOld,nOsc,nOscOld,
     &            1.0d0,work(ip_gbis1),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(ipT4),nOsc**2*nOscOld)
      Call getmem('Gbis1','Free','Real',ip_Gbis1,ngdim**4)
      Call getmem('Gbis1','Allo','Real',ip_Gbis1,ngdim**4)
      Call DGEMM_('T','N',
     &            nOsc**3,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,work(ip_gbis1),nOsc**3)
      Call DGEMM_('T','T',
     &            nOscOld**2,nOsc,nOscOld,
     &            1.0d0,work(ip_gprm0),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,Work(ipT4),nOscOld**2)
      Call DGEMM_('T','T',
     &            nOsc*nOscOld,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,work(ip_gprm0),nOsc*nOscOld)
      Call DGEMM_('T','N',
     &            nOsc**2,nOsc,nOscOld,
     &            1.0d0,work(ip_gprm0),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(ipT4),nOsc**2)
      Call getmem('Gprm0','Free','Real',ip_Gprm0,ngdim**3)
      ngdim=nosc
      Call getmem('Gprm0','Allo','Real',ip_Gprm0,ngdim**3)
      call dcopy_(nOsc**3,Work(ipT4),1,work(ip_gprm0),1)
C!
      Call DGEMM_('T','T',
     &            nOscOld**3,nOsc,nOscOld,
     &            1.0d0,work(ip_gbis0),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,Work(ipT4),nOscOld**3)
      Call DGEMM_('T','T',
     &            nOsc*nOscOld**2,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBaseInv),nOsc,
     &            0.0d0,work(ip_gbis0),nOsc*nOscOld**2)
      Call DGEMM_('T','N',
     &            nOsc**2*nOscOld,nOsc,nOscOld,
     &            1.0d0,work(ip_gbis0),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(ipT4),nOsc**2*nOscOld)
      Call getmem('Gbis0','Free','Real',ip_Gbis0,ngdim**4)
      Call getmem('Gbis0','Allo','Real',ip_Gbis0,ngdim**4)
      Call DGEMM_('T','N',
     &            nOsc**3,nOsc,nOscOld,
     &            1.0d0,Work(ipT4),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,work(ip_gbis0),nOsc**3)
      call GetMem('T4','Free','Real',ipT4,nOscOld**4)
      End If
      Call GetMem('temp','Allo','Real',iptemp,nOscOld*nOscOld)
      Call GetMem('temp2','Allo','Real',iptemp2,nOscOld*nOscOld)
C!
      call dcopy_(nOscOld*nOscOld,Work(ipG0),1,Work(iptemp),1)
      call dcopy_(nOscOld**2,0.0d0,0,Work(iptemp2),1)
      call dcopy_(nOscOld,1.0d0,0,Work(iptemp2),nOscOld+1)
      Call Dool(Work(iptemp),nOscOld,nOscOld,Work(iptemp2),
     &      nOscOld,nOscOld,det)
      Call DGEMM_('N','N',
     &            nOscOld,nOsc,nOscOld,
     &            1.0d0,Work(iptemp2),nOscOld,
     &            Work(ipBase),nOscOld,
     &            0.0d0,Work(iptemp),nOscOld)
      Call GetMem('temp2','Free','Real',iptemp2,nOscOld*nOscOld)
      Call GetMem('temp2','Allo','Real',iptemp2,nOsc*nOsc)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOscOld,
     &            1.0d0,Work(ipBase),nOscOld,
     &            Work(iptemp),nOscOld,
     &            0.0d0,Work(iptemp2),nOsc)
      call dcopy_(nOsc**2,0.0d0,0,Work(ipG0),1)
      call dcopy_(nOsc,1.0d0,0,Work(ipG0),nOsc+1)
      Call Dool(Work(iptemp2),nOsc,nOsc,Work(ipG0),
     &      nOsc,nOsc,det)
C!
      Call GetMem('temp','Free','Real',iptemp,nOscOld*nOscOld)
      Call GetMem('temp2','Free','Real',iptemp2,nOsc*nOsc)

      call GetMem('Smat','Free','Real',ipSmat,3*NumOfAt*nOscOld)
C!
C!----      Set up Hamilton matrix for the first state.
      If ( Forcefield ) Then
c            Base2 = 0.0d0
      call dcopy_(nOsc*nOsc,0.0d0,0,Work(ipBase2),1)
      Do i = 1,nOsc
      Work(ipBase2+i-1+nOsc*(i-1)) = 1.0d0
      End Do

      Call SetUpHmat2(energy1,energy2,Work(ipC1),
     &      Work(ipW1),det1,Work(ipr01),Work(ipr02),
     &      max_mOrd,max_nOrd,
     &      max_nOrd,max_mInc,max_nInc,max_nInc,iWork(ipmMat),
     &      iWork(ipnMat),iWork(ipmInc),
     &              iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &       Work(ipH1),Work(ipS1),Work(ipHess1),Work(ipG1),
     &       Work(ipBase2),Work(ipr01),nnsiz,
     &       nDimTot,nOsc)
      Call SolveSecEq(Work(ipH1),nDimTot,Work(ipU1),
     &      Work(ipS1),Work(ipE1))
      Write(6,'(20f10.1)')
     &      ((Work(ipE1+i-1)-Work(ipE1))*219474.63d0,i=2,nOsc)
      Call SetUpHmat2(energy1,energy2,Work(ipC2),
     &      Work(ipW2),det2,Work(ipr01),Work(ipr02),
     &      max_mOrd,max_nOrd,
     &              max_nOrd,max_mInc,max_nInc,max_nInc,
     &  iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),
     &              iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &         Work(ipH2),Work(ipS2),Work(ipHess2),Work(ipG2),
     &  Work(ipBase2),Work(ipr02),nnsiz,
     &       nDimTot,nOsc)
      Call SolveSecEq(Work(ipH2),nDimTot,Work(ipU2),
     &      Work(ipS2),Work(ipE2))
      Write(6,'(20f10.1)')
     &      ((Work(ipE2+i-1)-Work(ipE2))*219474.63d0,i=2,nOsc)
      Else
      Call SetUpHmat(energy1,Work(ipr1),iWork(ipipow),
     &      Work(ipvar),Work(ipyin1),Work(ipr00),trfName1,
     &             max_term,Work(ipC1),Work(ipW1),det1,
     &  Work(ipr01),Work(ipC2),Work(ipW2),det2,Work(ipr02),
     &             max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,
     &           max_nInc,iWork(ipmMat),iWork(ipnMat),iWork(ipmInc),
     &  iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &              Work(ipH1),Work(ipS1),Work(ipG1),Work(ipG2),
     &  Work(ipG0),work(ip_gprm1),work(ip_gprm2),
     &             work(ip_gprm0),
     &             work(ip_gbis1),work(ip_gbis2),work(ip_gbis0),
     &             Work(ipC),Work(ipW),det0,Mass,Work(ipr00),
     &  Work(ipBase),Work(ipr0),Work(ipr1),Work(ipr2),nnsiz,
     &   nterm,nvar,ndata,nosc,ndimtot,numofat)
      Call SolveSort(Work(ipH1),Work(ipU1),Work(ipS1),
     &      Work(ipE1),Work(ipW1),Work(ipW2),
     &      Work(ipW1),Work(ipC1),Work(ipC2),Work(ipC1),
     &      Work(ipr01),Work(ipr02),Work(ipr01),
     &              iWork(ipmInc),iWork(ipmDec),iWork(ipmMat),
     &   mdim1,mdim2,Work(ipOccNumMat1),nOsc,nDimTot)
      Call SetUpHmat(energy2,Work(ipr2),iWork(ipipow),
     &      Work(ipvar),Work(ipyin2),Work(ipr00),trfName2,
     &             max_term,Work(ipC1),Work(ipW1),det1,Work(ipr01),
     &  Work(ipC2),Work(ipW2),det2,Work(ipr02),
     &             max_mOrd,max_nOrd,max_nOrd,max_mInc,max_nInc,
     &             max_nInc,iWork(ipmMat),iWork(ipnMat),
     &  iWork(ipmInc),iWork(ipnInc),iWork(ipmDec),iWork(ipnDec),
     &             Work(ipH2),Work(ipS2),Work(ipG1),Work(ipG2),
     &  Work(ipG0),work(ip_gprm1),work(ip_gprm2),
     &             work(ip_gprm0),
     &             work(ip_gbis1),work(ip_gbis2),work(ip_gbis0),
     &             Work(ipC),Work(ipW),det0,Mass,Work(ipr00),
     &  Work(ipBase),Work(ipr0),Work(ipr1),Work(ipr2),nnsiz,
     &    nterm,nvar,ndata,nosc,ndimtot,numofat)
      Call SolveSort(Work(ipH2),Work(ipU2),Work(ipS2),
     &      Work(ipE2),Work(ipW1),Work(ipW2),
     &      Work(ipW2),Work(ipC1),Work(ipC2),Work(ipC2),
     &      Work(ipr01),Work(ipr02),Work(ipr02),
     &              iWork(ipnInc),iWork(ipnDec),iWork(ipnMat),
     &   ndim1,ndim2,Work(ipOccNumMat2),nOsc,nDimTot)
      Call GetMem('yin1','Free','Real',ipyin1,ndata)
      Call GetMem('yin2','Free','Real',ipyin2,ndata)

      End If
      Write(6,*) 'T0',T0*HarToRcm
      k2 = 1
      l_TermMat_1=nDimTot-1
      l_TermMat_2=nDimTot-1

      Call GetMem('TermMat','Allo','Real',
     &      ipTermMat,(l_TermMat_1+1)*(l_TermMat_2+1))

      Do jOrd = 0,nDimTot-1
      k1 = 1
      Do iOrd = 0,nDimTot-1
      Work(ipTermMat+iOrd+(l_TermMat_1+1)*jOrd) =
     &      T0+(Work(ipE2+k2-1)-Work(ipE1+k1-1))
      k1 = k1+1
      End Do
      k2 = k2+1
      End Do
      Call GetMem('H1','Free','Real',ipH1,nDimTot*nDimTot)
      Call GetMem('H2','Free','Real',ipH2,nDimTot*nDimTot)
      Call GetMem('S1','Free','Real',ipS1,nDimTot*nDimTot)
      Call GetMem('S2','Free','Real',ipS2,nDimTot*nDimTot)
c       Call GetMem('G0','Free','Real',ipG0,nOsc*nOsc)
      End If                    ! END: Vibrational Hessian (variational)
      Call getmem('Gprm0','Free','Real',ip_Gprm0,ngdim**3)
      Call getmem('Gbis0','Free','Real',ip_Gbis0,ngdim**4)
C!
      Call GetMem('Hess1','Free','Real',ipHess1,nOsc*nOsc)
      Call GetMem('Hess2','Free','Real',ipHess2,nOsc*nOsc)
      Call GetMem('AtCoord','Free','Real',ipAtCoord,3*NumOfAt)
      Call GetMem('G0','Free','Real',ipG0,nOsc*nOsc)
      Call GetMem('G1','Free','Real',ipG1,nOsc*nOsc)
      Call getmem('Gprm1','Free','Real',ip_Gprm1,ngdim**3)
      Call getmem('Gbis1','Free','Real',ip_Gbis1,ngdim**4)
      Call GetMem('G2','Free','Real',ipG2,nOsc*nOsc)
      Call getmem('Gprm2','Free','Real',ip_Gprm2,ngdim**3)
      Call getmem('Gbis2','Free','Real',ip_Gbis2,ngdim**4)
C!
C!-----------------------------------------------------------------------!
C!
C!      Intensity calculations
C!
C!-----------------------------------------------------------------------!
C!
C!PAM:      If max_dip=0, there are no transition dipole gradients, but the
C!      array TranDipGradInt must be formally all. anyway:
      If (max_dip.eq.0) then
      call GetMem('TranDipGradInt','Allo','Real',
     &      ipTranDipGradInt,3*nOsc)

      endif

      If ( Matel ) Then
      l_IntensityMat_1=ndimtot-1
      l_IntensityMat_2=ndimtot-1
      call GetMem('IntensityMat','Allo','Real',
     &      ipIntensityMat,(l_IntensityMat_1+1)*(l_IntensityMat_2+1))
c           call GetMem('IntensityMat','Allo','Real',ipIntensityMat,ndimtot*ndimtot)
      If ( Forcefield ) Then
      call dcopy_(nOsc*nOsc,0.0d0,0,Work(ipBase2),1)
c            Base2 = 0.0d0
      Do i = 1,nOsc
      Work(ipBase2+i-1+nOsc*(i-1)) = 1.0d0
      End Do
      Call Intensity2(Work(ipIntensityMat),Work(ipTermMat),
     &      T0,max_term,
     &        Work(ipU1),Work(ipU2),Work(ipE1),Work(ipE2),
     &  Work(ipC1),Work(ipW1),det1,Work(ipr01),
     &        Work(ipC2),Work(ipW2),det2,Work(ipr02),Work(ipC),
     &  Work(ipW),det0,Work(ipr00),m_max,n_max,max_dip,
     &        iWork(ipm_plot),iWork(ipn_plot),TranDip,
     &  Work(ipTranDipGradInt),
     &        Work(ipharmfreq1),Work(ipx_anharm1),
     &  Work(ipharmfreq2),Work(ipx_anharm2),
     &        Work(ipr0),Work(ipr1),Work(ipr2),Work(ipBase2),
     &      l_IntensityMat_1,l_IntensityMat_2,
     &  l_TermMat_1,l_TermMat_2, nOsc,
     &      nDimTot,l_n_plot,l_m_plot)
      call GetMem('TranDipGradInt','Free','real',
     &      ipTranDipGradInt,3*nOsc)
      Else
      Call Intensity(Work(ipIntensityMat),Work(ipTermMat),
     &      T0,max_term,iWork(ipipow),Work(ipvar),
     &              Work(ipt_dipin1),Work(ipt_dipin2),
     &   Work(ipt_dipin3),trfName1,Work(ipU1),Work(ipU2),
     &   Work(ipE1),Work(ipE2),
     &              Work(ipC1),Work(ipW1),det1,Work(ipr01),
     &   Work(ipC2),Work(ipW2),det2,
     &              Work(ipr02),work(ipC),Work(ipW),det0,
     &   Work(ipr00),
     &              m_max,n_max,max_dip,iWork(ipm_plot),
     &   iWork(ipn_plot),
     &              Work(ipharmfreq1),Work(ipharmfreq2),
     &   Work(ipr0),Work(ipr1),Work(ipr2),Work(ipBase),
     &      l_IntensityMat_1,l_IntensityMat_2,
     &   l_TermMat_1,l_TermMat_2,
     &      nOsc,nDimTot,nPolyTerm,ndata,nvar,MaxNumAt,
     &   l_n_plot,l_m_plot)
      Call GetMem('t_dipin1','Free','Real',ipt_dipin1,ndata)
      Call GetMem('t_dipin2','Free','Real',ipt_dipin2,ndata)
      Call GetMem('t_dipin3','Free','Real',ipt_dipin3,ndata)

      Call GetMem('var','Free','Real',ipvar,ndata*nvar)
      Call GetMem('ipow','Free','Inte',ipipow,nPolyTerm*nvar)

      End If
      Else
      l_IntensityMat_1=max_mOrd
      l_IntensityMat_2=max_nOrd

      call GetMem('IntensityMat','Allo','Real',
     &      ipIntensityMat,(l_IntensityMat_1+1)*(l_IntensityMat_2+1))
      Call IntForceField(Work(ipIntensityMat),
     &      Work(ipTermMat),T0,max_term,FC00,
     &              Work(ipC1),Work(ipW1),det1,Work(ipr01),
     &    Work(ipC2),Work(ipW2),det2,Work(ipr02),
     &              Work(ipC),Work(ipW),det0,Work(ipr00),
     &              m_max,n_max,max_dip,
     &              Trandip,Work(ipTranDipGradInt),
     &              Work(ipharmfreq1),Work(ipx_anharm1),
     &   Work(ipharmfreq2),Work(ipx_anharm2),
     &              iWork(ipmMat),iWork(ipmInc),iWork(ipmDec),
     &   iWork(ipnMat),iWork(ipnInc),iWork(ipnDec),OscStr,nsize,
     &  max_mOrd,max_nOrd, nDimTot,nOsc)
      call GetMem('TranDipGradInt','Free','Real',
     &      ipTranDipGradInt,3*nOsc)
      End If

C!
C!----      Write results to log.
      Write(6,*) ' Write intensity data to log file.'
      Call XFlush(6)
      Call WriteInt(Work(ipIntensityMat),Work(ipTermMat),
     &      iWork(ipmMat),iWork(ipnMat),
     &              Work(ipOccNumMat1),Work(ipOccNumMat2),MatEl,
     &   ForceField,Work(ipE1),Work(ipE2),
     &              T0,Work(ipharmfreq1),Work(ipharmfreq2),
     &   Work(ipx_anharm1),Work(ipx_anharm2),
     &      l_IntensityMat_1,l_IntensityMat_2,
     &      l_TermMat_1,l_TermMat_2,
     &      nDimTot,nOsc)
C!
      Call GetMem('C','Free','Real',ipC,nOsc*nOsc)
      Call GetMem('W','Free','Real',ipW,nOsc*nOsc)
      Call GetMem('W1','Free','Real',ipW1,nOsc*nOsc)
      Call GetMem('W2','Free','Real',ipW2,nOsc*nOsc)

      Call GetMem('C1','Free','Real',ipC1,nOsc*nOsc)
      Call GetMem('C2','Free','Real',ipC2,nOsc*nOsc)

      Call GetMem('r00','Free','Real',ipr00,nOsc)
      Call GetMem('r01','Free','Real',ipr01,nOsc)
      Call GetMem('r02','Free','Real',ipr02,nOsc)

      call GetMem('nInc','Free','Inte',ipnInc,(nTabDim+1)*nOsc)
      call GetMem('nDec','Free','Inte',ipnDec,(nTabDim+1)*nOsc)
      call GetMem('mInc','Free','Inte',ipmInc,(mTabDim+1)*nOsc)
      call GetMem('mDec','Free','Inte',ipmDec,(mTabDim+1)*nOsc)

      Call GetMem('harmfreq1','Free','Real',ipharmfreq1,nOsc)
      Call GetMem('harmfreq2','Free','Real',ipharmfreq2,nOsc)

      Call GetMem('x_anharm1','Free','Real',ipx_anharm1,nOsc*nOsc)
      Call GetMem('x_anharm2','Free','Real',ipx_anharm2,nOsc*nOsc)
      Call GetMem('TermMat','Free','Real',
     &  ipTermMat,(l_TermMat_1+1)*(l_TermMat_2+1))
      call GetMem('IntensityMat','Free','Real',
     &  ipIntensityMat,(l_IntensityMat_1+1)*(l_IntensityMat_2+1))

      call GetMem('U1','Free','Real',ipU1,l_l*l_l)
      call GetMem('E1','Free','Real',ipE1,l_l)
      call GetMem('U2','Free','Real',ipU2,l_l*l_l)
      call GetMem('E2','Free','Real',ipE2,l_l)
      call GetMem('mMat','Free','Inte',ipmMat,(mTabDim+1)*nOsc)
      call GetMem('nMat','Free','Inte',ipnMat,(nTabDim+1)*nOsc)

      Call GetMem('OccNumMat1','Free','Real',ipOccNumMat1,l_l*3)
      Call GetMem('OccNumMat2','Free','Real',ipOccNumMat2,l_l*3)

      Call GetMem('m_plot','Free','Inte',ipm_plot,l_m_plot)
      Call GetMem('n_plot','Free','Inte',ipn_plot,l_n_plot)

      Call GetMem('NormModes','Free','Inte',
     &    ipNormModes,l_NormModes)
      Call GetMem('Base','Free','Real',ipBase,nOscOld*nOsc)
      Call GetMem('BaseInv','Free','Real',
     &      ipBaseInv,nOsc*nOscOld)

c I don't understand why is it here??
      Call GetMem('TranDipGrad','Free','Real',
     &  ipTranDipGrad,3*NumOfAt)
      Call GetMem('Base2','Free','Real',ipBase2,nOsc*nOsc)
      Call GetMem('Lambda','Free','Real',ipLambda,nOsc)
      Call GetMem('r1','Free','Real',ipr1,nOsc)
      Call GetMem('r2','Free','Real',ipr2,nOsc)
      Call GetMem('r0','Free','Real',ipr0,nOsc)
C!
 999  Continue
C!
      iReturn=0
      End
