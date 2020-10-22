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
* Copyright (C) 1989-1992, Roland Lindh                                *
*               1990, IBM                                              *
*               1995, Anders Bernhardsson                              *
************************************************************************
      subroutine McKinley(ireturn)
************************************************************************
*                                                                      *
*  Object: Driver for the one and two electron integral second order   *
*          derivative program McKinley.                                *
*                                                                      *
*  Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA     *
*          July '89 - May '90                                          *
*                                                                      *
*          Roland Lindh, Dept. of Theoretical Chemistry, University of *
*          Lund, SWEDEN. Modified to gradient calculations September   *
*          1991 - February 1992.                                       *
*                                                                      *
*          Anders Bernhardsson, Dept. of Theoretical Chemistry,        *
*          University of  Lund, SWEDEN.                                *
*          Modified to  second order derivatives October '94 -         *
*          '95                                                         *
************************************************************************
      use Real_Spherical
      use Basis_Info
      use Temporary_Parameters
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "disp2.fh"
#include "cputime.fh"
#include "print.fh"
#include "etwas.fh"
cpcm_solvent
#include "rctfld.fh"
cpcm_solvent end
c      Parameter (nLines=12)
      Character*120 Lines
      Logical DoRys, Run_MCLR
      Real*8, Allocatable:: Hess(:), Temp(:), GradN(:)
#include "warnings.fh"
*                                                                      *
************************************************************************
*                                                                      *
*     Call McKinley_banner()
      Call CWTime(TCpu1,TWall1)
      iRout=1
      call dcopy_(9,[0.0d0],0,CpuStat,1)
*                                                                      *
************************************************************************
*                                                                      *
*     Print program header
*                                                                      *
************************************************************************
*                                                                      *
c      Lines(1)=_MOLCAS_VERSION_
c#ifdef _DEMO_
c      Lines(2)='DEMO VERSION'
c#else
c      Lines(2)=' '
c#endif
c      Lines(3)=' '
c      Lines(4)=Vrsn
c      Lines(5)='A Vectorized Direct Integral Program for derivatives'
c      Lines(6)='of Cartesian and Spherical Harmonic Gaussians'
c      Lines(7)='Written by Anders Bernhardsson and Roland Lindh '
c      Lines(8)='Backtransformation of the 2nd order density matrix '//
c     &         'from MO to SO by Per-AAke Malmqvist'
c      Lines(9)='Dept. of Theoretical Chemistry, '//
c     &          'Chemical Centre, Lund (Sweden)'
c      Lines(10)=' '
c      Lines(11)=' '
c      Lines(12)='Compiled at '//
c     &           _BUILD_DATE_
c      lLine=Len(Lines(1))
C     Call Banner(Lines,nLines,lLine)
*                                                                      *
************************************************************************
*                                                                      *
*     Set error conditions
*
      Call XuFlow()
      Call ErrSet(209,1,1,2,1,209)
*                                                                      *
************************************************************************
*                                                                      *
*     Check if a numerical procedure will be used
*
      Call SuperMac()
*                                                                      *
************************************************************************
*                                                                      *
*
*     Get the input information as Seward dumped on INFO.
*     Set up some info
*     Read input
*
      nDiff=2
      DoRys=.True.
      Call IniSew(DoRys,nDiff)
cpcm_solvent
c check if there is a reaction field
c     write(6,*)'In mckinley PCM',pcm
      Call Init_RctFld(.False.,iCharge_ref)
cpcm_solvent end
      nsAtom=0
      Do  iCnttp = 1, nCnttp
            nsAtom=nsAtom+dbsc(iCnttp)%nCntr
      End Do
      Call Inputh(Run_MCLR)
      iPrint=nPrint(iRout)
      nGrad=0
      Do i=0,nIrrep-1
        nGrad=nGrad+lDisp(i)
      End Do
      Call OpnFls_Mckinley()
*                                                                      *
************************************************************************
*                                                                      *
*     Allocate area for hessian etc
*
      nHess=nGrad*(nGrad+1)/2
*
      Call mma_allocate(Hess,nHess,Label='Hess')
      Hess(:)=Zero
      Call mma_allocate(Temp,nHess,Label='Temp')
      Temp(:)=Zero
*                                                                      *
************************************************************************
*                                                                      *
*    Calculate the second order derivatives of the one electron        *
*    integrals and contract with D.                                    *
*                                                                      *
************************************************************************
*                                                                      *
      If (lHss) Then
         If (iPrint.ge.6) Then
         Write(6,*)
         Write(6,'(A,A,A)')
     &    'The 2nd order derivatives of the one-electron',
     &    ' integrals are calculated and contracted with',
     &    ' the one-electron density matrix. '
         Write(6,*)
         End If
         Call Timing(dum1,Time,dum2,dum3)
         Call Drvh2(Hess,Temp,nHess,show)
         Call DrvEtc(nGrad)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*        Compute contribution from the nuclear repulsion.              *
*                                                                      *
************************************************************************
*                                                                      *
      If (lHss) Then
         Call DrvN2(Temp,nGrad)
         If (SHOW) Call HssPrt(Temp,nHess)
         Call DaXpY_(nHess,One,Temp,1,Hess,1)
         If (Show) Call HssPrt(Hess,nHess)
      End If
      If (lGrd) Then
          Call mma_allocate(GradN,nGrad,Label='GradN')
          Call DrvN1_mck(GradN,nGrad)
          iopt=0
          irc=-1
          Call dWrMCK(iRC,iOpt,'NUCGRAD',1,GradN,1)
          If (irc.ne.0) Call SysAbendMsg('mckinley','Error in writing',
     &                                   'Option=NUCGRAD')
          Call mma_deallocate(GradN)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*       Calculate the first order derivatives of the one electron      *
*       integrals and store on disc in file MCKINT                     *
*                                                                      *
************************************************************************
*                                                                      *
      If (iPrint.ge.6) Then
      Write(6,*)
      Write(6,'(A,A)')
     &          'The 1st order derivatives of the one-electron ',
     &          'integrals are calculated and stored on disk'
      Write(6,*)
      End If
      Call Drvh1_mck(nGrad,Nona)
*
************************************************************************
*                                                                      *
*      Calculate two electron integrals. First order is contracted     *
*      to Fock matrixes and MO (IJKl) on the fly. Second order         *
*      derivatives are contracted with P.                              *
*      Derivatives are stored in MCKINT.                               *
*                                                                      *
************************************************************************
*                                                                      *
      nhess=ngrad*(ngrad+1)/2
      Call Timing(dum1,Time,dum2,dum3)
      CPUStat(nOneel)=CPUStat(nOneel)+Time
      If (.Not.Onenly) Then
*
          nIsh(:)=0
          nAsh(:)=0
*
          Call PrepP()
*
          iOpt = 0
          iRC = -1
          Call iWrMck(iRC,iOpt,'NISH',1,nIsh,iDummer)
          If (iRC.ne.0) Then
             Write (6,*) 'Mckinley: Error writing to MckInt!'
             Call Abend()
          End If
          iOpt = 0
          iRC = -1
          Call iWrMck(iRC,iOpt,'NASH',1,nAsh,iDummer)
          If (iRC.ne.0) Then
             Write (6,*) 'Mckinley: Error writing to MckInt!'
             Call Abend()
          End If
*
*

          Call Drvg2(Temp,nhess, lGrd,lHss)
*
          Call CloseP
*
          If (lHss) Then
             Call GADSum(Temp,nHess)
             Call DScal_(nHess,Half,Temp,1)
             If (Show) Call HssPrt(Temp,nHess)
*
*----------- Accumulate contribution to the hessian!
*
             Call DaXpY_(nhess,One,Temp,1,Hess,1)
*
             If (Show) Then
                Call Banner('Complete static Hessian',1,23)
                Call HssPrt(Hess,nHess)
             End If
             Call WrHDsk(Hess,ngrad)
          End If
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*...  Close 'MCKINT' file
      iRc=-1
      iOpt=0
      Call ClsMck(iRC,iOpt)
      If ( iRc.ne.0 ) Then
         Write (6,*) 'McKinley: Error closing MCKINT!'
         Call Abend()
      End If
      Call mma_deallocate(Temp)
      Call mma_deallocate(Hess)
*
      Call ClsSew
*                                                                      *
************************************************************************
*                                                                      *
*     Epilogue
*
      Lines='All data is written to disk, and could be accessed '//
     &         'through the MCLR or RASSI program.'
      lLine=Len(Lines)
      Call Banner(Lines,1,lLine)
*
      Call CWTime(TCpu2,TWall2)
      Call SavTim(5,TCpu2-TCpu1,TWall2-TWall1)
*
C     Call DaTimm
      Call Timing(Time,dum,dum,dum)
      CPUStat(nTotal)=Time
      If (iPrint.ge.6) Call Sttstc
      If (Test) Then
         ireturn=_RC_INPUT_ERROR_
      Else
         Call Request_MCLR_Run(Run_MCLR,ireturn,iPrint)
      End If
*
      Return
      End
      Subroutine Request_MCLR_Run(Run_MCLR,ireturn,iPrint)
      Logical Run_MCLR
      Character*16 StdIn
#include "warnings.fh"
*
      If (Run_MCLR) Then
*
*        McKinley will automatically generate the input for MCLR
*        and signal to AUTO (iRC=2) to run the input file Stdin.x.
*
         If (iPrint.ge.6) Then
         Write (6,*)
         Write (6,*)
     &     ' McKinley requests the MCLR module to be executed!'
         Write (6,*)
         End If
*
         LuInput=11
         LuInput=IsFreeUnit(LuInput)
         Call StdIn_Name(StdIn)
         Call Molcas_Open(LuInput,StdIn)
         Write (LuInput,'(A)') ' &MCLR &End'
         Write (LuInput,'(A)') 'End of Input'
         Close(LuInput)
         ireturn=_RC_INVOKED_OTHER_MODULE_
      Else
         ireturn=_RC_ALL_IS_WELL_
      End if
*
      Return
*
      Write (6,*)
      Write (6,*) ' Error opening Stdin.x'
      Write (6,*)
      Call Quit(_RC_INPUT_EMIL_ERROR_)
*
      End
