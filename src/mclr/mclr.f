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
* Copyright (C) 1997, Anders Bernhardsson                              *
************************************************************************
      subroutine MCLR(ireturn)
************************************************************************
*                                                                      *
*               #     #  #####     #       ######                      *
*               ##   ## #     #    #       #     #                     *
*               # # # # #          #       #     #                     *
*               #  #  # #          #       ######                      *
*               #     # #          #       #   #                       *
*               #     # #     #    #       #    #                      *
*               #     #  #####     ####### #     #                     *
*                                                                      *
*     A linear response function program for general RASSCF states     *
*                                                                      *
*      OK right now we can just handle CASSCF, but the only thing      *
*      we need is yet another stupid PhD student who wants to code     *
*      the preconditioning elements between active and active orbitals *
*                                                                      *
*      The spin correction part of the code is a experimental          *
*      test platform, it is not finished and is not working,           *
*      If someone wants to help me finishing it please contact me      *
*                                                                      *
*                                                                      *
*   "Real programmers don't comment their codes.                       *
*    It was hard to write, and it must be hard to understand..."       *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "warnings.fh"
#include "WrkSpc.fh"
#include "machine.fh"
#include "SysDef.fh"

#include "blksize.fh"
#include "sa.fh"
#include "Pointers.fh"
#include "glbbas_mclr.fh"
#include "lbbas1.fh"
#include "detdim.fh"
#include "csfbas_mclr.fh"
#include "dmrginfo_mclr.fh"
#include "csfsd.fh"

      Logical Reduce_Prt
      External Reduce_Prt
*
*     first a few things for making Markus Happy
      logical converged(8)
      logical DoCholesky
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu1,TWall1)
*                                                                      *
************************************************************************
*                                                                      *
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=iPL-1
*                                                                      *
************************************************************************
*                                                                      *
*
      call qenter('MCLR')
*
      iAllo=0
c      idp=rtoi
      nrec=MBl_wa/rtob
*
      Call DecideOnCholesky(DoCholesky)
      Call get_iScalar('nSym',nSymX)

      If (DoCholesky .and. nSymX.gt.1) then
       write(6,*)'** Cholesky or RI/DF not implemented with symmetry **'
       Call Quit(_RC_INPUT_ERROR_)
      EndIf

*
      Call Init_Data
*
*

      doDMRG = .false.
      doMCLR = .false.
#ifdef _DMRG_
!      !> read dmrg parameters for mclr part
!      call read_dmrg_parameter_for_mclr()
!      if(doDMRG)then
!        doMCLR=.true.
!        write(*,*)"ndets_RGLR : ",ndets_RGLR
!        write(*,*)"ncsfs_RGLR : ",ncsfs_RGLR
!        write(*,*)"nstates_RGLR ",nstates_RGLR
!        write(*,*)"RGras2 : ",RGras2
!        write(*,*)"LRras2 : ",LRras2
!       open(unit=117,file="mclr_dets.initial")
!      end if
#endif

*                                                                      *
************************************************************************
*                                                                      *
*     open files
*
      Call OpnFls_MCLR(iPL)
      Call IpInit()

*                                                                      *
************************************************************************
*                                                                      *
*     Read input
*
      Call InpCtl_MCLR(iPL)

*                                                                      *
************************************************************************
*                                                                      *
*     Transform integrals and calculate fock matrixes etc.
*
      if(doDMRG)then  ! yma
        call dmrg_spc_change_mclr(RGras2(1:8),nAsh)
        call dmrg_spc_change_mclr(RGras2(1:8),nrs2)
      end if
      ntAsh=0
      ntAtri=0
      ntAsqr=0
      nnA=0
      Do iSym=1,nSym
         ntAsh=ntAsh+nAsh(iSym)
         ntAtri=ntAtri+nAsh(iSym)*(nAsh(iSym)+1)/2
         ntAsqr=ntAsqr+nAsh(iSym)*nAsh(iSym)
         nA(iSym)=nna
         nnA=nnA+nAsh(isym)
      end do

      Call Start_MCLR()

*
      nisp=max(8,nDisp)
*

*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
*                                                                      *
************************************************************************
*                                                                      *
*     File pointers
*
      Call GetMem('KapFile','ALLO','INTE',ifpK,nisp)
      Call ICOPY(nisp,[-1],0,iWork(ifpk),1)
      Call GetMem('SigFile','ALLO','INTE',ifpS,nisp)
      Call ICOPY(nisp,[-1],0,iWork(ifps),1)
      Call GetMem('RHSFile','ALLO','INTE',ifpRHS,nisp)
      Call ICOPY(nisp,[-1],0,iWork(ifprhs),1)
      If (iMethod.eq.2) Then
         Call GetMem('CIFile','ALLO','INTE',ifpCI,nisp)
         Call ICOPY(nisp,[-1],0,iWork(ifpci),1)
         Call GetMem('SCFile','ALLO','INTE',ifpSC,nisp)
         Call ICOPY(nisp,[-1],0,iWork(ifpsc),1)
         Call GetMem('RHCCIe','ALLO','INTE',ifpRHSCI,nisp)
         Call ICOPY(nisp,[-1],0,iWork(ifprhsci),1)
      Else
         ifpCI=1
         ifpSC=1
         ifpRHSCI=1
      End If

*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
*     Calculate response
*
*     Output is stored on disk
*
*

      If (SPINPOL) Then
         Call WfCtl_SP(iWork(ifpK),iWork(ifpS),iWork(ifpCI),
     &                 iWork(ifpSC),iWork(ifpRHS),iWork(ifpRHSCI))
*     Else if (ELECHESS) Then
*        Call WfCtl_PCG(iWork(ifpK),iWork(ifpS),iWork(ifpCI),
*                       iWork(ifpSC),iWork(ifpRHS),iWork(ifpRHSCI))
*        Call Abend()
      Else if(iMCPD) Then!pdft
         Call WfCtl_PDFT(iWork(ifpK),iWork(ifpS),iWork(ifpCI),
     &                 iWork(ifpSC),iWork(ifpRHS),
     &                 converged,iPL)
      Else if (SA) Then
         Call WfCtl_SA(iWork(ifpK),iWork(ifpS),iWork(ifpCI),
     &                 iWork(ifpSC),iWork(ifpRHS),
     &                 converged,iPL)
      Else If (TimeDep) Then
         Call WfCtl_td(iWork(ifpK),iWork(ifpS),iWork(ifpCI),
     &                 iWork(ifpSC),iWork(ifpRHS),iWork(ifpRHSCI),
     &                 converged)
      Else
         Call WfCtl_Hess(iWork(ifpK),iWork(ifpS),iWork(ifpCI),
     &                   iWork(ifpSC),iWork(ifpRHS),iWork(ifpRHSCI),
     &                   converged)
      End If


*                                                                      *
************************************************************************
*                                                                      *
*     Contract response to hessian etc
*
      If(.not.(TwoStep.and.(StepType.eq.'RUN1'))) Then
         If (PT2.or.SA.or.iMCPD) Then
            Call Out_PT2(iWork(ifpK),iWork(ifpCI))
         Else If (TimeDep) Then
            Call Output_td(iWork(ifpK),iWork(ifpS),
     &                     iWork(ifpCI),iWork(ifpSC),
     &                     iWork(ifpRHS),iWork(ifpRHSCI),converged)
         Else
            Call Output_mclr(iWork(ifpK),iWork(ifpS),
     &                       iWork(ifpCI),iWork(ifpSC),
     &                       iWork(ifpRHS),iWork(ifpRHSCI),converged)
            If (mckinley) Call isoloop(Double)
         End If
*
         If (RASSI)   Call OutRAS   (iWork(ifpK),iWork(ifpCI))
*
         If (TimeDep) Call OutRAS_td(iWork(ifpK),iWork(ifpCI))
      End If
*                                                                      *
************************************************************************
*                                                                      *

*     Deallocate memory
*
      nDum=1
*
*     Arrays in csf.f
      If (iMethod.eq.2) Then
         CALL GETMEM('KICTS ','Free','INTEGER',KICTS(1),nDum)
         CALL GETMEM('KICONF','Free','INTEGER',KICONF(1),nDum)
         CALL GETMEM('KDTOC','Free','REAL   ',KDTOC,nDum)
         CALL GETMEM('KCFTP','Free','INTEGER',KCFTP,nDum)
         CALL GETMEM('KDFTP','Free','INTEGER',KDFTP,nDum)
      End if

*     Free arrays allocated by memstr.f
      If (iMethod.eq.2) Call FreeStr
*     Arrays in inpone.f
      Call GetMem('ONEHAM','Free','REAL',kint1,nDum)
*     Arrays in incsfsd.f
      If (iAllo.eq.1) Then
         Call Getmem('KICONF2','Free','INTEGER',kiconf(2),nDum)
         Call Getmem('KICTS2','Free','INTEGER',Kicts(2),nDum)
      End If
*     Arrays in detctl.f
      If (iMethod.eq.2) Then
         Call Getmem('TwoOff','Free','INTE',KpINT2,nDum)
         Call Getmem('OneOff','Free','INTE',KpINT1,nDum)
      End if

*     Array in rdab.f
      If (iMethod.eq.1) Then
         Call Free_Work(ipCMO)
      End If
*     Arrays in rdjobip_td.f and rdjobiph.f
      If (iMethod.eq.2) Then
         Call GetMem(' G2 ','Free','Real',ipG2t,nDum)
         If (Timedep) Call GetMem(' G2 ','Free','Real',ipG2sq,nDum)
         Call GetMem(' G1 ','Free','Real',ipG1t,nDum)
         Call Free_Work(ipCMO)
      End If

*     Arrays allocated in stpert.f
      Call GetMem('CONN','Free','Real',ipHss,nDum)
*     Arrays allocated in fckmat.f
      Call GetMem('FASQMO','Free','Real',ipFAMO,ndens2)
      If (iMethod.eq.2) Call GetMem('K2Int','Free','Real',k2int,nDum)
      Call GetMem('fisqMO','Free','Real',ipfiMO,ndens2)
      Call GetMem('f0sqMO','Free','Real',ipf0sqMO,ndens2)

      If (iMethod.eq.2) Then
         Call GetMem('RHCCIe','Free','INTE',ifpRHSCI,nisp)
         Call GetMem('SCFile','Free','INTE',ifpSC,nisp)
         Call GetMem('CIFile','Free','INTE',ifpCI,nisp)
      End If
      Call GetMem('RHSFile','Free','INTE',ifpRHS,nisp)
      Call GetMem('SigFile','Free','INTE',ifpS,nisp)
      Call GetMem('KapFile','Free','INTE',ifpK,nisp)

*
*     Close files
*
      Call ClsFls_MCLR()
*
      If (NewCho) Then
        Call Cho_X_Final(irc)
        Do i=1,nSym
          Call DACLOS(LuAChoVec(i))
          Call DACLOS(LuIChoVec(i))
        End Do
        Call DACLOS(LuChoInt(1))
        Call DACLOS(LuChoInt(2))
        nOrbBas =0
        Do iSym=1,nSym
          nOrbBas  = nOrbBas  + nOrb(iSym)*nBas(iSym)
        End Do
        Call GetMem('CMO_inv','Free','Real',ip_CMO_inv,nOrbBas)
      End If

      If(TwoStep.and.(StepType.eq.'RUN1')) irc=ipclose(-1)
*                                                                      *
************************************************************************
*                                                                      *
      If (.not.fail) Then
         If (iPL.ge.2) Then
         Write(6,*)
         Write(6,'(6X,A,A)') 'The response parameters are ',
     &                      'written to the file RESP.'
         End If
         ireturn=_RC_ALL_IS_WELL_
      Else
         ireturn=_RC_NOT_CONVERGED_
      End If
      Call qStat(' ')
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu3,TWall3)
      If (iPL.ge.3) Then
        Write(6,*)
        Write(6,'(2X,A)') 'Timings'
        Write(6,'(2X,A)') '-------'
        Write(6,*)
        Write(6,'(2X,A)')
     &      '- - - - - - - - - - - - - - - - -'//
     &      ' - - - - - - - - - - - - - - - - -'
        Write(6,'(2X,A,T44,A,A,A)')
     &      ' ',' ','    CPU time','     elapsed'
        Write(6,'(2X,A,T44,A,2F12.2)')
     &      '1) Initialization',':',TCpu2-TCpu1,TWall2-TWall1
        Write(6,'(2X,A,T44,A,2F12.2)')
     &      '2) Response calculation',':',TCpu3-TCpu2,TWall3-TWall2
        Write(6,'(2X,A)')
     &      '- - - - - - - - - - - - - - - - -'//
     &      ' - - - - - - - - - - - - - - - - -'
        Write(6,'(2X,A,T44,A,2F12.2)')
     &      'Total',':',TCpu3-TCpu1,TWall3-TWall1
        Write(6,'(2X,A)')
     &      '- - - - - - - - - - - - - - - - -'//
     &      ' - - - - - - - - - - - - - - - - -'

      EndIf
*                                                                      *
************************************************************************
*                                                                      *
      If (iPL.ge.3) Then
          Call FastIO('STATUS')
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call QExit('MCLR')
      Return
      End
