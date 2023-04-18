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
      Use Basis_Info, only: Basis_Info_Free
      Use Center_Info, only: Center_Info_Free
      use Symmetry_Info, only: Symmetry_Info_Free
      use Arrays, only: Hss, FAMO, FAMO_SpinP, FAMO_SpinM, SFock,
     &                  G2mm, G2mp, G2pp, Fp, Fm, G1p, G1m,
     &                  CMO_Inv, CMO,
     &                  Int1, pINT1, INT2, pINT2, G2t, G2sq, G1t,
     &                  FIMO, F0SQMO
      use Str_Info, only: DFTP, CFTP, DTOC, CNSM
      use negpre, only: SS
      use PDFT_Util, only :Do_Hybrid,WF_Ratio,PDFT_Ratio
*     Added for CMS NACs
      use Definitions, only: iwp, u6
      Implicit Real*8 (a-h,o-z)
#include "Input.fh"
#include "warnings.h"
#include "stdalloc.fh"
#include "machine.fh"
#include "SysDef.fh"

#include "sa.fh"
#include "Pointers.fh"
#include "detdim.fh"
#include "dmrginfo_mclr.fh"
#include "csfsd.fh"
      Integer, Allocatable:: ifpK(:), ifpS(:), ifpRHS(:),
     &            ifpCI(:), ifpSC(:), ifpRHSCI(:)

#include "Files_mclr.fh"

      Logical Reduce_Prt
      External Reduce_Prt
      Integer get_MBl_wa
      External get_MBl_wa
*
*     first a few things for making Markus Happy
      logical converged(8)
      logical DoCholesky

* Additional things for CMS-NACs Optimization
      Character*8 Method
      integer(kind=iwp) :: LuInput, istatus, LuSpool2
      character(len=16) :: StdIn
      character(len=180) :: Line
      character(len=128) :: FileName
      logical(kind=iwp) :: Exists
      integer(kind=iwp), external :: isFreeUnit
      Logical :: CalcNAC_Opt = .False., MECI_via_SLAPAF = .False.

*   This used to be after the CWTIME() functional call                 *
************************************************************************
*                                                                      *
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=iPL-1
*                                                                      *
************************************************************************
*                                                                      *
!This is where you should put the information for CMS
      Call Get_cArray('Relax Method',Method,8)
      if(Method.eq.'MSPDFT  ') then
          Call Get_iArray('cmsNACstates    ', cmsNACstates,2)
          Call Get_iArray('NACstatesOpt    ', NACstates,2)
          Call Get_lscalar('CalcNAC_Opt     ', CalcNAC_Opt)
          call Get_lscalar('MECI_via_SLAPAF ', MECI_via_SLAPAF)

          if(MECI_via_SLAPAF.eqv..FALSE.) then
              if ( (cmsNACstates(1).ne.NACstates(1)).or.
     &    (cmsNACstates(2).ne.NACstates(2)) ) Then
                    NACstates(1) = cmsNACstates(1)
                    NACstates(2) = cmsNACstates(2)
              end if
          end if

          if ( (cmsNACstates(1).ne.NACstates(1)).or.
     & (cmsNACstates(2).ne.NACstates(2)) ) Then
             if (iPL >= 3) then
                 write(u6,*)
                 write(u6,*) 'MS-PDFT Potentials for root(s)'
                 write(u6,*) cmsNACstates(1), cmsNACstates(2)
                 write(u6,*) 'MCLR Lag. Mult. for root'
                 write(u6,*) NACstates(1), NACstates(2)
                 write(u6,*) 'MS-PDFT roots and MCLR roots do not match'
                 write(u6,*) 'MCLR requests MCPDFT to be run before'
                 write(u6,*) 'it states again'
                 write(u6,*)
             end if

             LuInput = 11
             LuInput = IsFreeUnit(LuInput)
             call StdIn_Name(StdIn)
             call Molcas_open(LuInput,StdIn)

             write(LuInput,'(A)') '>ECHO OFF'
             write(LuInput,'(A)') '>export MCLR_OLD_TRAP=$MOLCAS_TRAP'
             write(LuInput,'(A)') '>export MOLCAS_TRAP=ON'

             write(LuInput,'(A)') ' &MCPDFT &END'
             write(LuInput,'(A)') ' KSDFT=T:PBE'
             write(LuInput,'(A)') ' MSPDft'
             write(LuInput,'(A)') ' GRAD'
             if (NACstates(2).ne.0) then
               write(LuInput,'(A)') 'NAC'
               write(LuInput,'(I5,1X,I5)') NACstates(1),NACstates(2)
               if(CalcNAC_Opt) Write(LuInput,'(A)') 'MECI'
             end if
             write(LuInput,'(A)') 'End of Input'
             write(LuInput,'(A)') ' '

             FileName = 'MCLRINP'
             call f_inquire(Filename,Exists)

             if (Exists) then
               LuSpool2 = 77
               LuSpool2 = IsFreeUnit(LuSpool2)
               call Molcas_Open(LuSpool2,Filename)
               do
                  read(LuSpool2,'(A)',iostat=istatus) Line
                  if (istatus > 0) call Abend()
                  if (istatus < 0) exit
                  write(LuInput,'(A)') Line
               end do
               close(LuSpool2)
             else
               write(LuInput,'(A)') ' &MCLR &End'
             end if

             write(LuInput,'(A)') '>export MOLCAS_TRAP=$MCLR_OLD_TRAP'
             write(LuInput,'(A)') '>ECHO ON'
             close(LuInput)
             call Finish(_RC_INVOKED_OTHER_MODULE_)
           End if
      End if

*                                                                      *
************************************************************************
*                                                                      *
*     Call MCLR_banner()
      Call CWTime(TCpu1,TWall1)
*                                                                      *
************************************************************************
*                                                                      *
      iAllo=0
c      idp=rtoi
      nrec=get_MBl_wa()/rtob
*
      Call DecideOnCholesky(DoCholesky)
      Call get_iScalar('nSym',nSymX)

      If (DoCholesky .and. nSymX.gt.1) then
       write(6,*)'** Cholesky or RI/DF not implemented with symmetry **'
       Call Quit(_RC_INPUT_ERROR_)
      EndIf

*
      Call Init_Data()
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
      nacpar=(nnA+1)*nnA/2
      nacpr2=(nacpar+1)*nacpar/2

      Call Start_MCLR()
*
      nisp=max(8,nDisp)
*                                                                      *
************************************************************************
*                                                                      *
      Call CWTime(TCpu2,TWall2)
*                                                                      *
************************************************************************
*                                                                      *
*     File pointers
*
      Call mma_allocate(ifpK,nisp,Label='ifpK')
      ifpK(:)=-1
      Call mma_allocate(ifpS,nisp,Label='ifpS')
      ifpS(:)=-1
      Call mma_allocate(ifpRHS,nisp,'ifpRHS')
      ifpRHS(:)=-1
      If (iMethod.eq.2) Then
         Call mma_allocate(ifpCI,nisp,Label='ifpCI')
         Call mma_allocate(ifpSC,nisp,Label='ifpSC')
         Call mma_allocate(ifpRHSCI,nisp,Label='ifpRHSCI')
      Else
         Call mma_allocate(ifpCI,   1,Label='ifpCI')
         Call mma_allocate(ifpSC,   1,Label='ifpSC')
         Call mma_allocate(ifpRHSCI,   1,Label='ifpRHSCI')
      End If
      ifpCI(:)=-1
      ifpSC(:)=-1
      ifpRHSCI(:)=-1
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
*     Calculate response
*
*     Output is stored on disk
*
      If (SPINPOL) Then
         Call WfCtl_SP(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI)
*     Else if (ELECHESS) Then
*        Call WfCtl_PCG(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI)
*        Call Abend()
      Else if(iMCPD) Then!pdft

        Do_Hybrid=.false.
        CALL qpg_DScalar('R_WF_HMC',Do_Hybrid)
        If(Do_Hybrid) Then
         CALL Get_DScalar('R_WF_HMC',WF_Ratio)
         PDFT_Ratio=1.0d0-WF_Ratio
        End If

        if(iMSPD) then
          if(Do_Hybrid) then
           CALL WarningMessage(2,
     &     'Hybrid MS-PDFT gradient not supported yet')
           CALL Quit(_RC_EXIT_EXPECTED_)
          end if
          Call WfCtl_MSPD(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,converged,iPL)
        else
          Call WfCtl_PDFT(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,converged,iPL)
        end if
      Else if (SA.or.PT2) Then
         Call WfCtl_SA(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,converged,iPL)
      Else If (TimeDep) Then
         Call WfCtl_td(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI,
     &                 converged)
      Else
         Call WfCtl_Hess(ifpK,ifpS,ifpCI,ifpSC,ifpRHS,ifpRHSCI,
     &                   converged)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Contract response to hessian etc
*
      If(.not.(TwoStep.and.(StepType.eq.'RUN1'))) Then
         If (PT2.or.SA.or.iMCPD) Then
            Call Out_PT2(ifpK,ifpCI)
            If (PT2) Close (LUPT2) !! this file is opend in wfctl_sa
         Else If (TimeDep) Then
            Call Output_td(ifpK,ifpS,ifpCI,ifpSC,
     &                     ifpRHS,ifpRHSCI,converged)
         Else
            Call Output_mclr(ifpK,ifpS,ifpCI,ifpSC,
     &                       ifpRHS,ifpRHSCI,converged)
            If (mckinley) Call isoloop(Double)
         End If
*
         If (RASSI)   Call OutRAS   (ifpK,ifpCI)
*
         If (TimeDep) Call OutRAS_td(ifpK,ifpCI)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call Basis_Info_Free()
      Call Center_Info_Free()
      Call Symmetry_Info_Free()
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate memory
*
*     Arrays in csf.f
      If (iMethod.eq.2) Then
         Call mma_deallocate(DTOC)
         Call mma_deallocate(CFTP)
         Call mma_deallocate(DFTP)
      End if
      Do i = 1, MXCNSM
         If (Allocated(CNSM(i)%ICONF))
     &                Call mma_deallocate(CNSM(i)%ICONF)
         If (Allocated(CNSM(i)%ICTS))
     &                Call mma_deallocate(CNSM(i)%ICTS)
      End Do

*     Free arrays allocated by memstr.f
      If (iMethod.eq.2) Call FreeStr()
*     Arrays in inpone.f
      Call mma_deallocate(INT1)
*     Arrays in detctl.f
      If (iMethod.eq.2) Then
         Call mma_deallocate(pINT2)
         Call mma_deallocate(pINT1)
      End if

*     Array in rdab.f
      If (iMethod.eq.1) Then
         Call mma_deallocate(CMO)
      End If
*     Arrays in rdjobip_td.f and rdjobiph.f
      If (iMethod.eq.2) Then
         Call mma_deallocate(G2t)
         If (Timedep) Call mma_deallocate(G2sq)
         Call mma_deallocate(G1t)
         Call mma_deallocate(CMO)
      End If

*     Arrays allocated in stpert.f
      Call mma_deallocate(Hss)
      If (SPINPOL) Then
         Call mma_deallocate(FAMO_SpinP)
         Call mma_deallocate(FAMO_SpinM)
         Call mma_deallocate(G2mm)
         Call mma_deallocate(G2mp)
         Call mma_deallocate(G2pp)
         Call mma_deallocate(Fp)
         Call mma_deallocate(Fm)
         Call mma_deallocate(G1p)
         Call mma_deallocate(G1m)
         Call mma_deallocate(SFock)
      End If
*     Arrays allocated in fckmat.f
      Call mma_deallocate(FAMO)
      Call mma_deallocate(Int2)
      Call mma_deallocate(FIMO)
      Call mma_deallocate(F0SQMO)

      Call mma_deallocate(ifpRHSCI)
      Call mma_deallocate(ifpSC)
      Call mma_deallocate(ifpCI)
      Call mma_deallocate(ifpRHS)
      Call mma_deallocate(ifpS)
      Call mma_deallocate(ifpK)
      If (Allocated(SS)) Call mma_deallocate(SS)
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
        Call mma_deallocate(CMO_Inv)
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
      Return
      End
