************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Alaska_Super_Driver(iRC)
      Implicit Real*8 (a-h,o-z)
      Character*8 Method
      Logical Do_Cholesky, Numerical, Do_DF, Do_ESPF, StandAlone, Exist,
     &        Do_Numerical_Cholesky, Do_1CCD, MCLR_Ready, Reduce_Prt
      External Reduce_Prt
      Character FileName*128, Line*180, KSDFT*16
      Character*16 StdIn
      Integer Columbus
      Character(Len=16) mstate1,mstate2
#include "warnings.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "nac.fh"
#include "alaska_root.fh"
#include "para_info.fh"
*                                                                      *
************************************************************************
*                                                                      *
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
*                                                                      *
************************************************************************
*                                                                      *
*
*     Check the input for the numerical keyword
*
*
*     copy STDINP to LuSpool
*
      LuSpool=37
      Call SpoolInp(LuSpool)
      Call Chk_Numerical(LuSpool,Numerical)
*
*     Call the numerical procedure if numerical option is available.
*     Otherwise hope that the analytic code know how to handle the
*     case.
*
      Call Get_cArray('Relax Method',Method,8)
      Call Get_iScalar('Columbus',Columbus)
*                                                                      *
************************************************************************
*                                                                      *
*     Default for Cholesky or RI/DF is to do the numerical procedure.
*     However, for pure DFT we have analytic gradients.
*
      Call DecideOnCholesky(Do_Cholesky)
      Call DecideOnDF      (Do_DF)
      Call DecideOn1CCD    (Do_1CCD)
      Call Get_iScalar('NSYM',nSym)
*
*     Default for MBPT2 is the numerical procedure but if variational
*     densities are calculated analytical gradients shall be used.
      iMp2Prpt = 0
      If(Method .eq. 'MBPT2   ') Then
         Call Get_iScalar('mp2prpt',iMp2Prpt)
*
*        Make sure that the analytic procedure is used if possible.
*
         If (nSym.eq.1 .and. iMp2Prpt.ne.2 .and..NOT.Numerical) Then
            Call WarningMessage(2,'Error in Alaska_Super_Driver')
            Write (6,*) 'Alaska: the MBPT2 module was run without the'
     &                      //' Grdt option!'
            Write (6,*) '   Correct the input and restart the'
     &                      //' calculation!'
            Call Abend()
         End If
      End If
*
      Do_Numerical_Cholesky = Do_Cholesky .or. Do_DF
      Call Get_iScalar('agrad',iForceAnalytical)
      If(iForceAnalytical .eq. 1) Do_Numerical_Cholesky=.False.
*
      ExFac=0.0D0
      If (Method.eq.'KS-DFT  '.and.Do_Numerical_Cholesky) Then
         Call Get_cArray('DFT functional',KSDFT,16)
         ExFac=Get_ExFac(KSDFT)
*
         If (Do_DF                 .or.                  ! RI/DF
     &       (Do_Cholesky.and.Do_1CCD.and.nSym.eq.1)      ! 1C-CD
     &                                ) Then
            Do_Numerical_Cholesky=.False.
         End If

      End If
*
      If (( Do_DF                    .or.
     &     (Do_Cholesky.and.Do_1CCD.and.nSym.eq.1))) Then
*
         If( (Method .eq. 'KS-DFT  ') .or.
     &       (Method .eq. 'UHF-SCF ') .or.
     &       (Method .eq. 'RHF-SCF ') ) Then
            Do_Numerical_Cholesky= .False.
     &
         Else If(Method.eq.'MBPT2   ' .and.nSym.eq.1) Then
            Do_Numerical_Cholesky = .False.
*           Write(6,*) 'Do Numerical', Do_Numerical_Cholesky
         End If
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Call DecideOnESPF(Do_ESPF)
*     No ESPF for NAC, right?
*     If (isNAC) Do_ESPF=.False.
*                                                                      *
************************************************************************
      if(Method .eq. 'DMRGSCFS')then
        Call Get_iScalar('SA ready',iGo)
      end if
*                                                                      *
      If (Numerical              .OR.
     &    Do_Numerical_Cholesky  .OR.
     &    Method .eq. 'RASSCFSA' .OR.
     &    Method .eq. 'GASSCFSA' .OR.
     &  ((Method .eq. 'DMRGSCFS').and.(iGo.ne.2)) .OR.
     &    Method .eq. 'CASPT2'   .OR.
     &  ((Method .eq. 'MBPT2').and.(iMp2Prpt.ne.2)) .OR.
     &    Method .eq. 'CCSDT'    ) Then
         If (isNAC) Then
           Call Store_Not_Grad(0,NACstates(1),NACstates(2))
           Call WarningMessage(2,'Numerical nonadiabatic coupling not'
     &                           //' implemented')
           If (Auto) Then
             iRC=0
             Return
           Else
             Call Abend()
           End If
         EndIf
*                                                                      *
************************************************************************
*                                                                      *
*        Numerical gradients to be used!
*        Alaska will automatically generate the input for CASPT2_Gradient
*        and signal to AUTO (iRC=2) to run the input file Stdin.x.
*
         If (iPL.ge.3) Then
            Write (6,*)
            Write (6,*) ' Alaska requests the Numerical_Gradient'
     &                //' module to be executed!'
            Write (6,*)
         End If
*
         LuInput=11
         LuInput=IsFreeUnit(LuInput)
         Call StdIn_Name(StdIn)
         Call Molcas_Open(LuInput,StdIn)
*
         Write (LuInput,'(A)') '>ECHO OFF'
         Write (LuInput,'(A)') '>export AL_OLD_TRAP=$MOLCAS_TRAP'
         Write (LuInput,'(A)') '>export MOLCAS_TRAP=ON'
*
         Write (LuInput,'(A)') ' &NUMERICAL_GRADIENT &End'
         Write (LuInput,'(A)') 'End of Input'
         Write (LuInput,'(A)') '>export MOLCAS_TRAP=$AL_OLD_TRAP'
         Write (LuInput,'(A)') '>ECHO ON'
         Close(LuInput)
         Call Finish(_RC_INVOKED_OTHER_MODULE_)
*                                                                      *
************************************************************************
*                                                                      *
*     These do not work in parallel. Warn and stop early, better than
*     crash or give wrong results
*
      Else If (Do_Cholesky.and.(Method.eq.'CASSCFSA')
     &        .and.(nProcs.gt.1)) Then
         Call WarningMessage(2,'Error in Alaska_Super_Driver')
         Write (6,*) 'RI SA-CASSCF analytical gradients do not work'
     &             //' correctly in parallel (yet).'
         Call Abend()
      Else If (Do_Cholesky.and.(Method.eq.'MBPT2')
     &        .and.(nProcs.gt.1)) Then
         Call WarningMessage(2,'Error in Alaska_Super_Driver')
         Write (6,*) 'RI MBPT2 analytical gradients do not work'
     &             //' correctly in parallel (yet).'
         Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
      Else If (Method.eq.'CASSCFSA' .or.
     &        (Method.eq.'DMRGSCFS' .and. iGo.ne.2)) Then
*                                                                      *
************************************************************************
*                                                                      *
*        State-Average CASSCF / DMRGSCF
*
         Call Get_iScalar('SA ready',iGo)
         Call Get_iScalar('Relax CASSCF root',iRlxRoot)

         If (iRlxRoot.eq.0) iRlxRoot=1
         If (isNAC) Then
            Write(mstate1,'(1X,I7,",",I7)') NACStates(1),NACStates(2)
         Else
            Write(mstate1,'(I16)') iRlxRoot
         End If

*        iGo=-1 non-equivalent multi state SA-CASSCF
*        iGo=0  equivalent multi state SA-CASSCF
*        iGo=2  single root SA-CASSCF
         mstate2=''
         if(iGo.ne.2)then
           Call Get_cArray('MCLR Root',mstate2,16)
         end if

* If an explicit root was requested in MCLR and none in ALASKA,
* go for it
         If (DefRoot) Then
            If (mstate2(1:1).eq.'+') Then
               mstate1=mstate2
               If (Index(mstate2,',').ne.0) Then
                  Read(mstate2,'(1X,I7,1X,I7)')
     &               NACStates(1),NACStates(2)
                  ForceNAC=.true.
               Else
                  Read(mstate2,'(1X,I15)') iRlxRoot
               End If
            End If
         End If
         mstate1(1:1)=mstate2(1:1)
         MCLR_Ready=(iGO.eq.1).and.(mstate1.eq.mstate2)

         If (MCLR_Ready.or.(iGO.gt.1)) Then
            Call Alaska(LuSpool,iRC)
*
*        Add ESPF contribution
*
            If (Do_ESPF) Then
               StandAlone=.False.
               Call ESPF(iReturn,StandAlone)
               If (iReturn.ne.0) Then
                  Call WarningMessage(2,'Error in Alaska_Super_Driver')
                  Write (6,*) 'Alaska: ESPF finish with non-zero return'
     &                      //' code!'
                  Call Abend()
               End If
            End If
*           Reset iGO to 0 to allow for new MCLR/ALASKA calculations
            If (iGo.eq.1) iGo=0
            Call Put_iScalar('SA ready',iGo)
         Else If (iGO.eq.-1) Then
            Call WarningMessage(2,'Error in Alaska_Super_Driver')
            Write (6,*) 'Gradients not implemented for SA-CASSCF'//
     &                  ' with non-equivalent weights!'
            Call Abend()
         Else
            If (iPL.ge.3) Then
               Write (6,*)
               Write (6,*)
     &           ' Alaska requests MCLR to be run before it starts'
     &           //' again!'
               Write (6,*)
            End If
*
            LuInput=11
            LuInput=IsFreeUnit(LuInput)
            Call StdIn_Name(StdIn)
            Call Molcas_open(LuInput,StdIn)
*
            Write (LuInput,'(A)') '>ECHO OFF'
            Write (LuInput,'(A)') '>export AL_OLD_TRAP=$MOLCAS_TRAP'
            Write (LuInput,'(A)') '>export MOLCAS_TRAP=ON'
*
            Write (LuInput,'(A)') ' &MCLR &End'
            If (isNAC) Then
              Write (LuInput,'(A)') 'NAC'
              Write (LuInput,'(I5,1X,I5)') NACstates(1),NACstates(2)
            EndIf
            Write (LuInput,'(A)') 'End of Input'
            Write (LuInput,'(A)') ' '
*
            FileName='ALASKINP'
            Call f_inquire(Filename,Exist)
*
            If (Exist) Then
               LuSpool2 = 77
               LuSpool2 = IsFreeUnit(LuSpool2)
               Call Molcas_Open(LuSpool2, Filename)
*
 100           Continue
               Read (LuSpool2,'(A)',End=900) Line
               Write(LuInput,'(A)') Line
               Go To 100
 900           Continue
*
               Close(LuSpool2)
*
            Else
*
               Write (LuInput,'(A)') ' &Alaska &End'
*              Write (LuInput,'(A)') 'Show'
               Write (LuInput,'(A)') 'CutOff'
               Write (LuInput,'(A)') '1.0D-7'
               Write (LuInput,'(A)') 'End of Input'
*
            End If
*
            Write (LuInput,'(A)') '>RM -FORCE $Project.MckInt'
            Write (LuInput,'(A)') '>export MOLCAS_TRAP=$AL_OLD_TRAP'
            Write (LuInput,'(A)') '>ECHO ON'
            Close(LuInput)
            Call Finish(_RC_INVOKED_OTHER_MODULE_)
*
           End If

************************************************************************
*                                                                      *
      Else If (Method.eq.'MCPDFT') Then
*                                                                      *
************************************************************************
*                                                                      *
*        MC-PDFT calculation
*
         Do_ESPF = .False.
         Call Get_iScalar('SA ready',iGo)
         Call Get_iScalar('Relax CASSCF root',iRlxRoot)

         !Andrew - I need to identify the root and make sure it is not a
         !state averaged calculation.  iGo=1 means do MCLR

         !iGo=99 means the potentials were not calculated during the
         !MCPDFT step, which is required for analytic gradients.
         If(iGO.eq.99) then
            Call WarningMessage(2,'Error in Alaska_Super_Driver')
            Write (6,*) 'MC-PDFT was run without the GRADient'//
     &                  ' keyword.  Analytic gradients require'//
     &                  ' this keyword.  Please use the GRADient'//
     &                  ' keyword in the preceeding MC-PDFT step.'
            Call Abend()
         End if

         If (iRlxRoot.eq.0) iRlxRoot=1
!         If (isNAC) Then
!            Write(mstate1,'(1X,I7,",",I7)') NACStates(1),NACStates(2)
!         Else
            Write(mstate1,'(I16)') iRlxRoot
!         End If

*        iGo=-1 non-equivalent multi state SA-CASSCF
*        iGo=0  equivalent multi state SA-CASSCF
*        iGo=2  single root SA-CASSCF
         mstate2=''
         if(iGo.ne.2)then
           Call Get_cArray('MCLR Root',mstate2,16)
         end if

* If an explicit root was requested in MCLR and none in ALASKA,
* go for it
         If (DefRoot) Then
            If (mstate2(1:1).eq.'+') Then
               mstate1=mstate2
               If (Index(mstate2,',').ne.0) Then
                  Read(mstate2,'(1X,I7,1X,I7)')
     &               NACStates(1),NACStates(2)
                  ForceNAC=.true.
               Else
                  Read(mstate2,'(1X,I15)') iRlxRoot
               End If
            End If
         End If
         mstate1(1:1)=mstate2(1:1)
         MCLR_Ready=(iGO.eq.1).and.(mstate1.eq.mstate2)

         If (MCLR_Ready.or.(iGO.gt.1)) Then
            Call Alaska(LuSpool,iRC)
*
*        Add ESPF contribution
*
!            If (Do_ESPF) Then
!               StandAlone=.False.
!               Call ESPF(iReturn,StandAlone)
!               If (iReturn.ne.0) Then
!                  Call WarningMessage(2,'Error in Alaska_Super_Driver')
!                  Write (6,*) 'Alaska: ESPF finish with non-zero return'
!     &                      //' code!'
!                  Call Abend()
!               End If
!            End If
*           Reset iGO to 0 to allow for new MCLR/ALASKA calculations
            If (iGo.eq.1) iGo=0
            Call Put_iScalar('SA ready',iGo)
         Else If (iGO.eq.-1) Then
            Call WarningMessage(2,'Error in Alaska_Super_Driver')
            Write (6,*) 'Gradients not implemented for SA-CASSCF'//
     &                  ' with non-equivalent weights!'
            Call Abend()
         Else
            If (iPL.ge.3) Then
               Write (6,*)
               Write (6,*)
     &           ' Alaska requests MCLR to be run before it starts'
     &           //' again!'
               Write (6,*)
            End If
*
            LuInput=11
            LuInput=IsFreeUnit(LuInput)
            Call StdIn_Name(StdIn)
            Call Molcas_open(LuInput,StdIn)
*
            Write (LuInput,'(A)') '>ECHO OFF'
            Write (LuInput,'(A)') '>export AL_OLD_TRAP=$MOLCAS_TRAP'
            Write (LuInput,'(A)') '>export MOLCAS_TRAP=ON'
*
            Write (LuInput,'(A)') ' &MCLR &End'
            Write (LuInput,'(A)') ' PRINT = 100'
!            If (isNAC) Then
!              Write (LuInput,'(A)') 'NAC'
!              Write (LuInput,'(I5,1X,I5)') NACstates(1),NACstates(2)
!            EndIf
            Write (LuInput,'(A)') 'End of Input'
            Write (LuInput,'(A)') ' '
*
            FileName='ALASKINP'
            Call f_inquire(Filename,Exist)
*
            If (Exist) Then
               LuSpool2 = 77
               LuSpool2 = IsFreeUnit(LuSpool2)
               Call Molcas_Open(LuSpool2, Filename)
*
 105           Continue
               Read (LuSpool2,'(A)',End=905) Line
               Write(LuInput,'(A)') Line
               Go To 105
 905           Continue
*
               Close(LuSpool2)
*
            Else
*
               Write (LuInput,'(A)') ' &Alaska &End'
*              Write (LuInput,'(A)') 'Show'
               Write (LuInput,'(A)') 'CutOff'
               Write (LuInput,'(A)') '1.0D-7'
               Write (LuInput,'(A)') 'End of Input'
*
            End If
*
            Write (LuInput,'(A)') '>RM -FORCE $Project.MckInt'
            Write (LuInput,'(A)') '>export MOLCAS_TRAP=$AL_OLD_TRAP'
            Write (LuInput,'(A)') '>ECHO ON'
            Close(LuInput)
            Call Finish(_RC_INVOKED_OTHER_MODULE_)
*
           End If


************ columbus interface ****************************************

       Else If (method.eq.'MR-CISD ' .and. Columbus.eq.1) Then

*
*       COLUMBUS MR-CI gradient:
*           effective density matrix on RUNFILE 'densao_var'
*           effective fock matrix  on RUNFILE   'FockOcc'
*           effective  D2          on GAMMA
*           G_TOC indices          on binary file gtoc
*           RUNFILE read in drvh1
*           GAMMA read in pget0
*           gtoc read in drvg1

*
*
           Call Get_iScalar('Relax CASSCF root',iRlxRoot)
           Call Alaska(LuSpool,iRC)
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
*        Read the root, as it could be a CASSCF state-specific excited
*        state calculation

*
         iRlxRoot=0
         Call Qpg_iScalar('Relax CASSCF root',Exist)
         If (Exist) Call Get_iScalar('Relax CASSCF root',iRlxRoot)
         If (iRlxRoot.eq.0) iRlxRoot=1
*
*        Go ahead and compute the gradients
*
         Call Alaska(LuSpool,iRC)
*
*        Add ESPF contribution
*
         If (Do_ESPF) Then
            StandAlone=.False.
            Call ESPF(iReturn,StandAlone)
            If (iReturn.ne.0) Then
               Call WarningMessage(2,'Error in Alaska_Super_Driver')
               Write (6,*) 'Alaska: ESPF finish with non-zero return'
     &                   //' code!'
               Call Abend()
            End If
         End If

*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
* Read the gradient and store in the GRADS file
* It is done here because the gradient could have been modified by ESPF
* and we do not want to pass the root to ESPF (yet)
*
      Call Get_Grad(ipGrad,nGrad)
      If (isNAC) Then
        Call Store_Grad(Work(ipGrad),nGrad,0,NACstates(1),NACstates(2))
      Else
        Call Store_Grad(Work(ipGrad),nGrad,iRlxRoot,0,0)
      End If
      Call GetMem('Grad','Free','Real',ipGrad,nGrad)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
