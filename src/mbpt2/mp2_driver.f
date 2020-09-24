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
* Copyright (C) 1998, Roland Lindh                                     *
*               2004,2012, Thomas Bondo Pedersen                       *
*               2007, Francesco Aquilante                              *
************************************************************************
      subroutine MP2_Driver(ireturn)
************************************************************************
*                                                                      *
*     New MP2 driver for Molcas                                        *
*                                                                      *
*                                                                      *
*    Author: R. Lindh                                                  *
*            Dept. of Chemical Physics                                 *
*            University of Lund, Sweden                                *
*            November 7, 1998                                          *
*                                                                      *
*    Modified:                                                         *
*                                                                      *
*       - code using Cholesky vectors directly                         *
*         October 2004, T. B. Pedersen                                 *
*         Dept. of Theoretical Chemistry                               *
*         University of Lund, Sweden                                   *
*                                                                      *
*       - code for the "Scaled Opposite-Spin" (SOS) MP2                *
*         May 2007, F. Aquilante                                       *
*         Dept. of Theoretical Chemistry                               *
*         University of Lund, Sweden                                   *
*                                                                      *
*       - code for Laplace-SOS-MP2 for Cholesky/DF and LDF             *
*         November-December 2012, T. B. Pedersen                       *
*         Centre for Theoretical and Computational Chemistry           *
*         Dept. of Chemistry                                           *
*         University of Oslo, Norway                                   *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "Molcas.fh"
      External Cho_X_GetTol, Seconds
#include "cddos.fh"
#include "trafo.fh"
#include "mxdim.fh"
#include "corbinf.fh"
#include "orbinf2.fh"
#include "mbpt2aux.fh"
#include "files_mbpt2.fh"
#include "timtra.fh"

#include "real.fh"
#include "WrkSpc.fh"
#include "print_mbpt2.fh"
#include "chomp2_cfg.fh"
#include "mp2grad.fh"
#include "namact.fh"
      Common /ctimt/cpubas
      Real*8 E2BJAI, ESCF, REFC, Seconds
      Integer nIsh(8), nAsh(8), nFro_tra(8), nDel_tra(8)
      Logical Ready, Direct, Debug, Exist
      Logical Conventional
      Data Debug/.False./
      Character*8 Method, Method1
      Integer  Cho_X_GetTol
*                                                                      *
************************************************************************
*                                                                      *
      cpubas=seconds()
      Call Set_Data
*
************************************************
*     Check so it is a RHF-SCF reference that is
*     being used.
* TBP, November 2012: do not quit, just issue a warning!
      Call Get_cArray('Relax Method',Method1,8)
      If((Method1(1:7) .ne. 'RHF-SCF').and.
     &   (Method1(1:5) .ne. 'MBPT2')) Then
         Write(6,*)
         Call WarningMessage(1,
     &            'MP2 implementation intended for RHF references only')
         Write(6,'(A,A)')
     &   'MBPT2 WARNING: Reference function according to RunFile:',
     &   Method1
         Write(6,'(A)')
     &   'I''ll assume you know what you''re doing and continue'
         Write(6,*)
         Call xFlush(6)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Figure out if it's a cholesky run, DF run, LDF run
      DoCholesky=.false.
      DoDF=.false.
      DoLDF=.false.
      Call DecideOnCholesky(DoCholesky)
      Call DecideOnDF(DoDF)
      Call DecideOnLocalDF(DoLDF)
*                                                                      *
************************************************************************
*                                                                      *
*     read COMFILE Interface file from SCF and allocate memory...
*
      CALL TIMING(TCPE(1),TCPT,TIOE(1),TIOT)
      Call RdMBPT(mAdCMO,lthCMO,mAdEOr,lthEOr)
*
      Call GetMem('EOcc  ','Allo','Real',mAdEOc,lthEOr)
      Call GetMem('EExt  ','Allo','Real',mAdEEx,lthEOr)
      Call FZero(Work(mAdEOc),lthEOr)
      Call FZero(Work(mAdEEx),lthEOr)
*                                                                      *
************************************************************************
*                                                                      *
*     read Program Input Cards...
*
      Call RdInp(Work(mAdCMO),Work(mAdEOr),Work(mAdEOc),Work(mAdEEx),
     &           iTst,ESCF)
*                                                                      *
************************************************************************
*                                                                      *
      mAdCMO_t=mAdCMO
      mAdEOr_t=mAdEOr
      Call DelGHOST_MBPT(mAdCMO,mAdCMO_t,lthCMO,mAdEOr,mAdEOr_t,lthEOr)
*                                                                      *
************************************************************************
*                                                                      *
      E2BJAI=Zero
      REFC=Zero
      PE2=Zero
      PVECL2=Zero
*
      Wref=Zero
*                                                                      *
************************************************************************
*                                                                      *
*     Write out input parameters
*
      Call PrInp_MBPT2(Work(mAdEOc),Work(mAdEEx))
*                                                                      *
************************************************************************
*                                                                      *
*     Copy pointers to orbital energies to chomp2_dec.fh
*     Needed for amplitude Cholesky decomposition.
*
      If (DoCholesky) Then
         Call ChoMP2_SetPtsOen(mAdEOc,mAdEEx)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Check, if there is an ORDINT file...
*
      Call f_Inquire(FNINTA,Exist)
      Call DecideOnDirect(.True.,Exist,Direct,DoCholesky)
      CALL TIMING(TCPE(2),TCPT,TIOE(2),TIOT)
      CALL TIMING(TCPE(3),TCPT,TIOE(3),TIOT)
*
*
      If (DoT1amp) Then
        l_T1=nOcc(1)*nExt(1)
        nOccT=nOcc(1)
        Do iSym=2,nSym
           l_T1=l_T1+nOcc(iSym)*nExt(iSym)
           nOccT=nOccT+nOcc(iSym)
        End Do
        Call GetMem('T1amp','Allo','Real',jT1amp,l_T1)
        Call Thouless_T1(Work(mAdCMO),nSym,nBas,nFro,nOcc,nExt,
     &                   Work(jT1amp))
        t1nrm=ddot_(l_T1,Work(jT1amp),1,Work(jT1amp),1)
        t1dg=sqrt(t1nrm/nOccT)
        write(6,'(A,F8.4)') '       T1 diagnostic : ',t1dg
        write(6,*)
        Call Set_iOff(nSym,nOcc,nExt,jT1amp-1,iOffT1)
      End If
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
      If (DoLDF) Then ! LDF
         Conventional=.false.
         Ready=.false.
         If (Laplace .and. SOS_MP2) Then
            Call WarningMessage(2,
     &                       'LDF-Laplace-SOS-MP2 not implemented yet!')
            Call SysHalt('mp2_driver')
         Else
            Call WarningMessage(2,
     &                          'Only LDF-Laplace-SOS-MP2 implemented!')
            Call SysHalt('mp2_driver')
         End If
      Else If (DoCholesky .and. ChoAlg.gt.0 .and. .not.SOS_mp2
     &                                      .and. .not.FNOMP2
     &                                      .and. .not.LovMP2) Then
         Conventional = .false.
         Ready = .false.
         Call ChoMP2_Drv(irc,E2BJAI,Work(mAdCMO),
     &                   Work(mAdEOc),Work(mAdEEx))
         If (irc .ne. 0) Then
            Write(6,*) 'MP2_Driver: ChoMP2_Drv returned ',irc
            Call SysAbendMsg('MP2_Driver',
     &                       'Non-zero return code from ChoMP2_Drv',
     &                       ' ')
         Else
            Ready = .true.
         End If
      Else If (DoCholesky .and. SOS_mp2) Then ! CD/DF SOS-MP2
         Conventional=.false.
         Ready=.false.
         If (Laplace) Then
            Call ChoMP2_Drv(irc,E2BJAI,Work(mAdCMO),
     &                      Work(mAdEOc),Work(mAdEEx))
            If (irc .ne. 0) Then
               Write(6,*) 'MP2_Driver: ChoMP2_Drv returned ',irc
               Call SysAbendMsg('MP2_Driver',
     &                          'Non-zero return code from ChoMP2_Drv',
     &                          ' ')
            Else
               Ready=.true.
            End If
         Else
            Call Cho_SOSmp2_Drv(irc,E2BJAI,Work(mAdCMO),
     &                          Work(mAdEOc),Work(mAdEEx))
            If (irc .ne. 0) Then
               Write(6,*) 'SOS-MP2_Driver: Cho_SOSmp2_Drv returned ',irc
               Call SysAbendMsg('SOS-MP2_Driver',
     &                       'Non-zero return code from Cho_SOSmp2_Drv',
     &                          ' ')
            Else
               Ready=.true.
            End If
         End If
!     CD/DF Frozen Natural Orbital MP2
      Else If (DoCholesky .and. FNOMP2) Then
         Conventional = .false.
         Ready = .false.
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Write(6,'(A)') ' Start FNO-MP2 section '
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Write(6,'(A,8I4)')
         Write(6,'(A,I3,A)')
     &  ' NOs specified: ',int(vkept*100),'% of the total virtual space'
         Call FNOMP2_Drv(irc,E2BJAI,Work(mAdCMO),
     &                   Work(mAdEOc),Work(mAdEEx))
         If (irc .ne. 0) Then
            Write(6,*) 'MP2 driver: FNOMP2_Drv returned ',irc
            Call SysAbendMsg('MP2 driver',
     &                       'Non-zero return code from FNOMP2_Drv',
     &                       ' ')
         Else
            Ready = .true.
         End If
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Write(6,'(A)') ' End FNO-MP2 section '
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Write(6,'(A,8I4)')
         Write(6,'(A,8I4)')
      Else If (DoCholesky .and. LovMP2) Then ! CD/DF Localized O-V MP2
         Conventional = .false.
         Ready = .false.
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Write(6,'(A)') ' Start LovMP2 section '
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Write(6,'(A,8I4)')
         Call LovMP2_Drv(irc,E2BJAI,Work(mAdCMO),
     &                   Work(mAdEOc),Work(mAdEEx),NamAct,nActa,
     &                   ThrLov,DoMP2,all_Vir)
         If (irc .ne. 0) Then
            Write(6,*) 'MP2 driver: LovMP2_Drv returned ',irc
            Call SysAbendMsg('MP2 driver',
     &                       'Non-zero return code from LovMP2_Drv',
     &                       ' ')
         Else
            Ready = .true.
         End If
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Write(6,'(A)') ' End LovMP2 section '
         Write(6,'(A)')
     &   '-------------------------------------------------------'
         Write(6,'(A,8I4)')
         Write(6,'(A,8I4)')
      Else ! conventional (possibly with Cholesky)
         Conventional = .true.
         If (DoCholesky) Then
            If (ChoAlg .eq. 0 .and. iPL.ge.2) Then
               Write(6,*) 'Conventional algorithm used.'
               Write(6,*)
               Write(6,*) 'Integrals generated from Cholesky vectors '
     &                   //'(Algorithm 0):'
            Else
               Call SysHalt('mp2_driver')
            End If
         Else
            If (iPL.ge.2) Write(6,*) 'Conventional algorithm used...'
         End If
         If ( iTst.ne.0 ) Go To 100
*
*------- Use the driver from the CASPT2 code.
*
         Call DaName_MF_wa(LuIntM,FnIntM)

         If(DoDens) Then

*--------   Same as just above but keeping frozen and deleted orbitals in the
*--------   transformation
            Do i=1,8
               nIsh(i) = nOrb(i)+nDel(i)
            EndDo
            itmp=0
            Call ICopy(8,[itmp],0,nFro_tra,1)
            Call ICopy(8,[itmp],0,nDel_tra,1)
         Else
            Call ICopy(8,nOcc,1,nIsh,1)
         End If

         itmp=0
         Call ICopy(8,[itmp],0,nAsh,1)
         If(.not.DoDens) Then
* PAM Jan 2013: Set correct nOrb:
      do i=1,8
        norb(i)=norb(i)-nfro(i)
      end do
            Call SetUp_CASPT2_Tra(nSym,nBas,nOrb,nIsh,nAsh,
     &                            nFro,nDel,mAdCMO,lthCMO,
     &                            LuIntM,LuHlf1,LuHlf2,LuHlf3)
* End of patch
         Else
            Call SetUp_CASPT2_Tra(nSym,nBas,nIsh,nIsh,nAsh,
     &                            nFro_tra,nDel_tra,mAdCMO,lthCMO,
     &                            LuIntM,LuHlf1,LuHlf2,LuHlf3)
         End If
         If(.NOT.DoCholesky) then
           iRC=-1
           iOpt=0
           Call OpnOrd(iRC,iOpt,FnIntA,LuIntA)
           If (iRC.ne.0) Then
             Write (6,*) 'mp2_driver: error opening MOLINT'
             Call Abend
           End If
         EndIf
*
CGG         Call TraCtl
         iType=1  ! Means that TraCtl_Drv is called by MP2
CGG      DoExch2=.True. ! Do generate Exch-2 integrals (not really used).
         Call TraCtl_Drv(iType,.True.,8)
*
         Call DaClos(LuHlf1)
         Call DaClos(LuHlf2)
         Call DaClos(LuHlf3)
*
*
         CALL TIMING(TCPE(3),TCPT,TIOE(3),TIOT)
*        PRINT TRANSFORMED INTEGRALS (USED ONLY FOR DEBUGGING PURPOSES)
         If (Debug) Then
            iPrc=0
            Call RDInt2_MP2(iPrc)
         ENDIF
*                                                                      *
************************************************************************
*                                                                      *
         If(DoDens) then
*
*           Obtain variational density if requested by user
*
            ipEOcc=mAdEOc
            ipEVir=mAdEEx
            ipCMO=mAdCMO
            Call Mp2Dens_drv(E2BJAI,REFC)
         Else
*
*           COMPUTE CORRELATION ENERGY CONTRIBUTION
*
            CALL BJAI(IADOUT,Work(mAdEOc),Work(mAdEEx),E2BJAI,REFC)
         End If
         CALL TIMING(TCPE(4),TCPT,TIOE(4),TIOT)

         Ready=.True.

*        Close files (for Cholesky, OrdInt was never opened)
         If (.not. DoCholesky) Then
            iRc=-1
            iOpt=0
            Call ClsOrd(iRc,iOpt)
            If (iRc.ne.0) Then
               Write (6,*) 'MP2_Driver: Error closing ORDINT'
               Call Abend
            End If
         End If
         Call DaClos(LuIntM)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Write out resulting energy, etc.
*
      If (Ready) Then

         Wref=One/(One+Wref)   ! Note: this is Cref**2

         If (DoLDF) Then
            Call WarningMessage(2,'LDF should not be implemented....')
            Call SysHalt('mp2_driver')
         Else If (DoCholesky .and. ChoAlg.gt.0 .and. .not.SOS_mp2
     &                                        .and. .not.LovMP2
     &                                        .and. .not.FNOMP2) Then
            If (iPL.ge.2)
     &      Write (6,'(3(/6X,A,F20.10,A)//6X,A,F20.10,A//6X,A,F15.5)')
     &    ' SCF energy                           =',ESCF,       ' a.u.',
     &    ' Second-order correlation energy      =',E2BJAI,     ' a.u.',
     &    ' ( Opposite-Spin contribution         =',-EOSMP2,    ' )',
     &    ' Total energy                         =',E2BJAI+ESCF,' a.u.',
     &    ' Reference weight ( Cref**2 )         =',Wref
         ElseIf (DoCholesky .and. LovMP2) Then
            ESSMP2=E2BJAI+EOSMP2
            If (iPL.ge.2)
     &      Write (6,'(4(/6X,A,F20.10,A)//6X,A,F20.10,A//6X,A,F15.5)')
     &    ' SCF energy                           =',ESCF,       ' a.u.',
     &    ' Second-order correlation energy      =',E2BJAI,     ' a.u.',
     &    ' ( Opposite-Spin contribution         =',-EOSMP2,    ' )',
     &    ' ( Same-Spin contribution             =', ESSMP2,    ' )',
     &    ' Total energy                         =',E2BJAI+ESCF,' a.u.',
     &    ' Reference weight ( Cref**2 )         =',Wref
         ElseIf (DoCholesky .and. FNOMP2) Then
            ESSMP2=E2BJAI+EOSMP2-XEMP2
            If (iPL.ge.2)
     &      Write (6,'(5(/6X,A,F20.10,A)//6X,A,F20.10,A//6X,A,F15.5)')
     &    ' SCF energy                           =',ESCF,       ' a.u.',
     &    ' Second-order correlation energy      =',E2BJAI,     ' a.u.',
     &    ' ( Opposite-Spin contribution         =',-EOSMP2,    ' )',
     &    ' ( Same-Spin contribution             =', ESSMP2,    ' )',
     &    ' ( Truncation error estimate          =', XEMP2,     ' )',
     &    ' Total energy                         =',E2BJAI+ESCF,' a.u.',
     &    ' Reference weight ( Cref**2 )         =',Wref
         ElseIf (DoCholesky .and. SOS_mp2) Then ! CD/DF SOS-MP2 energy
            E2BJAI = C_os*E2BJAI
            If (iPL.ge.2) Then
               If (Laplace) Then
                  Write(6,'(/,6X,A,I4)')
     &    ' Number of Laplace grid points:',Laplace_nGridPoints
                  Write (6,'(3(/6X,A,F20.10,A)//6X,A,F20.10,A)')
     &    ' Opposite-Spin (OS) scaling factor    =',C_os,       '     ',
     &    ' SCF energy                           =',ESCF,       ' a.u.',
     &    ' L-SOS 2nd-order correlation energy   =',E2BJAI,     ' a.u.',
     &    ' Total L-SOS-MP2 energy               =',E2BJAI+ESCF,' a.u.'
               Else
                  Write (6,'(3(/6X,A,F20.10,A)//6X,A,F20.10,A)')
     &    ' Opposite-Spin (OS) scaling factor    =',C_os,       '     ',
     &    ' SCF energy                           =',ESCF,       ' a.u.',
     &    ' SOS 2nd-order correlation energy     =',E2BJAI,     ' a.u.',
     &    ' Total SOS-MP2 energy                 =',E2BJAI+ESCF,' a.u.'
               End If
            End If
         Else
            WRef=REFC**2
            If (iPL.ge.2)
     &      Write (6,'(2(/6X,A,F20.10,A)//6X,A,F20.10,A/6X,A,F15.5)')
     &    ' SCF energy                           =',ESCF,       ' a.u.',
     &    ' Second-order correlation energy      =',E2BJAI,     ' a.u.',
     &    ' Total energy                         =',E2BJAI+ESCF,' a.u.',
     &    ' Reference weight ( Cref**2 )         =',Wref
         End If
         If (iPL.ge.2) Then
           ETot=E2BJAI+ESCF
           Write(6,*)
           Call PrintResult(6,'(6X,A,T50,F19.10)',
     &                      'Total MBPT2 energy',0,'',[ETot],1)
         End If
         Call xFlush(6)
         If (iPL.ge.2) Write (6,*)
         Call Store_Energies(1,E2BJAI+ESCF,1)
         Method='MBPT2   '
         Call Put_cArray('Relax Method',Method,8)
         If(DoDens) Call Prpt()
         If (iPL.ge.2) Write (6,*)
      Else
         Write (6,*) ' Energy evaluation not completed.'
      End If
*
*     Shanks-type series convergence acceleration
*
      Call Compute_Shanks(ESCF,E2BJAI+ESCF,Work(mAdEOr),lthEOr,
     &                    nBas,nFro,nOcc,nSym,E0,Shanks1_E)
*
      If (iPL.ge.2) Then
         WRITE(6,'(6X,A,A,(F20.10),A)')
     &    ' Zeroth-order energy (E0)   ',
     &    '          =',E0,' a.u.'
         write(6,*)

         WRITE(6,'(6X,A,A,(F20.10),A)')
     &    ' Shanks-type energy S1(E)   ',
     &    '          =',Shanks1_E,' a.u.'
         write(6,*)
         write(6,*)
      End If
*
*       PRINT PROCESSING AND TIMING INFORMATION
*
      If (Conventional.and.iPL.ge.2) Then
         WRITE (6,'(//6X,A//6X,A/6X,A/)')
     &' Data processing and timing information:',
     &' Section                                              time(sec)',
     &'                                                    CPU  Elapsed'
         WRITE(6,'(6X,A,A,2(2X,F8.2))') 'Input data processing      ',
     &    '                  ',TCPE(2)-TCPE(1),TIOE(2)-TIOE(1)
         WRITE(6,'(6X,A,A,2(2X,F8.2))') 'Transformation of integrals',
     &    '                  ',TCPE(3)-TCPE(2),TIOE(3)-TIOE(2)
         WRITE(6,'(6X,A,A,2(2X,F8.2))') 'MBPT2 calculations (BJAI)  ',
     &    '                  ',TCPE(4)-TCPE(3),TIOE(4)-TIOE(3)
         WRITE(6,'(6X,A,A,2(2X,F8.2))') 'Total MBPT2 calculations   ',
     &    '                  ',TCPE(4)-TCPE(1),TIOE(4)-TIOE(1)
         WRITE(6,*)
         WRITE(6,*)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      iTol = Cho_X_GetTol(8)
      Call Add_Info('E_MP2',[E2BJAI+ESCF],1,iTol)
      Call Add_Info('HF ref weight',[WREF],1,iTol)
*                                                                      *
************************************************************************
*                                                                      *
 100  Continue
*
*---- Close up calculation
*
      If (DoT1amp) Call GetMem('T1amp','Free','Real',jT1amp,l_T1)
      Call GetMem('EExt  ','Free','Real',mAdEEx,lthEOr)
      Call GetMem('EOcc  ','Free','Real',mAdEOc,lthEOr)
      Call GetMem('EOrb  ','Free','Real',mAdEOr,lthEOr)
      Call GetMem('CMO   ','Free','Real',mAdCMO,lthCMO)
*
      Call qStat(' ')
      ireturn=0
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
