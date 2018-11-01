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
      Subroutine SuperMac()
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Character*8 Method
      Character*16 StdIn
      Logical Do_Cholesky, Do_ESPF, Numerical, Found
#include "warnings.fh"
#include "temperatures.fh"
*
      Call Get_cArray('Relax Method',Method,8)
*
      Numerical = Method .eq. 'RASSCFSA'    .or.
     &            Method .eq. 'CASPT2'      .or.
     &            Method .eq. 'UHF-SCF'     .or.
     &            Method .eq. 'MBPT2'       .or.
     &            Method .eq. 'CCSDT'       .or.
     &            Method .eq. 'CASSCFSA'    .or.
     &            Method .eq. 'KS-DFT'
*
      If (Method.eq.'CASSCF') Then
         Call Get_iScalar('NumGradRoot',irlxroot)
         Numerical=irlxroot.ne.1
      End If
*
      Call DecideOnCholesky(Do_Cholesky)
      If (Do_Cholesky) Numerical=.true.
*
      Call Qpg_dArray('GeoPC',Found,nData)
      Numerical = Numerical .or. (Found.and.nData.gt.0)
      Call DecideOnESPF(Do_ESPF)
      Numerical = Numerical .or. Do_ESPF
*
*     Analytical PCM frequencies are currently not implemented
*
      Call Qpg_dArray('PCM info',Found,nData)
      Numerical = Numerical .or. (Found.and.nData.gt.0)
*
      Call Qpg_iScalar('DNG',Found)
      If (Found) Then
         Call Get_iScalar('DNG',iDNG)
         Numerical = Numerical .or. (iDNG.eq.1)
      End If
*
      If (.Not.Numerical) Return
*                                                                      *
************************************************************************
*                                                                      *
* Create a backup runfile before running the numerical differentiation
*
      Call fCopy('RUNFILE','RUNBACK',iErr)
      If (iErr.ne.0) Call Abend()
*
      Call GetMem('Scr1','Allo','Inte',ipScr1,7)
      iWork(ipScr1)=0
      iWork(ipScr1+1)=0
      iWork(ipScr1+2)=-99
      iWork(ipScr1+3)=0
      Call Put_iArray('Slapaf Info 1',iWork(ipScr1),7)
      Call GetMem('Scr1','Free','Inte',ipScr1,7)
      LuInput=11
      LuInput=IsFreeUnit(LuInput)
      Call StdIn_Name(StdIn)
      Call Molcas_open(LuInput,StdIn)
*                                                                      *
************************************************************************
*                                                                      *
      iPL = iPrintLevel(-1)
*
      Write (LuInput,'(A)') '>ECHO OFF'
      Write (LuInput,'(A)') '>export MCK_OLD_TRAP=$MOLCAS_TRAP'
      Write (LuInput,'(A)') '>export MCK_OLD_MAXITER=$MOLCAS_MAXITER'
      Write (LuInput,'(A)') '> export MOLCAS_TRAP=ON'
      Write (LuInput,'(A)') '> export MOLCAS_MAXITER=500'
*
*     If SA-CASSCF run the MCLR code so that the reference dipole moment
*     is variational.
*
      If (Method .eq. 'RASSCFSA'.or.Method .eq. 'CASSCFSA')
     &   Write (LuInput,'(A)') '&MCLR'
*
      Write (LuInput,'(A)') '> DO WHILE <'
      Write (LuInput,'(A)') '> IF (ITER NE 1) <'
*
      Call Lu2Lu('SEWARINP',LuInput)
      Write (LuInput,*)
*
      If (Method .eq. 'RASSCFSA'.or.Method .eq. 'CASSCFSA' .or.
     &    Method .eq. 'CASSCF') Then
         Call Lu2Lu('RASSCINP',LuInput)
      Else If (Method .eq. 'CASPT2') Then
         Call Lu2Lu('RASSCINP',LuInput)
         Write (LuInput,'(A)')
         Call Lu2Lu('CASPTINP',LuInput)
      Else If (Method .eq. 'MBPT2') Then
         Call Lu2Lu('SCFINP',LuInput)
      Else If (Method .eq. 'CCSDT') Then
         Call Lu2Lu('SCFINP',LuInput)
         Write (LuInput,'(A)')
         Call Lu2Lu('CCSDTINP',LuInput)
      Else If (Method .eq. 'KS-DFT' .or. Method.eq.'RHF-SCF') Then
         Call Lu2Lu('SCFINP',LuInput)
      End If
*
      Write (LuInput,'(A)') '> END IF <'
*
*     To make sure MBPT2 is run with the Grdt option,
*     run it always (outside the the IF)
*
      If (Method .eq. 'MBPT2') Then
         Write (LuInput,'(A)')
         Call Lu2Lu('MBPT2INP',LuInput)
      End If
*
      Write (LuInput,'(A)')
      Write (LuInput,'(A)') '&Slapaf &End'
      Write (LuInput,'(A)') 'Numerical'
      Write (LuInput,'(A)') 'Iterations'
      Write (LuInput,'(A)') '0'
      Write (LUInput,'(A)') 'THERmochemistry'
*
      Call Get_iScalar('Rotational Symmetry Number',iSigma)
      Write (LUInput,'(I3)') iSigma
      Write (LUInput,'(A)') '1.0'
      Do i=1,NDefTemp
        Write (LUInput,'(F7.2)') DefTemp(i)
      End Do
      Write (LUInput,'(A)') 'End of PT'
      Write (LuInput,'(A)') 'End of Input'
      Write (LuInput,'(A)') '> END DO <'
      Write (LuInput,'(A)') '> export MOLCAS_TRAP=$MCK_OLD_TRAP'
      Write (LuInput,'(A)') '> export MOLCAS_MAXITER=$MCK_OLD_MAXITER'
      Write (LuInput,'(A)') '>ECHO ON'
      Close(LuInput)
*                                                                      *
************************************************************************
*                                                                      *
      Call Finish(_RC_INVOKED_OTHER_MODULE_)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
      Subroutine Lu2Lu(Filename,LuInput)
      Character FileName*(*), Line*180
      Logical Exist
#include "warnings.fh"
*
      Call f_inquire(Filename,Exist)
      If (.Not.Exist) Then
         Write (6,*) 'SuperMac: Missing ',Filename
         Call Finish(_RC_INTERNAL_ERROR_)
      End If
*
      LuSpool2 = 77
      LuSpool2 = IsFreeUnit(LuSpool2)
      Call Molcas_Open(LuSpool2, Filename)
*
 100  Continue
         Read (LuSpool2,'(A)',End=900) Line
         Write(LuInput,'(A)') Line
         Go To 100
 900  Continue
*
      Close(LuSpool2)
*
      Return
      End
