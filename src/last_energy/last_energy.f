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
      Subroutine Last_Energy(ireturn)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "timtra.fh"
#include "real.fh"
      Character*8 Method
      Logical Do_ESPF, StandAlone, FoundLastEn
*                                                                      *
************************************************************************
*                                                                      *
      ireturn=99
*
*     Get information regarding the last method used
*
      Call qpg_cArray('LastEnergyMethod',FoundLastEn,lengthlast)
      If (FoundLastEn) Then
         Call Get_cArray('LastEnergyMethod',Method,8)
      Else
         Call Get_cArray('Relax Method',Method,8)
      EndIf
*
      Call DecideOnESPF(Do_ESPF)
*
*
      If (Method(5:7) .eq. 'SCF'    .OR.
     &    Method(1:6) .eq. 'KS-DFT' .OR.
     &    Method(1:6) .eq. 'CASSCF' .OR.
     &    Method(1:6) .eq. 'RASSCF' .OR.
     &    Method(1:6) .eq. 'CASPT2' .OR.
     &    Method(1:5) .eq. 'MBPT2'  .OR.
     &    Method(1:5) .eq. 'CCSDT'  .OR.
     &    Method(1:4) .eq. 'CHCC'   .OR.
     &    Method(1:6) .eq. 'MCPDFT' .OR.
     &    Method(1:4) .eq. 'CHT3') Then
         Continue
      Else
         Write (6,'(A,A,A)') 'Last Energy for ',Method,
     &                    ' is not implemented yet.'
         Call Abend()
      End If
*
*                                                                      *
************************************************************************
*                                                                      *
*     Compute integrals
*
      Call StartLight('seward')
      Call Disable_Spool()
      Call Seward(ireturn)
      If (iReturn .ne. 0) Then
         Write(6,*) 'Last_Energy failed ...'
         Write(6,*)'Seward returned with return code, rc = ',iReturn
         Call Abend()
      End If
*
*     Compute ESPF
*
      If (Do_ESPF) Then
         Call StartLight('espf')
         Call Disable_Spool()
         StandAlone=.True.
         Call ESPF(ireturn,StandAlone)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*)'Espf returned with return code, rc = ',iReturn
            Call Abend()
         End If
      End If
*
*     Compute the wave function
*
      If (Method(5:7) .eq. 'SCF' .OR.
     &    Method(1:6) .eq. 'KS-DFT' .OR.
     &    Method(1:5) .eq. 'MBPT2' . OR.
     &    Method(1:4) .eq. 'CHCC' . OR.
     &    Method(1:4) .eq. 'CHT3') Then
         Call StartLight('scf')
         Call Disable_Spool()
         Call xml_open('module',' ',' ',0,'scf')
         Call SCF(iReturn)
         Call xml_close('module')
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'SCF returned with return code, rc = ',iReturn
            Call Abend()
         End If
      Else If (Method(1:6) .eq. 'RASSCF' .OR.
     &             Method(1:6) .eq. 'CASSCF' .OR.
     &             Method(1:6) .eq. 'CASPT2' .OR.
     &             Method(1:6) .eq. 'MCPDFT' .OR.
     &             Method(1:5) .eq. 'CCSDT') Then
         Call StartLight('rasscf')
         Call Disable_Spool()
         Call RASSCF(ireturn)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'RASSCF returned with return code, rc = ',
     &                  iReturn
            Call Abend()
         End If
      End If

*
      If (Method(1:5) .eq. 'MBPT2') Then
         Call StartLight('mbpt2')
         Call Disable_Spool()
         Call MP2_Driver(ireturn)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'MBPT2 returned with return code, rc = ',
     &                  iReturn
            Call Abend()
         End If
      End If
*
      If (Method(1:5) .eq. 'CCSDT') Then
         Call StartLight('motra')
         Call Disable_Spool()
         Call Motra(ireturn)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'Motra returned with return code, rc = ',
     &                  iReturn
            Call Abend()
         End If
*
         Call StartLight('ccsdt')
         Call Disable_Spool()
         Call CCSDT(ireturn)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'CCSDT returned with return code, rc = ',
     &                  iReturn
            Call Abend()
         End If
      End If
*
      If (Method(1:4) .eq. 'CHCC' .OR.
     &    Method(1:4) .eq. 'CHT3') Then
         Call StartLight('chcc')
         Call Disable_Spool()
         Call CHCC(ireturn)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'CHCC returned with return code, rc = ',
     &                  iReturn
            Call Abend()
         End If
      End If
*
      If (Method(1:4) .eq. 'CHT3') Then
         Call StartLight('cht3')
         Call Disable_Spool()
         Call CHT3(ireturn)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'CHT3 returned with return code, rc = ',
     &                  iReturn
            Call Abend()
         End If
      End If
*
      If (Method(1:6) .eq. 'CASPT2') Then
         Call StartLight('caspt2')
         Call Disable_Spool()
         Call CASPT2(ireturn)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'CASPT2 returned with return code, rc = ',
     &                  iReturn
            Call Abend()
         End If
      End If
*
      If (Method(1:6) .eq. 'MCPDFT') Then
         Call StartLight('mcpdft')
         Call Disable_Spool()
         Call mcpdft(ireturn)
         If (iReturn .ne. 0) Then
            Write(6,*) 'Last_Energy failed ...'
            Write(6,*) 'MC-PDFT returned with return code, rc = ',
     &                  iReturn
            Call Abend()
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
       nfld_tim  = 0
       nfld_stat = 0
*                                                                      *
************************************************************************
*                                                                      *
       Return
       End
