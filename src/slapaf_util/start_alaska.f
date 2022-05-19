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
      Subroutine Start_Alaska()
      use Slapaf_Parameters, only: Request_Alaska, Request_RASSI
      Implicit Real*8 (a-h,o-z)
#include "print.fh"
#include "nadc.fh"
      Character*100 ProgName, Get_ProgName
      Character*128 FileName
      Character*16 StdIn, JOB1, JOB2
      Character*8 Method
      Character(LEN=180):: Line
      Logical Exists
*                                                                      *
************************************************************************
*                                                                      *
*     Get the name of the module
*
      ProgName=Get_ProgName()
      Call Upcase(ProgName)
      Call LeftAd(ProgName)
*
      iEnd = 1
 99   If (ProgName(iEnd:iEnd).ne.' ') Then
         iEnd=iEnd+1
         Go To 99
      End If
      iEnd = Min(iEnd-1,5)
      FileName=progname(1:iend)//'INP'
*
      LuInput=11
      LuInput=IsFreeUnit(LuInput)
      Call StdIn_Name(StdIn)
      Call Molcas_Open(LuInput,StdIn)
*
      If (Request_RASSI) Then
*                                                                      *
************************************************************************
*                                                                      *
*     Request computation of overlaps.
*
         If (nPrint(1).ge.6) Then
            Write (6,*)
            Write (6,*)
     &        ' Slapaf requests the computation of overlaps first!'
            Write (6,*)
         End If
*
         Call Get_cArray('Relax Method',Method,8)
         If ((Method .eq. 'CASPT2').or.(Method .eq. 'RASPT2')) Then
           JOB1='JOBMIX'
         Else
           JOB1='JOBIPH'
         End If
         Call f_inquire('JOBAUTO',Exists)
         JOB2=JOB1
         If (Exists) JOB2='JOBAUTO'
*
         Write (LuInput,'(A)') '>ECHO OFF'
         Write (LuInput,'(A)') '> export SL_OLD_TRAP=$MOLCAS_TRAP'
         Write (LuInput,'(A)') '> export MOLCAS_TRAP=ON'
         Write (LuInput,'(A)') ' &RASSI &End'
         Write (LuInput,'(A)') 'StOverlaps'
         Write (LuInput,'(A)') 'NrOfJobIphs'
         Write (LuInput,'(A)') '  2 all'
         Write (LuInput,'(A)') 'IphNames'
         Write (LuInput,'(A)') '  '//Trim(JOB1)
         Write (LuInput,'(A)') '  '//Trim(JOB2)
         Write (LuInput,'(A)') ' End of Input'
         If (JOB1.eq.'JOBMIX') Then
           Write (LuInput,'(A)') '> copy $Project.JobMix JOBAUTO'
         Else
           Write (LuInput,'(A)') '> copy $Project.JobIph JOBAUTO'
         End If
         Write (LuInput,'(A)') '> export MOLCAS_TRAP=$SL_OLD_TRAP'
*
      Else If (Request_Alaska) Then
*                                                                      *
************************************************************************
*                                                                      *
*     Request computation of gradients.
*
         If (nPrint(1).ge.6) Then
            Write (6,*)
            Write (6,*)
     &        ' Slapaf requests the computation of gradients first!'
            If (iState(2).eq.0) Then
               Write (6,*) 'Root: ',iState(1)
            Else
               Write (6,*) 'Roots: ',iState(1),',',iState(2)
            End If
            Write (6,*)
         End If
*
         Write (LuInput,'(A)') '>ECHO OFF'
         Write (LuInput,'(A)') '>export SL_OLD_TRAP=$MOLCAS_TRAP'
         Write (LuInput,'(A)') '> export MOLCAS_TRAP=ON'
         Write (LuInput,'(A)') ' &Alaska &End'
         Write (LuInput,'(A)') 'AUTO'
         If (iState(2).ne.0) Then
            Write (LuInput,'(A)') 'NAC'
            Write (LuInput,'(I5,1X,I5)') iState(1),iState(2)
            Write (LuInput,'(A)') 'NoCSF'
         End If
         Write (LuInput,'(A)') ' End of Input'
         Write (LuInput,'(A)') '> export MOLCAS_TRAP=$SL_OLD_TRAP'
*        Repeat ECHO OFF, because it may have been turned on by alaska
         Write (LuInput,'(A)') '>ECHO OFF'
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Common code.
*
*
      Call f_inquire(Filename,Exists)
      If (Exists) Then
         LuSpool = 77
         LuSpool = IsFreeUnit(LuSpool)
         Call Molcas_Open(LuSpool, Filename)
*
 100     Continue
         Read (LuSpool,'(A)',End=900) Line
         Write(LuInput,'(A)') Line
         Go To 100
 900     Continue
*
         Close(LuSpool)
*
      Else
*
         Write (LuInput,'(A)') ' &Slapaf &End'
         Write (LuInput,'(A)') ' End of Input'
*
      End If
      Write (LuInput,'(A)') '>ECHO ON'
*
      Close(LuInput)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
