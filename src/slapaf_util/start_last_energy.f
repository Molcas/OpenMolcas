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
      Subroutine Start_Last_Energy()
      Implicit Real*8 (a-h,o-z)
      Character*16 StdIn
      Character*8 Method
#include "print.fh"
#include "stdalloc.fh"
      Logical Saddle, FoundLastEn
      Real*8, Dimension (:,:), Allocatable :: CN
      Interface
        Subroutine Get_Coord_New(CN,nCoord)
        Real*8, Dimension (:,:), Allocatable :: CN
        Integer nCoord
        End Subroutine
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      If (nPrint(1).ge.6) Then
         Write (6,*)
         Write (6,*) ' Slapaf requests the last energy to be computed!'
         Write (6,*)
      End If
*
      LuInput=11
      LuInput=IsFreeUnit(LuInput)
      Call StdIn_Name(StdIn)
      Call Molcas_Open(LuInput,StdIn)

      Call Get_iScalar('Grad ready',iGO)
*
      Write (LuInput,'(A)') '>ECHO OFF'
      Write (LuInput,'(A)') '>export SL_OLD_TRAP=$MOLCAS_TRAP'
      Write (LuInput,'(A)') '>export MOLCAS_TRAP=ON'
*
      Call qpg_dArray('Saddle',Saddle,nSaddle)
      If (Saddle) Then
*
*         Clean up among the runfiles, etc.
*
          Write (LuInput,'(A)')
     &          '>COPY $Project$SubProject.RunFile $Project.RunFile'

          Call qpg_cArray('LastEnergyMethod',FoundLastEn,lengthlast)
          If (FoundLastEn) Then
             Call Get_cArray('LastEnergyMethod',Method,8)
          Else
             Call Get_cArray('Relax Method',Method,8)
          EndIf

          If (Method .eq. 'CASSCF'   .or.
     &        Method .eq. 'RASSCF'   .or.
     &        Method .eq. 'CASSCFSA' .or.
     &        Method .eq. 'RASSCFSA' .or.
     &        Method .eq. 'CASPT2'   .or.
     &        Method .eq. 'RASPT2'   ) Then

             Write (LuInput,'(A)')
     &             '>COPY $Project$SubProject.JobIph $Project.JobIph'
             Write (LuInput,'(A)') '>RM $Project.Reac.JobIph'
             Write (LuInput,'(A)') '>RM $Project.Prod.JobIph'
          End If
          Write (LuInput,'(A)') '>RM $Project.Reac.RunFile'
          Write (LuInput,'(A)') '>RM $Project.Prod.RunFile'
          Write (LuInput,'(A)') '>export SubProject='
          Write (LuInput,'(A)') '>export MOLCAS_SADDLE=0'
*
*         Put the final structure as the reference structure
*
          Call Get_Coord_New(CN,nCoord)
          Call NameRun('RUNREAC')
          Call Put_dArray('Ref_Geom',CN,3*nCoord)
          Call NameRun('RUNPROD')
          Call Put_dArray('Ref_Geom',CN,3*nCoord)
          Call NameRun('RUNFILE')
          Call mma_deallocate(CN)
      End If
*
      If (iGo.le.1) Then
*
*        Normal behaviour
*
         Write (LuInput,'(A)') ' &Last_Energy &End'
         Write (LuInput,'(A)') 'End of Input'
*
      Else
*
*        A CI ISC search has been performed. For the moment I care only
*        for the ISC. Once the rest is done I will take care of that as
*        well.
*
         Write (LuInput,'(A)') '>COPY $OldProject.Seward.Input '
     &                       //'State1.Seward.Input'
         Write (LuInput,'(A)') '>COPY $OldProject.Seward.Input '
     &                       //'State2.Seward.Input'
         Write (LuInput,'(A)')'>COPY $OldProject.RunFile State1.RunFile'
         Write (LuInput,'(A)')'>COPY $OldProject.RunFile State2.RunFile'
         Write (LuInput,'(A)') '>RM molcas.env'
         Write (LuInput,'(A)') '>export Project=State1'
         Write (LuInput,'(A)') ' &Last_Energy &End'
         Write (LuInput,'(A)') 'End of Input'
         Write (LuInput,'(A)') '>RM molcas.env'
         Write (LuInput,'(A)') '>export Project=State2'
         Write (LuInput,'(A)') ' &Last_Energy &End'
         Write (LuInput,'(A)') 'End of Input'
      End If
      Write (LuInput,'(A)') '>export MOLCAS_TRAP=$SL_OLD_TRAP'
*
      Write (LuInput,'(A)') '>ECHO ON'

      Close(LuInput)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
