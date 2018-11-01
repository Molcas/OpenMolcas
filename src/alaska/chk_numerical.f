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
      Subroutine Chk_Numerical(LuSpool,Numerical)
      Implicit Real*8 (A-H,O-Z)
      Logical Numerical,Is_Root_Set, DNG, KeepOld, Found
      Character KWord*180, Key*180, Get_Ln*180
      External Get_Ln
#include "nac.fh"
#include "alaska_root.fh"
*
      Call Qpg_iScalar('DNG',DNG)
      If (DNG) Then
         Call Get_iScalar('DNG',iDNG)
         Numerical = iDNG.eq.1
      Else
         Numerical = .False.
      End If
      LuWr=6
*
* Setting the defaults
*    iRoot     : Which root to optimize the geometry for
*    rDelta    : Displacements are chosen as r_nearest_neighbor * rDelta
      iRoot     = 1
      rDelta    = 0.0100D0
      NACstates(1)=0
      NACstates(2)=0
      isNAC     = .False.
      KeepOld   = .False.
      DefRoot   = .True.
      ForceNAC  = .False.
      iRlxRoot  = 1
      Auto      = .False.
*
* RASSCF will set the default root to the root that is relaxed.
*
      Is_Root_Set = .False.
      Call Qpg_iScalar('NumGradRoot',Is_Root_Set)
      If (Is_Root_Set) Then
         Call Get_iScalar('NumGradRoot',iRoot)
      End If
*
      Rewind(LuSpool)
      Call RdNLst_(LuSpool,'ALASKA',No_Input_OK)
      KWord=' &ALASKA'
 998  Read (LuSpool,'(A72)',END=997,ERR=988) Key
      KWord = Key
      Call UpCase(KWord)
      If (KWord(1:4) .eq. 'NUME') Then
         Numerical = .True.
         Goto 998
      Else If (KWord(1:4) .eq. 'ROOT') Then
         Key = Get_Ln(LuSpool)
         Call Get_I(1,iRoot,1)
         DefRoot = .False.
         Goto 998
      Else If (KWord(1:4) .eq. 'DELT') Then
         Key = Get_Ln(LuSpool)
         Call Get_F(1,rDelta,1)
         Goto 998
      Else If (KWord(1:4) .eq. 'NAC ') Then
         Key = Get_Ln(LuSpool)
         Call Get_I(1,NACstates,2)
         isNAC = .True.
         DefRoot = .False.
         Goto 998
      Else If (KWord(1:4) .eq. 'KEEP') Then
         KeepOld = .True.
         Goto 998
      Else If (KWord(1:4) .eq. 'AUTO') Then
         Auto = .True.
         Goto 998
      Else If (KWord(1:4) .eq. 'END ') Then
         Goto 997
      Else
         Goto 998
      End If
*
 988  Call WarningMessage(2,
     &               'Chk_Numerical: Error reading the input')
      Write (LuWr,'(A,A)') 'Last read line=',KWord
      Call Quit_OnUserError()
*
 997  Continue
*
      Call Get_iScalar('Grad ready',iGO)
      iGO = iAnd(iGO,Not(2**0))
      Call Put_iScalar('Grad ready',iGO)
*
* Put on the runfile which root and delta to use
*
      Call qpg_iScalar('Relax CASSCF root',Found)
      If (Found) Then
         Call Get_iScalar('Relax CASSCF root',iRoot0)
         Call Put_iScalar('NumGradRoot',iRoot)
         Call Put_iScalar('Relax CASSCF root',iRoot)
      Else
         iRoot0=0
      End If
      Call qpg_iScalar('Relax Original root',Found)
      If (.NOT.Found) Then
         Call qpg_iScalar('Relax Original root',iRoot)
      Else
         Call Get_iScalar('Relax Original root',iRoot1)
         If (iRoot1.eq.iRoot0) Then
            Call Put_iScalar('Relax Original root',iRoot)
         End If
      End If
      Call Put_dScalar('Numerical Gradient rDelta',rDelta)
*
* These are really input options for numerical_gradient
*
      Call Put_lScalar('Keep old gradient',KeepOld)
*
      Return
      End
