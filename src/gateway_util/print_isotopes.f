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
* Copyright (C) 2018, Ignacio Fdez. Galvan                             *
************************************************************************
      SubRoutine Print_Isotopes()
      use Period
      use Basis_Info, only: nCnttp, dbsc
      Implicit Real*8 (A-H,O-Z)
#include "print.fh"
#include "constants2.fh"
      Logical Changed
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      iPrint = nPrint(iRout)
      If (iPrint.eq.0) Return
      Call qEnter('Print_Isotopes')
      LuWr=6
*                                                                      *
************************************************************************
*                                                                      *
*     Determine whether any mass is different from default
*
      Changed = .False.
      Do i=1,nCnttp
        If (dbsc(i)%Aux.or.dbsc(i)%Frag) Cycle
        iAtom=dbsc(i)%AtmNr
        If (dbsc(i)%CntMass.ne.rMass(iAtom)) Then
          Changed = .True.
          Exit
        End If
      End Do
      If (.not.Changed.and.iPrint.le.5) Return
*                                                                      *
************************************************************************
*                                                                      *
      Write (LuWr,*)
      Call CollapseOutput(1,'   Isotope specification:')
      Write (LuWr,'(3X,A)') '   ----------------------'
      Write (LuWr,*)
      If (Changed) Then
        Write(LuWr,10) 'Center                     [     Default     ]'
        Write(LuWr,10) 'Type   Z    A    mass (Da) [   A    mass (Da)]'
        Write(LuWr,10) '---------------------------------------------'
      Else
        Write(LuWr,10) 'Center'
        Write(LuWr,10) 'Type   Z    A    mass (Da)'
        Write(LuWr,10) '--------------------------'
      End If
      Do i=1,nCnttp
        If (dbsc(i)%Aux.or.dbsc(i)%Frag) Cycle
        iAtom=dbsc(i)%AtmNr
        act_Mass=dbsc(i)%CntMass/UToAU
        def_Mass=rmass(iAtom)/UToAU
        If (act_Mass.ne.def_Mass) Then
          Write(LuWr,101) i,iAtom,nInt(act_Mass),act_Mass,
     &                            nInt(def_Mass),def_Mass
        Else
          Write(LuWr,100) i,iAtom,nInt(act_Mass),act_Mass
        End If
      End Do
      Call CollapseOutput(0,'   Isotope specification:')
      Write (LuWr,*)
*                                                                      *
************************************************************************
*                                                                      *
 10   Format(1X,A)
100   Format(I5,1X,I3,1X,I4,1X,F12.6)
101   Format(I5,1X,I3,1X,I4,1X,F12.6,1X,'[',I4,1X,F12.6,']')
      Call qExit('Print_Isotopes')
      Return
      End
