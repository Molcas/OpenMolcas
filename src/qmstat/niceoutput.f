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
      Subroutine NiceOutPut(EelP,Gam,Gamma,BetaBol)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "real.fh"
#include "warnings.fh"
#include "constants.fh"

      Parameter (Conver1=1.0d10*CONST_BOHR_RADIUS_IN_SI_)
      Parameter (Conver2=2.0d0*Pi/360.0d0)
      Character*3 EelP
      Character*40 Word1,Word2,Word3
      Logical Eq,Pr,It,Si,Cl,Qu
      External Len_TrimAO

*
*-- Enter.
*

*
*-- Make sure that the string is of correct length.
*
      ia=Len_TrimAO(EelP)
      If(Len_TrimAO(EelP).ne.3) then
        Write(6,*)'Illegal call to NiceOutPut'
        Call Quit(_RC_INTERNAL_ERROR_)
      Endif

*
*-- Check what type of output that is requested.
*
      Eq=.false.
      Pr=.false.
      It=.false.
      Si=.false.
      Cl=.false.
      Qu=.false.
      If(index(EelP,'E').ne.0) Eq=.true.
      If(index(EelP,'P').ne.0) Pr=.true.
      If(index(EelP,'I').ne.0) It=.true.
      If(index(EelP,'S').ne.0) Si=.true.
      If(index(EelP,'C').ne.0) Cl=.true.
      If(index(EelP,'Q').ne.0) Qu=.true.

*
*-- Start printing!
*
*
*-- With aid of concatenation, we here construct a header.
*
      Write(6,*)
      Write(6,*)
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - - - -'
     &//' - - - - - - - - - - - -'
      Write(6,*)'  *   *   *   *   *   *   *   *   *   *   *   *   *  '
     &//' *   *   *   *   *   * '
      Write(6,*)' - - - - - - - - - - - - - - - - - - - - - - - - - - '
     &//'- - - - - - - - - - - -'
      Write(6,*)
      If(It) then
        Write(Word1,*)'QMStat simulation commencing: '
      Else
        Write(Word1,*)'SampFile analysis commencing'
      Endif
      If(Cl) then
        Write(Word2,*)'All Classical '
      Elseif(Qu) then
        Write(Word2,*)'Combined Quantum-Classical '
      Else
        Write(Word2,*)' '
      Endif
      If(Eq) then
        Write(Word3,*)'Equilibration'
      Elseif(Pr) then
        Write(Word3,*)'Production'
      Else
        Write(Word3,*)' '
      Endif
      iM1=Len_TrimAO(Word1)
      iM2=Len_TrimAO(Word2)
      iM3=Len_TrimAO(Word3)
      Write(6,*)Word1(1:iM1)//Word2(1:iM2)//Word3(1:iM3)

*
*-- Now dump a lot of information.
*
      If(It) then
        Write(6,*)
        Write(6,*)
        Write(6,11)'*  Parameters of the calculation  *'
        Write(6,*)
        Write(6,12)'--Macroscopic quantities'
        Write(6,12)'  Temperature(K)      Pressure(Atm.)   Permitivity'
        Write(6,13)Temp,Pres,Diel
        Write(6,12)'--Maximal MC-Step parameters'
        Write(6,12)'  Translation(Ang.)   Rotation(deg.)   '
     &//'Cavity Radius(Ang.)'
        Write(6,13)delX*Conver1,delFi/Conver2,delR*Conver1
        Write(6,12)'--Configuration data'
        Write(6,12)'  Initial conf.       Writing conf.    MC-Steps'
        If(iNrIn.ge.0)
     &  Write(6,14)iNrIn,iNrUt,nMicro*nMacro
        If(iNrIn.lt.0)
     &  Write(6,15)'   Random/Input',iNrUt,nMicro*nMacro
        If(QmType(1:4).ne.'RASS') then
        Write(6,12)'--Hartree-Fock simulation data'
        Write(6,12)'  Total Occupation    Number of Orbitals'
        Write(6,16)iOcc1,iOrb(1)
        Else
        Write(6,12)'--Rassi state simulation data'
        Write(6,12)'  State interacting with solvent'
        If(.not.lCiSelect) then
        Write(6,17)nEqState
        Else
        Write(6,12)'   CI-select overlap option used'
        Endif
        Write(6,12)'  State threshold     Density threshold'
        If(MoAveRed.and.ContrStateB) then
        Write(6,18)ThrsCont,ThrsRedOcc
        Elseif(MoAveRed.and..not.ContrStateB) then
        Write(6,19)'   N/A',ThrsRedOcc
        Elseif(.not.MoAveRed.and.ContrStateB) then
        Write(6,20)ThrsCont,'N/A'
        Else
        Write(6,12)'   N/A                 N/A'
        Endif
        If(nLvlShift.ne.0) then
        Write(6,12)'  Level shift applied'
        Endif
        Endif
      Endif
      Write(6,*)
      Write(6,*)' - - - - - - - - - - - - - - - - - - - - - - - - - - '
     &//'- - - - - - - - - - - -'
      Write(6,*)'  *   *   *   *   *   *   *   *   *   *   *   *   *  '
     &//' *   *   *   *   *   * '
      Write(6,*)'- - - - - - - - - - - - - - - - - - - - - - - - - - -'
     &//' - - - - - - - - - - - -'
      Write(6,*)
      If(iT) then
      Write(6,*)
      Write(6,*)'Simulation progress.'
      Write(6,*)
      Endif

*
*-- Some formats
*
11    Format('               ',A)
12    Format('    ',A)
13    Format('    ',3(F10.4,'        '))
14    Format('    ',3(I8,'        '))
15    Format('    ',A,2(I8,'         '))
16    Format('    ',2(I5,'                '))
17    Format('    ',I5)
18    Format('     ',2(E11.4,'         '),' ')
19    Format('    ',A,'               ',E11.4)
20    Format('     ',E11.4,'           ',A)

*
*-- Tschuss
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real(Gam)
        Call Unused_real(Gamma)
        Call Unused_real(BetaBol)
      End If
      End


*-------------------------------------------------------------------------*
* A subroutine that emulates the len_trim of later Fortran versions, but  *
* that is missing in some Fortran 77 compilers.                           *
*-------------------------------------------------------------------------*
      Integer Function Len_TrimAO(String)
      Character*(*) String
      Do 15,i=Len(String),1,-1
        If(String(i:i).ne.' ') Go To 20
15    Continue
20    Continue
      Len_TrimAO=i
      Return
      End
