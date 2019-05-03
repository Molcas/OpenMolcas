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
      Subroutine ZMatReader(iZMUnit,LuWr,nAtoms,nXAtoms,nBasis,
     & nAskAtoms,iErr)
      Implicit Integer (i-n)
      Implicit Real*8 (a-h,o-z)
#include "g_zmatconv.fh"
      Character*80 Line, Blank
      Character*3  Command
      Character*24 Words(7)
* nAtoms : Total number of real atoms. Include X dummy atoms: NAT(i)= 0
* nXAtoms: Total number of ghost (Z) atoms:                   NAT(i)=-1
* nBasis : Nummer of atom types requiring Basis Set
C  ***  nAskAtoms.EQ.-1  =>  Seward ZMAT input  => Use "End of"
C  ***  nAskAtoms.NE.-1  =>  GateWay ZMAT input => Use nAskAtoms
      iErr = 0
      Blank = ' '
      nAtoms  = 0
      nXAtoms = 0
      nBasis  = 0
      Do i = 1, 100 ! MaxNat
        BasReq(i) = .False.
      EndDo
      Do i = 1, MaxAtoms
        Symbols(i) = '     '
        NAT(i) = 0
        iZmat(i,1) = 0
        iZmat(i,2) = 0
        iZmat(i,3) = 0
        Zmat(i,1) = 0.0d0
        Zmat(i,2) = 0.0d0
        Zmat(i,3) = 0.0d0
      EndDo

* Read Line (or COMMAND)
10    If ((nAtoms + nXAtoms).EQ.nAskAtoms) GoTo 100
      Read(iZMUnit,'(A)',Err=9906,End=9999) Line
      If ( Line(1:1).EQ.'*' ) GoTo 10
      If ( Line.EQ.Blank ) GoTo 100
      Command = Line(1:3)
      Call UpCase(Command)
      If (Command.EQ.'END') GoTo 100
      iErr = 0
      NA = 0
      NB = 0
      NT = 0
      Dist  = 0.0d0
      Beta  = 0.0d0
      Theta = 0.0d0

* Read Symbol           [ Symb ]
      Call Pick_Words(Line,7,Nwords,Words)
      If (Nwords.LT.1) GoTo 9993
      Call FoundAtomicNumber(LuWr,Words(1),NAtom,iErr)
      If (iErr.NE.0) GoTo 9998
      If (NAtom.GE. 0) nAtoms  = nAtoms  + 1
      If (NAtom.EQ.-1) nXAtoms = nXAtoms + 1
      NAT(nAtoms + nXAtoms) = NAtom
      Symbols(nAtoms + nXAtoms) = Trim(Words(1))
      If (NAtom.GT. 0) BasReq(NAtom)=.True.
      If ((nAtoms + nXAtoms).EQ.1) GoTo 10 ! Raed Only the First Atom

* Read Distance         [ Symb   NA Dist ]
      Call Pick_Words(Line,7,Nwords,Words)
      If (Nwords.LT.3) GoTo 9993
      Call Get_iNumber(Words(2),NA,iErr)
      If (iErr.NE.0) GoTo 9998
      If (NA.GE.(nAtoms + nXAtoms)) GoTo 9997
      Call Get_dNumber(Words(3),Dist,iErr)
      If (iErr.NE.0) GoTo 9998
      If (Dist.LE.0.0d0) GoTo 9996
      iZmat(nAtoms + nXAtoms, 1) = NA
      Zmat(nAtoms + nXAtoms, 1) = Dist
      If ((nAtoms + nXAtoms).EQ.2) GoTo 10 ! Raed Only the Second Atom

* Read Planar angle     [ Symb   NA Dist   NB Beta ]
      Call Pick_Words(Line,7,Nwords,Words)
      If (Nwords.LT.5) GoTo 9993
      Call Get_iNumber(Words(4),NB,iErr)
      If (iErr.NE.0) GoTo 9998
      If (NB.GE.(nAtoms + nXAtoms)) GoTo 9997
      Call Get_dNumber(Words(5),Beta,iErr)
      If (iErr.NE.0) GoTo 9998
      If (Beta.LE.0.0d0 .OR. Beta.GE.180.0d0) GoTo 9995
      iZmat(nAtoms + nXAtoms, 2) = NB
      Zmat(nAtoms + nXAtoms, 2) = Beta
      If (NA.EQ.NB) GoTo 9994
      If ((nAtoms + nXAtoms).EQ.3) GoTo 10 ! Raed Only the Second Atom

* Read Dihedral angle   [ Symb   NA Dist   NB Beta   NT Theta]
      Call Pick_Words(Line,7,Nwords,Words)
      If (Nwords.LT.7) GoTo 9993
      Call Get_iNumber(Words(6),NT,iErr)
      If (iErr.NE.0) GoTo 9998
      If (NT.GE.(nAtoms + nXAtoms)) GoTo 9997
      Call Get_dNumber(Words(7),Theta,iErr)
      If (iErr.NE.0) GoTo 9998
      iZmat(nAtoms + nXAtoms, 3) = NT
      Zmat(nAtoms + nXAtoms, 3) = Theta
      If (NA.EQ.NB .OR. NB.EQ.NT .OR. NA.EQ.NT) GoTo 9994
      GoTo 10

* Pre-check Basis Set consistency  BasReq: Atom requiring Basis Set
100   nBasis = 0
      Do i = 1, 100
        If (BasReq(i)) nBasis = nBasis + 1
      EndDo

      GoTo 9999

9906  iErr = 1
      Write(LuWr,*) ' [ZMatReader]: Unable to read z-matrix file !'
      GoTo 9999

9997  iErr = 1
      Write(LuWr,*) ' [ZMatReader]: Wrong index in line'
      Write(LuWr,*) '               ',Line
      GoTo 9999

9996  iErr = 1
      Write(LuWr,*) ' [ZMatReader]: Wrong distance in line'
      Write(LuWr,*) '               ',Line
      GoTo 9999

9995  iErr = 1
      Write(LuWr,*) ' [ZMatReader]: Wrong planar angle in line'
      Write(LuWr,*) '               ',Line
      GoTo 9999

9994  iErr = 1
      Write(LuWr,*) ' [ZMatReader]: Multiple index in line'
      Write(LuWr,*) '               ',Line
      GoTo 9999

9993  iErr = 1
      Write(LuWr,*) ' [ZMatReader]: Z-Matrix incomplete in line'
      Write(LuWr,*) '               ',Line
      GoTo 9999

9998  iErr = 1
      Write(LuWr,*) ' [ZMatReader]: Error in line'
      Write(LuWr,*) '               ',Line
      GoTo 9999

9999  Return
      End
