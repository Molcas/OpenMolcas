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
* Copyright (C) 2017, Valera Veryazov                                  *
************************************************************************
      Subroutine XMatReader(iZMUnit,LuWr,nAtoms,nXAtoms,nBasis,
     &  nAskAtoms, nxbas,xb_label,xb_bas,iErr)
      Implicit Integer (i-n)
      Implicit Real*8 (a-h,o-z)
#include "g_zmatconv.fh"
      Character*80 Line, Blank
      Character*3  Command
      Character*24 Words(7)
        character *(*) xb_label(*)
        character *(*) xb_bas(*)
       xb_label(1)=' '
       xb_bas(1)=' '
       nxbas=1
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
      Dist  = 0.0d0

*  Here we read number or a file.
      read(Line,*,err=666,end=666) NA
      IreadHere=1
      goto 667
666   continue
      Ireadhere=0
      iXU=iZMUnit+1
      call molcas_open(iXU,line)
      read(iXU,*) NA
667   continue
        if(Ireadhere.eq.1) then
            Read(iZMUnit,'(A)',Err=9906,End=9999) Line
        else
            Read(iXU,'(A)',Err=9906,End=9999) Line
        endif
      do i=1,NA
        if(Ireadhere.eq.1) then
            Read(iZMUnit,'(A)',Err=9906,End=9999) Line
        else
            Read(iXU,'(A)',Err=9906,End=9999) Line
        endif
      Call Pick_Words(Line,4,Nwords,Words)
      If (Nwords.LT.4) GoTo 9993
      Call FoundAtomicNumber(LuWr,Words(1),NAtom,iErr)
      If (iErr.NE.0) GoTo 9998
      If (NAtom.GE. 0) nAtoms  = nAtoms  + 1
      If (NAtom.EQ.-1) nXAtoms = nXAtoms + 1
      NAT(nAtoms + nXAtoms) = NAtom
      Symbols(nAtoms + nXAtoms) = Trim(Words(1))
      If (NAtom.GT. 0) BasReq(NAtom)=.True.

      Call Get_dNumber(Words(2),Dist,iErr)
      Zmat(nAtoms + nXAtoms, 1) = Dist
      Call Get_dNumber(Words(3),Dist,iErr)
      Zmat(nAtoms + nXAtoms, 2) = Dist
      Call Get_dNumber(Words(4),Dist,iErr)
      Zmat(nAtoms + nXAtoms, 3) = Dist

      enddo
      if(Ireadhere.eq.0) close(iXU)
* Pre-check Basis Set consistency  BasReq: Atom requiring Basis Set
100   nBasis = 0
      Do i = 1, 100
        If (BasReq(i)) nBasis = nBasis + 1
      EndDo

      GoTo 9999

9906  iErr = 1
      Write(LuWr,*) ' [XMatReader]: Unable to read x-matrix file !'
      GoTo 9999

9993  iErr = 1
      Write(LuWr,*) ' [XMatReader]: X-Matrix incomplete in line'
      Write(LuWr,*) '               ',Line
      GoTo 9999

9998  iErr = 1
      Write(LuWr,*) ' [XMatReader]: Error in line'
      Write(LuWr,*) '               ',Line
      GoTo 9999

9999  Return
      End
