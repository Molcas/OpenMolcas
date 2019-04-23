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
* Copyright (C) 2006, Giovanni Ghigo                                   *
************************************************************************
       Subroutine ZMatrixConverter_GW(LuRd,LuWr,LuOut,nAskAtoms,iErr)
************************************************************************
* Author: Giovanni Ghigo                                               *
*         Torino (Italy)  October-November 2006                        *
*                                                                      *
* This is an adaptation of Subroutine ZMatrixConverter for GateWay     *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "g_zmatconv.fh"

      nAtoms = 0
      nBase  = 0
      Do i = 1, Num_Elem
        Base(i)   = ' '
        BasAva(i) = .False.
        BasReq(i) = .False.
      EndDo
      Do i = 1, MaxAtoms
        Coords(i,1) = 0.0d0
        Coords(i,2) = 0.0d0
        Coords(i,3) = 0.0d0
      EndDo
      nBasis = 0
      iErr   = 0

* Reading input
      iErr = 0
      Call ZMatReader(LuRd,LuWr,nAtoms,nXAtoms,nBasis,nAskAtoms,iErr)
      If (iErr.NE.0) GoTo 9906
      Write(LuOut,*) nAtoms
      Write(LuOut,'(A)') 'Angstrom'

* Some checks
      If (nAtoms.EQ.0) then
        Write(LuWr,*) 'ERROR: No atom coordinates specified !'
        GoTo 9999
      EndIf

      Call Put_iScalar('N ZMAT',nAtoms+nXAtoms)
      Call Put_cArray('Symbol ZMAT',Symbols(1),(nAtoms+nXAtoms)*5)
      Call Put_iArray('Index ZMAT',iZmat,MaxAtoms*3)
      Call Put_iArray('NAT ZMAT',NAT,nAtoms+nXAtoms)

* Calculate coordinates
      torad = 3.14159265358979323846d0 / 180.0d0
*     Atom #1
      If (nAtoms+nXAtoms.EQ.1) GoTo 2000
*     Atom #2
      Coords(2,3)=Zmat(2,1)  ! Z(2)=R
      If (nAtoms+nXAtoms.EQ.2) GoTo 2000
*     Atom #3
      If (iZmat(3,1).EQ.1) then
        Coords(3,1)=Zmat(3,1)*SIN(Zmat(3,2)*torad) ! X(2)=R sin(A)
        Coords(3,3)=Zmat(3,1)*COS(Zmat(3,2)*torad) ! Z(3)=R cos(A)
      else
        Coords(3,1)=Zmat(3,1)*SIN(Zmat(3,2)*torad)
        Coords(3,3)=Coords(2,3)-Zmat(3,1)*COS(Zmat(3,2)*torad)
      EndIf
      If (nAtoms+nXAtoms.EQ.3) GoTo 2000
*     Atom #4 ->
      Do iAtom = 4, nAtoms+nXAtoms
        Call ZMatConv(LuWr,iAtom,iErr)
      EndDo
      If (iErr.NE.0) GoTo 9999

*     Check for superposed atoms
2000  Do i = 1, nAtoms+nXAtoms
        If (NAT(i).GT.0) then
          Do j = i+1, nAtoms+nXAtoms
            If (NAT(j).GT.0) then
              r = 0.0d0
              Do k =1, 3
                r = r + (Coords(i,k)-Coords(j,k))**2
              EndDo
              If (r.LT.0.0001d0) GoTo 9907
            EndIf
          EndDo
        EndIf
      EndDo

* Writing

      Do i = 1, nAtoms+nXAtoms
        If (NAT(i).GT.0) Write(LuOut,999) Symbols(i),(Coords(i,k),k=1,3)
      EndDo
999   Format(A5,1X,3(F12.6))
      GoTo 9999

9906  Write(LuWr,*) ' ERROR: Wrong input in Z-Matrix definition !'
      GoTo 9999
9907  Write(LuWr,*) ' ERROR: Superimposed atoms: ',i,j,'  r=',sqrt(r)
      GoTo 9999

9999  Return
      End
