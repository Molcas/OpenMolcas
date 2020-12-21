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
* Copyright (C) Valera Veryazov                                        *
************************************************************************
       Subroutine XMatrixConverter(LuRd,LuWr,mxAtom,STDINP,lSTDINP,
     & iglobal,nxbas,xb_label,xb_bas,iErr)
************************************************************************
* Author: Valera Veryazov                                              *
*                                                                      *
* This is an adaptation of GG Program ZMatrixConverter                 *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Character*180 STDINP(mxAtom*2)
      Character*180 aDebug
      Character*12 Angstring
#include "g_zmatconv.fh"
      Logical IfTest
      Character ll*5
      Character lll*4
        character *(*) xb_label(*)
        character *(*) xb_bas(*)

#ifdef _DEBUGPRINT_
      IfTest=.True.
#else
      IfTest=.False.
#endif

C  ***  H-Fm (Atomic numbers 1-100)
C  ***  X dummy atoms (NA = 0 )
C  ***  Z ghost atoms (NA =-1 )
C  ***  nAskAtoms.EQ.-1  =>  Seward ZMAT input
C  ***  nAskAtoms.NE.-1  =>  GateWay ZMAT input

C nAtoms : nr. of atoms passed to SEWARD (includes X dummy atoms).
C nXAtoms: nr. of ghost Z atoms (not passed to SEWARD but resumed
C          by OutZMat in SLAPAF).
C nBase  : number of BasisSets found in input.
C Base(i): BasisSet for atom with Atomic Number -i-.
C BasAva(i) & BasReq(i): Logical to check BasisSet-consistency.
C Coords(_,i): X, Y, Z, coordinates (in Angstrom) for atom -i-.

      nAtoms  = 0
      nXAtoms = 0
      nBase  = 0
      lSTDINP = 0
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
      Angstring = '  / Angstrom'

* Reading input
      iErr = 0
      Call BasisReader(LuWr,nBase,iglobal,nxbas,
     *     xb_label,xb_bas,iErr)
      If (IfTest) then
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'XMatrixConverter - From BasisReader :'
        Write(LuWr,*) '                   nBase=',nBase
        Do i = 1, Num_Elem
          If (BasAva(i)) Write(LuWr,'(I23,3X,A)') i, Base(i)
        EndDo
        Write(LuWr,*)
      EndIf
      If (iErr.NE.0) GoTo 9905
      iErr = 0
      Call XMatReader(LuRd,LuWr,nAtoms,nXAtoms,nBasis,-1,
     & nxbas,xb_label,xb_bas,iErr)
      If (IfTest) then
        Do i = 1, nAtoms+nXAtoms
          Write(LuWr,'(1X,A,I3,3(F10.6))')
     &    Symbols(i),NAT(i),(Zmat(i,j),j=1,3)
        EndDo
        Write(LuWr,*)
      EndIf
      If (iErr.NE.0) GoTo 9906

* Some checks
      If (nBase.EQ.0) then
        Write(LuWr,*) 'ERROR: No basis set specified !'
        GoTo 9999
      EndIf
      If (nAtoms.EQ.0) then
        Write(LuWr,*) 'ERROR: No atom coordinates specified !'
        GoTo 9999
      EndIf
      If (nBase.lt.nBasis) then
        Write(LuWr,*) 'ERROR: Wrong number of basis sets !'
        Write(LuWr,*) '       Available=',nBase,'  Required=',nBasis
        GoTo 9999
      EndIf
      Call BasisConsistency(LuWr,iErr)
      If (iErr.NE.0) then
        Write(LuWr,*) 'ERROR: Basis set inconsistency !'
        GoTo 9999
      EndIf
      Call Put_iScalar('N ZMAT',0)
      Do i=1,nAtoms+nXAtoms
      Coords(i,1)=Zmat(i,1)
      Coords(i,2)=Zmat(i,2)
      Coords(i,3)=Zmat(i,3)
      enddo

      If (IfTest) then
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'ZMatrixConverter - XYZCoords (Angstroms) :'
        Do i = 1, nAtoms+nXAtoms
          Write(LuWr,99) i,NAT(i),(Coords(i,j),j=1,3)
        EndDo
        Write(LuWr,*)
      EndIf

*     Check for superposed atoms
      Do i = 1, nAtoms+nXAtoms
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
c2000  Continue
      If (NAT(1).EQ.-1) then
        NATprev = -1
      else
        NATprev = -9999
      EndIf
      iSTDINP = 1
      Do i = 1, nAtoms+nXAtoms
        If (NAT(i).EQ.-1) GoTo 2100
        If (NAT(i).NE.NATprev) then
          If (i.NE.1 .and. NATprev.NE.-1) then
             Write(STDINP(iSTDINP),'(A)') 'End of basis'
             iSTDINP = iSTDINP + 1
          EndIf
          Write(STDINP(iSTDINP),'(A)') 'Basis set'
          iSTDINP = iSTDINP + 1
          If (NAT(i).GT.0) then
            Write(STDINP(iSTDINP),'(A)') Base(NAT(i))
          else
            Write(STDINP(iSTDINP),'(A)') 'X..... / InLine '
            iSTDINP = iSTDINP + 1
            Write(STDINP(iSTDINP),'(A)') '0.   0'
            iSTDINP = iSTDINP + 1
            Write(STDINP(iSTDINP),'(A)') '0    0'
          EndIf
          iSTDINP = iSTDINP + 1
          NATprev = NAT(i)
        EndIf
        ll=Symbols(i)
        write(lll,'(i4)') i
        ineedfix=1
        do i5=1,5
         if(index('0123456789',ll(i5:i5)).ne.0) ineedfix=0
        enddo
        if(ineedfix.eq.1) then
         i5=index(ll,' ')
         do i4=1,4
           if(lll(i4:i4).ne.' ') then
               ll(i5:i5)=lll(i4:i4)
               i5=i5+1
           endif
         enddo
        endif
        Write(STDINP(iSTDINP),'(A5,3F24.18,A)')
     &          ll,(Coords(i,j),j=1,3),Angstring
        iSTDINP = iSTDINP + 1
2100    Continue
      EndDo
      Write(STDINP(iSTDINP),'(A)') 'End of basis'
      lSTDINP = iSTDINP
      If (IfTest) then
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'XMatrixConverter - The input passed to SEWARD : '
        Do i = 1, iSTDINP
          aDebug = STDINP(i)
          Write(LuWr,*) aDebug
        EndDo
        Write(LuWr,*)
      EndIf
      GoTo 9999

99    Format(I3,1X,I3,1X,3(F12.6))

9905  Write(LuWr,*) ' ERROR: Wrong input in Bases Set definition !'
      GoTo 9999
9906  Write(LuWr,*) ' ERROR: Wrong input in Z-Matrix definition !'
      GoTo 9999
9907  Write(LuWr,*) ' ERROR: Superimposed atoms: ',i,j,'  r=',sqrt(r)
      GoTo 9999

9999  Continue

      End
