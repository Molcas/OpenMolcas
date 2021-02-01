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
*  ZMatrixConverter
*
*> @brief
*>   Routine for reading seward input in Z-Matrix format
*> @author Giovanni Ghigo
*>
*> @details
*> The input for seward is read and a string vector is generate (\p STDINP).
*> This vector contains a standard seward input and it will be read later as usual
*> by a modified copy of the section ``BASI`` code already present in ::RdCtl_Seward.
*> This new code is in the ::StdSewInput routine.
*> Only the standard basis present in the ``$MOLCAS/basis_library`` are allowed.
*>
*> @param[in]    LuRd    Input file unit number
*> @param[in]    LuWr    Output file unit number
*> @param[in]    mxAtom  Parameter
*> @param[out]   STDINP  String vector of seward standard input
*> @param[out]   lSTDINP Length of String vector \p STDINP
*> @param[out]   iErr    Error flag
************************************************************************
       Subroutine ZMatrixConverter(LuRd,LuWr,mxAtom,STDINP,lSTDINP,
     & iglobal,nxbas,xb_label,xb_bas,iErr)
************************************************************************
* Author: Giovanni Ghigo                                               *
*         Torino (Italy)  October-November 2006                        *
*                                                                      *
* This is an adaptation of Program ZMatrixConverter                    *
* A converter of Z-Matrix in cartesian coordinates in MolCAS format.   *
* Version 1.0                                                          *
* The input for seward is read and a string vector is generate. This   *
* vector is a standard seward input and it will be read later as usual *
* by a copy on the code already present in RdCtl_Seward.               *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Character*180 STDINP(mxAtom*2)
      Character*180 aDebug
      Character*12 Angstring
#include "constants.fh"
#include "g_zmatconv.fh"
      Logical IfTest
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
      Call BasisReader(LuWr,nBase,iglobal,
     &      nxbas,xb_label,xb_bas,iErr)
      If (IfTest) then
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'ZMatrixConverter - From BasisReader :'
        Write(LuWr,*) '                   nBase=',nBase
        Do i = 1, Num_Elem
          If (BasAva(i)) Write(LuWr,'(I23,3X,A)') i, Base(i)
        EndDo
        Write(LuWr,*)
      EndIf
      If (iErr.NE.0) GoTo 9905
      iErr = 0
      Call ZMatReader(LuRd,LuWr,nAtoms,nXAtoms,nBasis,-1,iErr)
      If (IfTest) then
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'ZMatrixConverter - From ZMatReader :'
        Write(LuWr,*) '                   nAtoms=',nAtoms,
     &     ', nXAtoms=',nXAtoms,', Tot=',nAtoms+nXAtoms
        Write(LuWr,*)
     &' Label NA    i   bond        j   angle       k   dihedral'
        Do i = 1, nAtoms+nXAtoms
          Write(LuWr,'(1X,A,I3,3(1X,I4,1X,F10.6))')
     &    Symbols(i),NAT(i),(iZmat(i,j),Zmat(i,j),j=1,3)
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

      Call Put_iScalar('N ZMAT',nAtoms+nXAtoms)
      Call Put_cArray('Symbol ZMAT',Symbols(1),(nAtoms+nXAtoms)*5)
      Call Put_iArray('Index ZMAT',iZmat,MaxAtoms*3)
      Call Put_iArray('NAT ZMAT',NAT,nAtoms+nXAtoms)

* Calculate coordinates
      torad = CONST_PI_ / 180.0d0
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
2000  Continue
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
        Write(STDINP(iSTDINP),'(A5,3F16.10,A)')
     &          Symbols(i),(Coords(i,j),j=1,3),Angstring
        iSTDINP = iSTDINP + 1
2100    Continue
      EndDo
      Write(STDINP(iSTDINP),'(A)') 'End of basis'
      lSTDINP = iSTDINP
      If (IfTest) then
        Write(LuWr,*)
        Write(LuWr,*) '------------------------------------------------'
        Write(LuWr,*) 'ZMatrixConverter - The input passed to SEWARD : '
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
