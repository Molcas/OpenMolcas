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
      Subroutine Process_Coor(R,Coor,nAtoms,nSym,iOper)
      Implicit Real*8 (a-h,o-z)
      Real*8 R(3), Coor(3,*)
      Integer iOper(0:nSym-1)
*
*     Local array
*
      Real*8 Q(3)
*                                                                      *
************************************************************************
*                                                                      *
C     Call RecPrt('Coor(Enter)',' ',Coor,3,nAtoms)
C     Call RecPrt('R',' ',R,3,1)
*                                                                      *
************************************************************************
*                                                                      *
*     Identify if this is a new center.
*
      Do iAtom = 1, nAtoms
         If (R(1).eq.Coor(1,iAtom) .and.
     &       R(2).eq.Coor(2,iAtom) .and.
     &       R(3).eq.Coor(3,iAtom) ) Return
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     Add this atom to the list
*
      nAtoms=nAtoms+1
      call dcopy_(3,R,1,Coor(1,nAtoms),1)
      iRef = nAtoms
C     Call RecPrt('Coor(updated)',' ',Coor,3,nAtoms)
C     Write (*,*) 'nSym=',nSym
C     Write (*,*) 'iOper=',iOper
*                                                                      *
************************************************************************
*                                                                      *
*     Add symmetry degenerate atoms to the list
*
      Do iSym = 1, nSym-1
C        Write (6,*) 'iOper(iSym)=',iOper(iSym)
         call dcopy_(3,R,1,Q,1)
         If (iAnd(iOper(iSym),1).ne.0) Q(1)=-Q(1)
         If (iAnd(iOper(iSym),2).ne.0) Q(2)=-Q(2)
         If (iAnd(iOper(iSym),4).ne.0) Q(3)=-Q(3)
C        Call RecPrt('Q',' ',Q,3,1)
         Do iAtom = iRef, nAtoms
            If (Q(1).eq.Coor(1,iAtom) .and.
     &          Q(2).eq.Coor(2,iAtom) .and.
     &          Q(3).eq.Coor(3,iAtom) ) Go To 100
         End Do
         nAtoms = nAtoms + 1
         call dcopy_(3,Q,1,Coor(1,nAtoms),1)
 100     Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
