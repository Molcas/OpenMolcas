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
      Subroutine Connect_Fragments(nSet,nAtoms,iTabBonds,nBondMax,
     &                             nBonds,Coor,iTabAtoms,nMax,iANr)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "info_slapaf.fh"
      Real*8 Coor(3,nAtoms)
      Integer nSet(nAtoms), iTabAtoms(2,0:nMax,nAtoms), iANr(nAtoms),
     &        iTabBonds(3,nBondMax)
      Integer, Allocatable, Dimension(:) :: iStack
#include "bondtypes.fh"
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      dR_Thr=rFuzz
      HH_Thr=5.0D0
*
      Call mma_allocate(iStack,nAtoms,label="iStack")
*
 1    Continue
      Not_Defined=-1
      Call ICopy(nAtoms,[Not_Defined],0,nSet,1)
*
      iX = 0
      nStack=1
      iAtom = 1
      iStack(nStack)=iAtom
      iSet = 1
      mBonds=0
      nSet(iAtom)=iSet
*
#ifdef _DEBUGPRINT_
      Write (6,*) '>>>> Connect Fragments <<<<<'
      Call RecPrt('Coor',' ',Coor,3,nAtoms)
#endif
*                                                                      *
************************************************************************
*                                                                      *
 99   Continue
*
*                                                                      *
************************************************************************
*                                                                      *
*     Pick up an atom from the top of the stack.
*
      iAtom=iStack(nStack)
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'iAtom=',iAtom
      Write (6,*) 'iSet=',iSet
      Write (6,*) 'nStack=',nStack
      Write (6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
#endif
*
*     Loop over all atoms to which this atom has a bond
*
*     Decrement the stack size one step
*
      nStack=nStack-1
      Do iNeighbor = 1, iTabAtoms(1,0,iAtom)
         jAtom=iTabAtoms(1,iNeighbor,iAtom)
         iBond=iTabAtoms(2,iNeighbor,iAtom)
#ifdef _DEBUGPRINT_
         Write (6,*) 'jAtom=',jAtom
         Write (6,*) 'iBond=',iBond
         Write (6,*) 'iTabBonds(3,iBond)=',iTabBonds(3,iBond)
         Write (6,*) 'nSet(jAtom)=',nSet(jAtom)
#endif
*        Skip if vdW bond.
         If (iTabBonds(3,iBond).eq.vdW_Bond) Go To 101
         If (nSet(jAtom).eq.Not_Defined) Then
*
*           If not defined add to the stack
            nStack = nStack + 1
            iStack(nStack)=jAtom
*           Assign which set the atom belongs to.
            nSet(jAtom)=iSet
         Else If (nSet(jAtom).ne.iSet) Then
*
*           Error if atom belongs to the wrong set.
            Call WarningMessage(2,' Error in Connect_Fragments')
            Write (6,*) 'Connect_Fragments:'
            Write (6,*) 'Inconsistent set indices!'
            Write (6,*) 'iSet,jSet=',iSet,nSet(jAtom)
            Call Abend()
         End If
 101     Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*     If there are more atoms in the stack process the next one.
*
      If (nStack.ne.0) Go To 99
*                                                                      *
************************************************************************
*                                                                      *
*---- If all atoms now belong to a set we are done. Otherwise
*     we have a new set.
*
*     Loop over all atoms. If an atom is found which is not assigned
*     to a set add it to the stack.
*
      Do iAtom = 1, nAtoms
         If (nSet(iAtom).eq.Not_Defined) Then
*           Add atom to the stack
            nStack=1
            iStack(nStack)=iAtom
*           Set the new set index and assign it to the atom.
            iSet = iSet + 1
            nSet(iAtom)=iSet
            Go To 99
         End If
      End Do
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     Now all atoms belong to a set
*                                                                      *
************************************************************************
************************************************************************
*                                                                      *
*     If all atoms belong to the same set then we are done.
*
      If (iSet.eq.1) Then
         Call mma_deallocate(iStack)
         Return
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Now we have to connect the disjoint sets.
*                                                                      *
************************************************************************
*                                                                      *
*     Find the shortest distance between every two sets
*
      Call Allocate_Work(ipSetDist,iSet*iSet)
      call dcopy_(iSet*iSet,[1.0D6],0,Work(ipSetDist),1)
      Do kAtom = 1, nAtoms
         kSet = nSet(kAtom)
         Do lAtom = kAtom+1, nAtoms
            lSet = nSet(lAtom)
            If (lSet.eq.kSet) Go To 501
*
*           For l < k:
*           In (k,l), the minimum distance between the sets k and l
*           In (l,k), the minimum distance excluding H-H
*
            iOff = (Min(kSet,lSet)-1)*iSet + Max(kSet,lSet)-1
            iOffHH = (Max(kSet,lSet)-1)*iSet + Min(kSet,lSet)-1
            RTest = (Coor(1,kAtom)-Coor(1,lAtom))**2
     &             +(Coor(2,kAtom)-Coor(2,lAtom))**2
     &             +(Coor(3,kAtom)-Coor(3,lAtom))**2
            If (iAnr(kAtom).ne.1 .or. iANr(lAtom).ne.1)
     &         Work(ipSetDist+iOff) = Min(Work(ipSetDist+iOff),RTest)
            Work(ipSetDist+iOffHH) = Min(Work(ipSetDist+iOffHH),RTest)
 501        Continue
         End Do
      End Do
*
*     Find the shortest distance between any two sets.
*
      rShort = 1.0D5
      Do kSet = 1, iSet
         Do lSet = kSet+1, iSet
*
*           H-H distances are only considered if the minimum distance
*           (excluding H-H) is greater than HH_Thr
            iOff = (kSet-1)*iSet + lSet-1
            If (Work(ipSetDist+iOff).gt.HH_Thr**2)
     &         iOff = (lSet-1)*iSet + kSet-1
            If (Work(ipSetDist+iOff).lt.rShort) Then
               rShort = Work(ipSetDist+iOff)
            End If
         End Do
      End Do
*
*     Add bonds between sets if they are shorter than rShort+dR_Thr
*
      Do kAtom = 1, nAtoms
         kSet = nSet(kAtom)
         Do lAtom = kAtom+1, nAtoms
            lSet = nSet(lAtom)
            If (lSet.eq.kSet) Go To 502
*
*           Again, add H-H bonds only if necessary
            iOff = (Min(kSet,lSet)-1)*iSet + Max(kSet,lSet)-1
            If ((Work(ipSetDist+iOff).le.HH_Thr**2) .and.
     &          (iANr(kAtom).eq.1 .and. iANr(lAtom).eq.1)) Go To 502
*
            RTest = (Coor(1,kAtom)-Coor(1,lAtom))**2
     &             +(Coor(2,kAtom)-Coor(2,lAtom))**2
     &             +(Coor(3,kAtom)-Coor(3,lAtom))**2
            RTest = Sqrt(RTest)
            If (RTest.le.Sqrt(rShort)+dR_Thr) Then
               Do iBond = 1, nBonds
*
*                 Look through the bond list and find if it is there.
                  If ( (iTabBonds(1,iBond).eq.kAtom .and.
     &                  iTabBonds(2,iBond).eq.lAtom ) .or.
     &                 (iTabBonds(1,iBond).eq.lAtom .and.
     &                  iTabBonds(2,iBond).eq.kAtom ) ) Then
                     iTabBonds(3,iBond) = Fragments_Bond
*
*                    If the bond is already in the list that is it.
                     Go To 503
                  End If
               End Do
*
*              Add the bond to the bond list
               If (nBonds+1.gt.nBondMax) Then
                  Call WarningMessage(2,' Error in Connect_Fragments')
                  Write (6,*) 'Connect_Fragments: nBonds+1.gt.nBondMax'
                  Call Abend()
               End If
*
               nBonds = nBonds + 1
               iTabBonds(1,nBonds) = lAtom
               iTabBonds(2,nBonds) = kAtom
               iTabBonds(3,nBonds) = Fragments_Bond
*
*              Update atoms list
*
               nNeighbor_k = iTabAtoms(1,0,kAtom)+1
               If (nNeighbor_k.gt.nMax) Then
                  Call WarningMessage(2,' Error in Connect_Fragments')
                  Write (6,*) 'Connect_Fragments: nNeighbor_k.gt.nMax'
                  Write (6,*) 'kAtom=',kAtom
                  Write (6,*) 'nNeighbor_k=',nNeighbor_k
                  Write (6,*) 'nMax=',nMax
                  Call Abend()
               End If
               iTabAtoms(1,0,kAtom) = nNeighbor_k
               iTabAtoms(1,nNeighbor_k,kAtom) = lAtom
               iTabAtoms(2,nNeighbor_k,kAtom) = nBonds
*
               nNeighbor_l = iTabAtoms(1,0,lAtom)+1
               If (nNeighbor_l.gt.nMax) Then
                  Call WarningMessage(2,' Error in Connect_Fragments')
                  Write (6,*) 'Connect_Fragments: nNeighbor_l.gt.nMax'
                  Write (6,*) 'lAtom=',lAtom
                  Write (6,*) 'nNeighbor_l=',nNeighbor_l
                  Write (6,*) 'nMax=',nMax
                  Call Abend()
               End If
               iTabAtoms(1,0,lAtom) = nNeighbor_l
               iTabAtoms(1,nNeighbor_l,lAtom) = kAtom
               iTabAtoms(2,nNeighbor_l,lAtom) = nBonds
*
 503           Continue
            End If
 502        Continue
         End Do
      End Do
      Call Free_Work(ipSetDist)
*
*     Try again!
      Go To 1
*                                                                      *
************************************************************************
*                                                                      *
      End
