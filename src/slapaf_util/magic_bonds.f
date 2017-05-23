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
      Subroutine Magic_Bonds(Coor,nAtoms,iTabBonds,nBondMax,nBonds,
     &                       iTabAtoms,nMax)
      Implicit real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Coor(3,nAtoms)
      Integer iTabBonds(3,nBondMax), iTabAtoms(2,0:nMax,nAtoms)
#include "bondtypes.fh"
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*
      Fi_limit = (165.0D0/180.D0)*Pi
      CosFi_limit=Cos(Fi_limit)
      Do iBond = 1, nBonds
         iBondType=iTabBonds(3,iBond)
         If (iBondType.eq.vdW_Bond) Go To 100
         If (iBondType.gt.Magic_Bond) Go To 100
*
         Do iCase = 1, 2
*
            If (iCase.eq.1) Then
               iAtom=iTabBonds(1,iBond)
               jAtom=iTabBonds(2,iBond)
            Else
               iAtom=iTabBonds(2,iBond)
               jAtom=iTabBonds(1,iBond)
            End If
*
            xij=Coor(1,iAtom)-Coor(1,jAtom)
            yij=Coor(2,iAtom)-Coor(2,jAtom)
            zij=Coor(3,iAtom)-Coor(3,jAtom)
            rij=Sqrt(xij**2+yij**2+zij**2)
*
            nNeighbor_j = iTabAtoms(1,0,jAtom)
            Do kNeighbor = 1, nNeighbor_j
               kAtom = iTabAtoms(1,kNeighbor,jAtom)
               jBond = iTabAtoms(2,kNeighbor,jAtom)
               If (jBond.ge.iBond) Go To 99
               jBondType=iTabBonds(3,jBond)
               If (jBondType.eq.vdW_Bond) Go To 99
               If (jBondType.gt.Magic_Bond) Go To 99
*
*              If this is close to a linear system we will form
*              a magic bond between iAtom and kAtom
*
               xkj=Coor(1,kAtom)-Coor(1,jAtom)
               ykj=Coor(2,kAtom)-Coor(2,jAtom)
               zkj=Coor(3,kAtom)-Coor(3,jAtom)
               rkj=Sqrt(xkj**2+ykj**2+zkj**2)
*
               CosFi=(xij*xkj+yij*ykj+zij*zkj)/(rij*rkj)
#ifdef _DEBUG_
               Write(6,*) 'iAtom,jAtom,kAtom=',iAtom,jAtom,kAtom
               Write(6,*) 'CosFi,CosFi_limit=',CosFi,CosFi_limit
#endif

*
*------------> Observe that this limit should be coordinated with
*              the limit for excluding torsions in torsion_list!
*
               If (CosFi.le.CosFi_Limit) Then
*
#ifdef _DEBUG_
                  Write (6,*) 'Forming a "magic" bond'
                  Write (6,*) 'iAtom,kAtom=',iAtom,kAtom
#endif
*                 Double check that this bond is not included already.
*                 If that is the case just update the bond type if it
*                 is not classified already as a covalent bond.
*
                  Do kBond = 1, nBonds
                     If ( (iTabBonds(1,kBond).eq.iAtom .and.
     &                     iTabBonds(2,kBond).eq.kAtom) .or.
     &                    (iTabBonds(1,kBond).eq.kAtom .and.
     &                     iTabBonds(2,kBond).eq.iAtom) .and.
     &                     iTabBonds(3,kBond).ne.Covalent_Bond) Then
                         iTabBonds(3,kBond)=jAtom + Magic_Bond
                         Go To 99
                     End If
                  End Do
*
                  If (nBonds+1.gt.nBondMax) Then
                     Call WarningMessage(2,'Error in Magic_Bonds')
                     Write (6,*) 'Magic_Bonds: nBonds.gt.nBondMax'
                     Write (6,*) 'iTabBonds:'
                     Do kBond = 1, nBonds
                        Write (6,*)
                        Write (6,*) 'kBond=',kBond
                        Write (6,*)
                        Write (6,*) 'Atoms=',iTabBonds(1,kBond),
     &                                       iTabBonds(2,kBond)
                        Write (6,*) 'Bondtype:',iTabBonds(3,kBond)
                     End Do
                     Call Abend()
                  End If
                  nBonds = nBonds + 1
                  iTabBonds(1,nBonds)=Max(iAtom,kAtom)
                  iTabBonds(2,nBonds)=Min(iAtom,kAtom)
                  iTabBonds(3,nBonds)=jAtom + Magic_Bond
*
                  nNeighbor_i=iTabAtoms(1,0,iAtom)+1
                  If (nNeighbor_i.gt.nMax) Then
                     Call WarningMessage(2,'Error in Magic_Bonds')
                     Write (6,*)
     &                  'Magic_Bonds: nNeighbor_i.gt.nMax'
                     Write (6,*) 'iAtom=',iAtom
                     Write (6,*) 'nNeighbor_i=',nNeighbor_i
                     Write (6,*) 'nMax=',nMax
                     Call Abend()
                  End If
                  iTabAtoms(1,0,iAtom)=nNeighbor_i
                  iTabAtoms(1,nNeighbor_i,iAtom)=kAtom
                  iTabAtoms(2,nNeighbor_i,iAtom)=nBonds
*
                  nNeighbor_k=iTabAtoms(1,0,kAtom)+1
                  If (nNeighbor_k.gt.nMax) Then
                     Call WarningMessage(2,'Error in Magic_Bonds')
                     Write (6,*)
     &                  'Magic_Bonds: nNeighbor_k.gt.nMax'
                     Write (6,*) 'kAtom=',kAtom
                     Write (6,*) 'nNeighbor_k=',nNeighbor_k
                     Write (6,*) 'nMax=',nMax
                     Call Abend()
                  End If
                  iTabAtoms(1,0,kAtom)=nNeighbor_k
                  iTabAtoms(1,nNeighbor_k,kAtom)=iAtom
                  iTabAtoms(2,nNeighbor_k,kAtom)=nBonds
               End If
*
 99            Continue
            End Do   ! kNeighbor
*
         End Do      ! iCase
*
 100     Continue
      End Do         ! iBond
*
      Return
      End
