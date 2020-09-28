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
      Subroutine Bond_Tester(Coor,nAtoms,iTab,nx,ny,nz,ix,iy,iz,iAtom,
     &                       iRow,iANr,Schlegel,iOptC,iTabBonds,nBonds,
     &                       nBondMax,iTabAtoms,nMax,ThrB,ThrB_vdW)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#define _VDW_
#include "ddvdt.fh"
      Real*8 Coor(3,nAtoms)
      Integer iTab(0:nMax,nx,ny,nz), iANr(nAtoms),
     &        iTabBonds(3,nBondMax), iTabAtoms(2,0:nMax,nAtoms)
      Logical Help, Schlegel
      Real*8 A(3), B(3)
#include "bondtypes.fh"
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     Check box indices for consistency
*
      If (ix.lt.1.or.ix.gt.nx) Return
      If (iy.lt.1.or.iy.gt.ny) Return
      If (iz.lt.1.or.iz.gt.nz) Return
      Nr=iTab(0, ix,iy,iz) ! nr of atoms in the box.
      If (Nr.eq.0) Return
*
      iRow=iTabRow(iANr(iAtom))
      nVal_i=0
      nn=iTabAtoms(1,0,iAtom)
      Do i = 1, nn
         If (iTabBonds(3,iTabAtoms(2,i,iAtom)).eq.Covalent_Bond)
     &       nVal_i=nVal_i+1
      End Do
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) 'Bond_Tester: iAtom,=',iAtom
      Write (6,*) '                Box(ix,iy,iz)=(',ix,iy,iz,')'
      Write (6,*) '                        nVal_i=',nVal_i
      Write (6,*)
#endif
*
*     Loop over all atoms in the box
*
      If (ThrB.lt.ThrB_vdw) ivdW=vdW_Bond
      Box: Do Ir = 1, Nr
         jAtom=iTab(Ir,ix,iy,iz)
C        If (iAtom.le.jAtom) Cycle Box
         If (iAtom.ge.jAtom) Cycle Box
         jRow =iTabRow(iANr(jAtom))
         Help = iRow.gt.3.or.jRow.gt.3
#ifdef _DEBUG_
         Write (6,*) ' jAtom, iAnr(jAtom)=',jAtom, iAnr(jAtom)
         Write (6,*) 'Help=',Help
#endif
*
         x = Coor(1,iAtom)-Coor(1,jAtom)
         y = Coor(2,iAtom)-Coor(2,jAtom)
         z = Coor(3,iAtom)-Coor(3,jAtom)
         A(:)=Coor(:,iAtom)-Coor(:,jAtom)
         ANorm = dot_product(A,A)
         rij2 = x**2 + y**2 + z**2
         r0 = rAv(iRow,jRow)
         alpha=aAv(iRow,jRow)
*
*--------Test if we have a bond iAtom-jAtom
*
         If (Schlegel.or.Help) Then
            Rab=Sqrt(rij2)
            RabCov=CovRad(iANr(iAtom))+CovRad(iANr(jAtom))
#ifdef _DEBUG_
            Write (6,*) 'Rab=',Rab
            Write (6,*) CovRad(iANr(iAtom)),CovRad(iANr(jAtom))
#endif
            If (Rab.le.1.25d0*RabCov) Then
*
*              covalent bond
*
               test=1.0D0
               test_vdW=0.0D0
*
*              Skip if we are looking for vdW bonds
*
               If (ThrB.gt.ThrB_vdW) Cycle Box
            Else If (Rab.gt.1.25d0*RabCov .and.
     &               Rab.le.2.00D0*RabCov) Then
*
*              vdW's bond
*
               test=0.0D0
               test_vdW=ThrB_VdW
            Else
*
*              No Bond!
*
               Cycle Box
            End If

         Else

            test=Exp(alpha*(r0**2-rij2))
            If (iAnd(iOptC,2048).eq.2048) Then
               r0_vdW=r_ref_vdW(iRow,jRow)
               test_vdW=Exp(-alpha_vdW*(r0_vdW-SQRT(rij2))**2)
            Else
               r0_vdW=0.0D0
               test_vdW=0.0D0
            End If
#ifdef _DEBUG_
            Write (6,*)
            Write (6,*) 'Bond_Tester: iAtom,jAtom=',iAtom,jAtom
            Write (6,*) 'Bond_Tester: iRow, jRow =',iRow,jRow
            Write (6,*) 'Bond_Tester: test=',test,ThrB
            Write (6,*) 'Bond_Tester: test_vdW=',test_vdW,ThrB_vdW
            Write (6,*) 'Bond_Tester: r0_vdW=',r0_vdW, SQRT(rij2)
#endif
*
*           If the valence force constant small but not too small
*           denote the bond as an vdW's bond.
            Test_vdW=Max(Test_vdW,Test)
*
*           If already valence bond skip if also vdW bond. We picked
*           up this bond before!
*
            If (test.ge.ThrB .and. Test_vdW.ge.ThrB_vdW) Cycle Box
*
*           If none skip
*
            If (test.lt.ThrB .and. test_vdW.lt.ThrB_vdW) Cycle Box
*
*           Some logic to see if vdw bond should be included.
*           Hydrogen-hydrogen is always included.
*
            If (iANr(iAtom).ne.1 .or. iANr(jAtom).ne.1) Then
*
*              Skip if any of the atoms has more than 6 valence bonds
*              and the other at least 1 valence bond.
               nVal_j=0
               nn=iTabAtoms(1,0,jAtom)
               Do i = 1, nn
                  If ( iTabBonds(3,iTabAtoms(2,i,jAtom)) .eq.
     &                 Covalent_Bond) nVal_j=nVal_j+1
               End Do
#ifdef _DEBUG_
               Write (6,*) 'nVal_j=',nVal_j
#endif
               If ((nVal_i.ge.6.and.nVal_j.ge.1) .or.
     &             (nVal_j.ge.6.and.nVal_i.ge.1)     ) Cycle Box
            End If
#define _OLD_CODE
#ifndef _OLD_CODE_
*
*           We need to exclude vdW bonds if there is an atom close to
*           being in between the two atoms being considered and
*           forming a covalent bond.
*
            If (test.lt.ThrB) Then ! only in case of vdW bond
*
*              Loop over all covalently bonded neighbors of atom iAtom
*
#ifdef _DEBUG_
               Write (6,*) 'Test validity of vdW bond'
#endif
               Do i = 1, iTabAtoms(1,0,iAtom)
                  kAtom=iTabAtoms(1,i,iAtom)
                  iBond=iTabAtoms(2,i,iAtom)
                  iType=iTabBonds(3,iBond)
                  If (iType.ne.Covalent_Bond) Cycle
                  B(:)=Coor(:,iAtom)-Coor(:,kAtom)
                  BNorm = dot_product(B,B)
                  CosPhi=dot_product(A,B)/Sqrt(ANorm*BNorm)
                  Phi=ACOS(CosPhi)
#ifdef _DEBUG_
                  Write (6,*) 'kAtom,Phi=',kAtom,Phi
#endif
                  If (Abs(Phi).lt.Pi/Four) Cycle Box
               End Do
            End If
#endif
*
         End If
*
         If (nBonds+1.gt.nBondMax) Then
            Write (6,*) 'Bond_Tester: nBonds+1.gt.nBondMax'
            Write (6,*) 'nBonds+1=',nBonds+1
            Write (6,*) 'nBondMax=',nBondMax
            Call Abend()
         End If
*
         nBonds = nBonds + 1
         iTabBonds(1,nBonds)=iAtom
         iTabBonds(2,nBonds)=jAtom
*
         If (test.ge.ThrB) Then
            ivdW=Covalent_Bond
         Else If (test_vdW.ge.ThrB_vdW) Then
            ivdW=vdW_Bond
         Else
            Write (6,*) 'Bond_Tester: Illegal operation'
            Call Abend()
            ivdW=99
         End If
         iTabBonds(3,nBonds)=ivdW
#ifdef _DEBUG_
         Write (6,*)
         Write (6,*) 'Add bond to the list.'
         Write (6,*) 'Atom pair:',iAtom,jAtom
         Write (6,*) 'Bond type:', Bondtype(Min(3,ivdW))
         Write (6,*)
#endif
*
         nNeighbor=iTabAtoms(1,0,iAtom)
         If (nNeighbor+1.gt.nMax) Then
            Write (6,*) 'Bond_Tester(1): nNeighbor+1.gt.nMax'
            Write (6,*) 'iAtom=',iAtom
            Write (6,*) 'nNeighbor=',nNeighbor
            Write (6,*) 'nMax=',nMax
            Call Abend()
         End If

*        Update the neighbor lists of atoms iAtom and jAtom
         nNeighbor=nNeighbor+1
         iTabAtoms(1,0,iAtom) = nNeighbor
         iTabAtoms(1,nNeighbor,iAtom)=jAtom
         iTabAtoms(2,nNeighbor,iAtom)=nBonds
*
         nNeighbor=iTabAtoms(1,0,jAtom)
         If (nNeighbor+1.gt.nMax) Then
            Write (6,*) 'Bond_Tester(2): nNeighbor+1.gt.nMax'
            Write (6,*) 'jAtom=',jAtom
            Write (6,*) 'nNeighbor=',nNeighbor
            Write (6,*) 'nMax=',nMax
            Call Abend()
         End If
         nNeighbor=nNeighbor+1
         iTabAtoms(1,0,jAtom) = nNeighbor
         iTabAtoms(1,nNeighbor,jAtom)=iAtom
         iTabAtoms(2,nNeighbor,jAtom)=nBonds
*
      End Do Box
*
      Return
      End
