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
      Subroutine PolRas(iDist,iDistIm,iDT,iFI,iFP,iFil,iCStart
     &                 ,iTriState,VMat,Smat,DiFac,Ract,icnum,Energy
     &                 ,NVarv,iSTC,Haveri,iQ_Atoms,ip_ExpVal,Poli)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm2.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension Poli(MxQCen,10),VMat(MxStOT),FFp(nPol*nPart,3)
      Dimension VpolMat(MxStOT),Smat(MxStOT),RoMatSt(MxStOT)
      Dimension EEigen(MxState)
      Dimension iDT(3),iFI(3),iFP(3),iFil(MxQCen,10)
      Logical JaNej,Haveri

*
*-- Allocate and initialize the eigenvector matrix with the unit matrix.
*
      Call GetMem('Coeff','Allo','Real',iSTC,nState**2)
      call dcopy_(nState**2,[ZERO],iZERO,Work(iSTC),iONE)
      call dcopy_(nState,[ONE],iZERO,Work(iSTC),nState+1)
*
*-- Define some numbers.
*
      nQMCent=(iQ_Atoms*(iQ_Atoms+1))/2
      nPolCent=nPart*nPol
      Rinv=1.0d0/Ract
      R2inv=Rinv**2
      PolFac=DiFac/Ract**3
      Egun=0.0d0

*
*-- Compute distances and vectors needed for solving the polarization
*   equations. Observe the use of a memory allocator before the call
*   to polprep. This is to conserve memory but without having to
*   rewrite the entire polprep routine (written originally with
*   static allocations). A future project is to rewrite polprep.
*
      Call Memory_PolPrep('Allo',ixx,iyy,izz,irr3,ixxi,iyyi,izzi
     &                   ,iGri,nPol,nPart)
      Call PolPrep(iDist,iDistIM,Work(ixx),Work(iyy),Work(izz)
     &            ,Work(irr3),Work(ixxi),Work(iyyi),Work(izzi)
     &            ,Work(iGri),iCNum,nPolCent)

*
*-- Polarization loop commencing.
*
      NVarv=0
7912  Continue
      NVarv=NVarv+1
      Energy=0.0d0
      Call PolSolv(iDT,iFI,iFP,Work(ixx),Work(iyy),Work(izz),Work(irr3)
     &            ,Work(ixxi),Work(iyyi),Work(izzi),Work(iGri),FFp
     &            ,iCNum,r2Inv,DiFac,nPolCent)
      Call DensiSt(RomatSt,Work(iSTC),nEqState,nState,nState)
      Call Polins(Energy,iCall,nPolCent,nQMCent,iFil,VpolMat,FFp
     &           ,PolFac,poli,xyzMyQ,xyzMyI,xyzMyP,iCstart,iQ_Atoms
     &           ,qTot,ChaNuc,RoMatSt,xyzQuQ,CT)

*
*-- Assemble the Hamiltonian matrix.
*
      Do 815, i=1,iTriState
        HmatState(i)=HmatSOld(i)+Vmat(i)+VpolMat(i)+SMat(i)
815   Continue
      Energy=0.5*Energy

*
*-- Diagonalize the bastard. Eigenvalues are sorted and the
*   relevant eigenvalue is added to total energy.
*
      Call GetMem('Scratch','Allo','Real',iScratch,nState**2)
      Call Diag_Driver('V','A','L',nState,HMatState,Work(iScratch)
     &                ,nState,Dummy,Dummy,iDum,iDum,EEigen,Work(iSTC)
     &                ,nState,1,-1,'J',nFound,iErr)
      If(lCiSelect) Call CiSelector(nEqState,nState,iSTC,nCIRef,iCIInd
     &                                                         ,dCIRef)
      Energy=Energy+EEigen(nEqState)
      Call GetMem('Scratch','Free','Real',iScratch,nState**2)

*
*-- Check if polarization loop has converged.
*
      Call HaveWeConv(iCNum,iCStart,iQ_Atoms,nPolCent,iDT,FFp,xyzMyI
     &               ,Egun,Energy,NVarv,JaNej,Haveri)
      If(Haveri) GoTo 8108
      If(.not.JaNej) GoTo 7912

8108  Continue

*
*-- Deallocate.
*
      Call Memory_PolPrep('Free',ixx,iyy,izz,irr3,ixxi,iyyi,izzi,iGri
     &                   ,nPol,nPart)

*
*-- If expectation values are extracted, make a detour.
*
      If(lExtr(6))Call Expectus('RASSI',HmatSOld,Vmat,VpolMat,Smat
     &                         ,MxStOT,iSTC,nState,lExtr(4),iExtr_Eig
     &                         ,ip_ExpVal)

*
*-- Is it dead? It's terminated!
*
      Return
      End
