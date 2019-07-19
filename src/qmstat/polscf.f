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
      Subroutine PolScf(iDist,iDistIm,iDT,iFI,iFP,iFil,iCStart
     &                 ,iTri,VMat,Smat,DiFac,Ract,icnum,Energy
     &                 ,NVarv,iMOC,Haveri,iQ_Atoms,ip_ExpVal,Poli)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm1.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension Poli(MxQCen,10),VMat(MxOT),FFp(nPol*nPart,3)
      Dimension VpolMat(MxOT),Smat(MxOT),RoMat(MxOT)
      Dimension EEigen(MxOrb)
      Dimension iDT(3),iFI(3),iFP(3),iFil(MxQCen,10)
      Logical JaNej,Haveri

*
*-- Allocate and initialize the eigenvector matrix with the unit matrix.
*
      iOrba=iOrb(1)
      Call GetMem('Coeff','Allo','Real',iMOC,iOrba**2)
      call dcopy_(iOrba**2,[ZERO],iZERO,Work(iMOC),iONE)
      call dcopy_(iOrba,[ONE],iZERO,Work(iMOC),iOrba+1)

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
*-- Compute some distances needed for the polarization. See polras for
*   some additional information on this.
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
      Call Densi_MO(Romat,Work(iMOC),1,iOcc1,iOrb(1),iOrb(1))
      If(Mp2DensCorr) Call DCorrCorr(Romat,DenCorrD,Trace_MP2
     &                              ,iOrb(1),iOcc1)
      Call Polink(Energy,iCall,nPolCent,nQMCent,iFil,VpolMat,FFp
     &          ,PolFac,Poli,iCstart,iTri,iQ_Atoms,qTot,ChaNuc,xyzMyQ
     &          ,xyzMyI,xyzMyP,Romat,xyzQuQ,CT)

*
*-- Construct the Fock-matrix from two-electron super-matrix and
*   one-electron matrix, with solvent perturbations added.
*
      Do 801, i=1,iTri
        FockM(i)=0.0d0
        Do 802, j=1,iTri
          FockM(i)=FockM(i)+RoMat(j)*Work(iSupM+iTri*(i-1)+j-1)
802     Continue
        OneEl=HHmat(i)+Vmat(i)+VpolMat(i)+Smat(i)
        FockM(i)=FockM(i)+OneEl
*---- See Szabo-Ostlund eq. 3.184.
        Energy=Energy+(FockM(i)+OneEl)*RoMat(i)
801   Continue
*-- Add potential-nuclear energy.
      Energy=Potnuc+Energy*0.5

*
*-- If energy is strange, scream!
*
      If(Energy.gt.0) then
        Write(6,*)
        Write(6,*)'  SCF energy is positive. Serious error somewhere.'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif

*
*-- Diagonalize the Fock-matrix. Eigenvalues are sorted.
*
      Call GetMem('Scratch','Allo','Real',iScratch,iOrba**2)
      Call Diag_Driver('V','A','L',iOrba,FockM,Work(iScratch)
     &                ,iOrba,Dummy,Dummy,iDum,iDum,EEigen,Work(iMOC)
     &                ,iOrba,1,-1,'J',nFound,iErr)
      Call GetMem('Scratch','Free','Real',iScratch,iOrba**2)

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
      If(lExtr(6))Call Expectus('SCF  ',HHmat,Vmat,VpolMat,Smat
     &                         ,MxOT,iMOC,iOrba,.false.,iOcc1
     &                         ,ip_ExpVal)

*
*-- The end is near, hold me!
*

      Return
      End
