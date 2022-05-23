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
      Subroutine ExScf(iCStart,nBaseQ,nBaseC,nCnC_C,iQ_Atoms
     &                ,nAtomsCC,Ax,Ay,Az,itri,Smat,SmatPure,InCutOff
     &                ,ipAOSum)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm1.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension Smat(MxOT),CorTemp(3),nCnC_C(MxBasC)
      Dimension SmatPure(MxOT)
      Dimension Inside(MxAt,3)
      Logical Inside,NearBy,InCutOff

*----------------------------------------------------------------------*
* Deduce how much the QM-molecule is translated from its position as   *
* definied in Seward.                                                  *
*----------------------------------------------------------------------*
      Ax=Cordst(1,1)-outxyz(1,1)
      Ay=Cordst(1,2)-outxyz(1,2)
      Az=Cordst(1,3)-outxyz(1,3)
*----------------------------------------------------------------------*
* Rotate solvent orbitals and make AO integration.                     *
*----------------------------------------------------------------------*
      Cut_ExSq1=Cut_Ex1**2
      Cut_ExSq2=Cut_Ex2**2
      nOrbSize=iOrb(1)*iOrb(2)
      nV2size=iOrb(2)*nBaseC
      nAObaseSize=nBaseQ*nBaseC
      nStorlek=iOrb(1)*nBaseC
      Call GetMem('RotOrb','Allo','Real',iV2,nV2size)
      Call GetMem('Sint','Allo','Real',ipAOint,nAObaseSize)
      Call GetMem('Sintpar','Allo','Real',ipAOintpar,nAObaseSize)
      Call GetMem('OvlMO','Allo','Real',iOvlMO,nOrbSize)
      Call GetMem('Intermed','Allo','Real',iInte,nStorlek)
      Call GetMem('OvlMOpure','Allo','Real',iOPure,nOrbSize)
      Call GetMem('OvlMOene','Allo','Real',iOvlMOE,nOrbSize)
      Call GetMem('AUX','Allo','Real',ipAUX,iOrb(1)**2)
      Call GetMem('AUXp','Allo','Real',ipAUXp,iOrb(1)**2)
      Call GetMem('AUXtri','Allo','Real',ipAUXtri,iTri)
*Jose****************************************************************
      If(lExtr(8)) then
         nAOqMOcl=nBaseQ*iOrb(2)
         iAOAOTri=nBaseQ*(nBaseQ+1)/2
         Call GetMem('qAOclMOOvl','Allo','Real',iAOMOOvl,nAOqMOcl)
         Call GetMem('qAOclMOOvlE','Allo','Real',iAOMOOvlE,nAOqMOcl)
         Call GetMem('AuxAOp','Allo','Real',ipAOAUX,nBaseQ**2)
         Call GetMem('AuxAOpTri','Allo','Real',ipAOAUXtri,iAOAOTri)
      Endif
*     Call GetMem('SumOvlAOQ','Allo','Real',ipAOSum,iAOAOTri)
*********************************************************************
      InCutOff=.false.
      Do 575, i=1,iTri
        Smat(i)=0
        SmatPure(i)=0
575   Continue
      nInsideCut=0
      Do 501, N=iCStart-1,nCent*(nPart-1),nCent
*
*-- Initialize.
*
        dist_sw=1D20
        r3=1D20
        Do 11, i=1,MxAt
          Inside(i,1)=.false.
          Inside(i,2)=.false.
          Inside(i,3)=.false.
11      Continue
        NearBy=.false.
*-----Loop over atoms.
        Do 511, inwm=1,iQ_Atoms
          Do 512, k=1,3
            CorTemp(k)=(Cordst(N+1,k)-Cordst(inwm,k))**2
512       Continue
          r2=CorTemp(1)+CorTemp(2)+CorTemp(3)
          dist_sw=min(dist_sw,r2)
          DH1=0.0d0 !Distances for the inner cut-off. Also include the
          DH2=0.0d0  !hydrogens.
          Do 513, k=1,3
            DH1=DH1+(Cordst(N+2,k)-Cordst(inwm,k))**2
            DH2=DH2+(Cordst(N+3,k)-Cordst(inwm,k))**2
513       Continue
          r3temp1=min(DH1,DH2)
          r3temp2=min(r3temp1,r2)
          r3=min(r3,r3temp2)
*----Check if this atom-atom pair inside
          If(r2.lt.Cut_ExSq1) then
            Inside(inwm,1)=.true.
            NearBy=.true.
          Endif
          If(DH1.lt.Cut_ExSq1) then
            Inside(inwm,2)=.true.
            NearBy=.true.
          Endif
          If(DH2.lt.Cut_ExSq1) then
            Inside(inwm,3)=.true.
            NearBy=.true.
          Endif
511     Continue
*
*-- Now make the cut-off test.
*
        If(.not.NearBy) Go To 501
        If(r3.lt.Cut_ExSq2) then !Inner cut-off.
          InCutOff=.true.
        Endif
        nInsideCut=nInsideCut+1

*
*-- Make the AO-AO overlap integration.
*
        Call AOIntegrate(iCStart,nBaseQ,nBaseC,Ax,Ay,Az,nCnC_C
     &          ,iQ_Atoms,nAtomsCC,ipAOint,ipAOintpar,iV2,N,lmax
     &          ,Inside)

*
*-- Transform to MO-MO overlap.
*
        Call Dgemm_('T','N',iOrb(1),nBaseC,nBaseQ,ONE,Work(iV1)
     &    ,nBaseQ,Work(ipAOint),nBaseQ,ZERO,Work(iInte),iOrb(1))
        Call Dgemm_('N','N',iOrb(1),iOrb(2),nBaseC,ONE,Work(iInte)
     &    ,iOrb(1),Work(iV2),nBaseC,ZERO,Work(iOvlMO),iOrb(1))
        call dcopy_(nOrbSize,[ZERO],iZERO,Work(iOvlMOE),iONE)
        Do 5751, i=1,iOrb(2)
          ind=iOrb(1)*(i-1)
          Call DaxPy_(iOrb(1),c_orbene(i),Work(iOvlMO+ind),iONE
     &              ,Work(iOvlMOE+ind),iONE)
5751    Continue
*
***Jose
      If(lExtr(8)) then
        Call Dgemm_('N','N',nBaseQ,iOrb(2),nBaseC,ONE,Work(ipAOint)
     &    ,nBaseQ,Work(iV2),nBaseC,ZERO,Work(iAOMOOvl),nBaseQ)
        call dcopy_(nAOqMOcl,[ZERO],iZERO,Work(iAOMOOvlE),iONE)
        Do 5761, i=1,iOrb(2)
          ind=nBaseQ*(i-1)
          Call DaxPy_(nBaseQ,c_orbene(i),Work(iAOMOOvl+ind),iONE
     &              ,Work(iAOMOOvlE+ind),iONE)
5761    Continue
      Endif
********************
*
*-- If you are interested, print some bla bla bla.
*
        If(iPrint.ge.29) then
          Write(6,*)
          Write(6,*)'OVERLAP BETWEEN QM-SYSTEM AND SOLVENT MOLECULE '
     &,N/nCent
          Write(6,*)'QM-MO  SOLV-MO  OVERLAP'
          call dcopy_(iOrb(1)*iOrb(2),Work(iOvlMO),iONE
     &              ,Work(iOPure),iONE)
          Do 54599, i=0,iOrb(1)-1
            Do 54699, j=0,iOrb(2)-1
              Write(6,8888)i+1,j+1,Work(iOPure+i+j*iOrb(1))
54699       Continue
54599     Continue
8888      Format(I3,'    ',I3,'       ',F12.10)
        Endif

*
*-- Construct the perturbation and accumulate pure overlap for
*   subsequent construction of higher order term.
*
        Call Dgemm_('N','T',iOrb(1),iOrb(1),iOrb(2),exrep2
     &            ,Work(iOvlMO),iOrb(1),Work(iOvlMOE),iOrb(1),ZERO
     &            ,Work(ipAUX),iOrb(1))
        Call SqToTri_Q(Work(ipAUX),Work(ipAUXtri),iOrb(1))
        Call DaxPy_(iTri,ONE,Work(ipAUXtri),iONE,Smat,iONE)
        Call Dgemm_('N','T',iOrb(1),iOrb(1),iOrb(2),ONE
     &            ,Work(iOvlMO),iOrb(1),Work(iOvlMO),iOrb(1),ZERO
     &            ,Work(ipAUXp),iOrb(1))
        Call SqToTri_Q(Work(ipAUXp),Work(ipAUXtri),iOrb(1))
        Call DaxPy_(iTri,ONE,Work(ipAUXtri),iONE,SmatPure,iONE)

*Jose*********************************
      If(lExtr(8)) then
        Call Dgemm_('N','T',nBaseQ,nBaseQ,iOrb(2),exrep2
     &            ,Work(iAOMOOvl),nBaseQ,Work(iAOMOOvlE),nBaseQ,ZERO
     &            ,Work(ipAOAUX),nBaseQ)
        Call SqToTri_Q(Work(ipAOAUX),Work(ipAOAUXtri),nBaseQ)
        Call DaxPy_(iAOAOTri,ONE,Work(ipAOAUXtri),iONE,Work(ipAOSum)
     &            ,iONE)
      Endif
**************************************
*
*-- This solvent molecule ends now!
*
501   Continue

      Call GetMem('RotOrb','Free','Real',iV2,nV2size)
      Call GetMem('Sint','Free','Real',ipAOint,nAObaseSize)
      Call GetMem('Sintpar','Free','Real',ipAOintpar,nAObaseSize)
      Call GetMem('OvlMO','Free','Real',iOvlMO,nOrbSize)
      Call GetMem('Intermed','Free','Real',iInte,nStorlek)
      Call GetMem('OvlMOpure','Free','Real',iOPure,nOrbSize)
      Call GetMem('OvlMOene','Free','Real',iOvlMOE,nOrbSize)
      Call GetMem('AUX','Free','Real',ipAUX,iOrb(1)**2)
      Call GetMem('AUXp','Free','Real',ipAUXp,iOrb(1)**2)
      Call GetMem('AUXtri','Free','Real',ipAUXtri,iTri)
*Jose****************************************************************
      If(lExtr(8)) then
        Call GetMem('qAOclMOOvl','Free','Real',iAOMOOvl,nAOqMOcl)
        Call GetMem('qAOclMOOvlE','Free','Real',iAOMOOvlE,nAOqMOcl)
        Call GetMem('AuxAOp','Free','Real',ipAOAUX,nBaseQ**2)
        Call GetMem('AuxAOpTri','Free','Real',ipAOAUXtri,iAOAOTri)
      Endif
*********************************************************************

      Return
      End
