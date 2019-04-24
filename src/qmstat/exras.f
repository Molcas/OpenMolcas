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
      Subroutine ExRas(iCStart,nBaseQ,nBaseC,nCnC_C,iQ_Atoms
     &                ,nAtomsCC,Ax,Ay,Az,itristate,SmatRas,SmatPure
     &                ,InCutOff,ipAOSum)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"
#include "qmcom.fh"
#include "qm2.fh"
#include "numbers.fh"
#include "WrkSpc.fh"

      Dimension SmatRas(MxStOT),CorTemp(3)
      Dimension nCnC_C(MxBasC),SmatPure(MxStOT)
      Dimension Inside(MxAt,3)
      Logical Inside,NearBy,InCutOff

*----------------------------------------------------------------------*
* Deduce how much the QM-molecule is translated from its position as   *
* definied in Seward.                                                  *
*----------------------------------------------------------------------*
      Ax=Cordst(1,1)-outxyzRAS(1,1)
      Ay=Cordst(1,2)-outxyzRAS(1,2)
      Az=Cordst(1,3)-outxyzRAS(1,3)

*
*-- Make some initializations.
*
      Cut_ExSq1=Cut_Ex1**2
      Cut_ExSq2=Cut_Ex2**2
      nOrbSize=iOrb(1)*iOrb(2)
      nV2size=iOrb(2)*nBaseC
      nAObaseSize=nBaseQ*nBaseC
      nStorlek=iOrb(1)*nBaseC
      Call GetMem('RotOrb','Allo','Real',iV2,nV2size)
      Call GetMem('Sint','Allo','Real',ipAOint,nAObaseSize)
      Call GetMem('Sintpar','Allo','Real',ipAOintpar,nAObaseSize)
      nHalf=nBaseQ*iOrb(2)
      nGross=nBaseQ*(nBaseQ+1)/2
      Call GetMem('HalfTrans','Allo','Real',iHalfpar,nHalf)
      Call GetMem('HalfPure','Allo','Real',iHalf,nHalf)
      Call GetMem('HalfOrbE','Allo','Real',iHalfE,nHalf)
      Call GetMem('Auxilary','Allo','Real',ipAux,nBaseQ**2)
      Call GetMem('AuxilaryP','Allo','Real',ipAuxp,nBaseQ**2)
      Call GetMem('GammaAO','Allo','Real',ipAOG,nGross)
      Call GetMem('Accumulate','Allo','Real',ipACC,nBaseQ**2)
      Call GetMem('Accumulate','Allo','Real',ipACCp,nBaseQ**2)
      Call GetMem('AccumulateT','Allo','Real',ipACCt,nGross)
      Call GetMem('AccumulateTP','Allo','Real',ipACCtp,nGross)
      Call GetMem('TEMP','Allo','Real',iTEMP,nRedMO*iOrb(2))
      call dcopy_(nBaseQ**2,[ZERO],iZERO,Work(ipACC),iONE)
      call dcopy_(nBaseQ**2,[ZERO],iZERO,Work(ipACCp),iONE)
*Jose****************************************************************
      If(lExtr(8)) then
         Call GetMem('qAOclMOOvl','Allo','Real',iAOMOOvl,nHalf)
         Call GetMem('qAOclMOOvlE','Allo','Real',iAOMOOvlE,nHalf)
         Call GetMem('AuxAOp','Allo','Real',ipAOAUX,nBaseQ**2)
         Call GetMem('AuxAOpTri','Allo','Real',ipAOAUXtri,nGross)
      Endif
*********************************************************************
      InCutOff=.false.
      Do 575, i=1,iTriState
        SmatRas(i)=0
        SmatPure(i)=0
575   Continue
      nInsideCut=0
      If(MoAveRed) then
        nDim1=nRedMO
        nDim2=iOrb(2)
      Else
        nDim1=nBaseQ
        nDim2=iOrb(2)
      Endif
      nDimP=nDim1*nDim2
      nDimT=nDim1*(nDim1+1)/2
*
*-- Start loop over all solvent molecules.
*
      Do 501, N=iCStart-1,nCent*(nPart-1),nCent

*
*-- Initialize
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
*-----See if these atom-atom pairs inside cut-off.
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
*-----Make some cut-off tests.
*
        If(.not.NearBy) Go To 501  !If all distances larger than
                              !cut-off, jump to new solvent.
        If(r3.lt.Cut_ExSq2) then !Inner cut-off. Set flag to true then
          InCutOff=.true.        !huge energy is added later.
        Endif                    !S*S matrix, however!
        nInsideCut=nInsideCut+1

*
*-----Start integrating.
*
        Call AOIntegrate(iCStart,nBaseQ,nBaseC,Ax,Ay,Az,nCnC_C
     &          ,iQ_Atoms,nAtomsCC,ipAOint,ipAOintpar,iV2,N,lmax
     &          ,Inside)

*
*-- Transform overlaps from solvent AO to solvent MO.
*
        Call Dgemm_('N','N',nBaseQ,iOrb(2),nBaseC,ONE
     &      ,Work(ipAOint),nBaseQ,Work(iV2),nBaseC,ZERO
     &      ,Work(iHalf),nBaseQ)
*Jose***************************************************************
        If(lExtr(8)) then
          call dcopy_(nHalf,Work(iHalf),iONE,Work(iAOMOOvl),iONE)
          call dcopy_(nHalf,[ZERO],iZERO,Work(iAOMOOvlE),iONE)
          Do 5731, k=1,iOrb(2)
              ind=nBaseQ*(k-1)
              Call DaxPy_(nBaseQ,c_orbene(k),Work(iAOMOOvl+ind),iONE
     &                  ,Work(iAOMOOvlE+ind),iONE)
5731      Continue
        Endif
********************************************************************
*
*-- If average natural orbital basis used, transform again. We also
*   define Dim1 and Dim2 as dimensions of matrices depending on
*   whether average natural orbitals have been used. This means that
*   some vectors may be larger than necessary, but since no
*   huge demand of memory is required, this should not cause problem.
*
        If(MoAveRed) then
          Call Dgemm_ ('T','N',nRedMO,iOrb(2),nBaseQ,ONE
     &        ,Work(ipAvRed),nBaseQ,Work(iHalf),nBaseQ,ZERO
     &        ,Work(iTEMP),nRedMO)
          call dcopy_(nDim1*nDim2,Work(iTEMP),iONE,Work(iHalf),iONE)
        Endif
*
*-- Hook on the orbital energy.
*
        call dcopy_(nDimP,[ZERO],iZERO,Work(iHalfE),iONE)
        Do 5751, k=1,iOrb(2)
          ind=nDim1*(k-1)
          Call DaxPy_(nDim1,c_orbene(k),Work(iHalf+ind),iONE
     &              ,Work(iHalfE+ind),iONE)
5751    Continue

*
*-- Construct auxilary matrix for the non-electrostatic operator
*   in AO-basis for QM region. Also construct matrix of pure
*   overlaps for the higher order terms.
*
        Call Dgemm_('N','T',nDim1,nDim1,iOrb(2),ONE
     &      ,Work(iHalf),nDim1,Work(iHalfE),nDim1,ZERO
     &      ,Work(ipAUX),nDim1)
        Call Dgemm_('N','T',nDim1,nDim1,iOrb(2),ONE
     &      ,Work(iHalf),nDim1,Work(iHalf),nDim1,ZERO
     &      ,Work(ipAUXp),nDim1)

*
*-- Accumulate.
*
        Call DaxPy_(nDim1**2,ONE,Work(ipAUX),iONE,Work(ipACC),iONE)
        Call DaxPy_(nDim1**2,ONE,Work(ipAUXp),iONE,Work(ipACCp),iONE)

*Jose*********************************
      If(lExtr(8)) then
        Call Dgemm_('N','T',nBaseQ,nBaseQ,iOrb(2),exrep2
     &            ,Work(iAOMOOvl),nBaseQ,Work(iAOMOOvlE),nBaseQ,ZERO
     &            ,Work(ipAOAUX),nBaseQ)
        Call SqToTri_Q(Work(ipAOAUX),Work(ipAOAUXtri),nBaseQ)
        Call DaxPy_(nGross,ONE,Work(ipAOAUXtri),iONE,Work(ipAOSum)
     &            ,iONE)
      Endif
**************************************
*
*-- The end for this solvent molecule.
*
501   Continue
*
*-- Now construct the matrix elements to the non-electrostatic
*   operator in RASSI basis.
*
      kaunter=0
      facx=0
      Do 5703, iS=1,nState
        Do 5704, jS=1,iS
          HighS=0
          kaunter=kaunter+1
*--------Collect the relevant part of the transistion density matrix.
          Call dCopy_(nDimT,Work(iBigT+nDimT*(kaunter-1))
     &              ,iONE,Work(ipAOG),iONE)
*--------Then transform according to theory.
          Call SqToTri_Q(Work(ipACC),Work(ipACCt),nDim1)
          Addition=Ddot_(nDimT,Work(ipAOG),iONE,Work(ipACCt),iONE)
          SmatRas(kaunter)=SmatRas(kaunter)+exrep2*Addition
*--------And include pure S*S for subsequent higher order overlap
*        repulsion.
          Call SqToTri_Q(Work(ipACCp),Work(ipACCtp),nDim1)
          HighS=Ddot_(nDimT,Work(ipAOG),iONE,Work(ipACCtp),iONE)
          SmatPure(kaunter)=SmatPure(kaunter)+HighS
5704    Continue
5703  Continue

*
*-- Deallocations.
*
      Call GetMem('RotOrb','Free','Real',iV2,nV2size)
      Call GetMem('Sint','Free','Real',ipAOint,nAObaseSize)
      Call GetMem('Sintpar','Free','Real',ipAOintpar,nAObaseSize)
      Call GetMem('HalfTrans','Free','Real',iHalfpar,nHalf)
      Call GetMem('HalfPure','Free','Real',iHalf,nHalf)
      Call GetMem('HalfOrbE','Free','Real',iHalfE,nHalf)
      Call GetMem('Auxilary','Free','Real',ipAux,nBaseQ**2)
      Call GetMem('AuxilaryP','Free','Real',ipAuxp,nBaseQ**2)
      Call GetMem('GammaAO','Free','Real',ipAOG,nGross)
      Call GetMem('Accumulate','Free','Real',ipACC,nBaseQ**2)
      Call GetMem('Accumulate','Free','Real',ipACCp,nBaseQ**2)
      Call GetMem('AccumulateT','Free','Real',ipACCt,nGross)
      Call GetMem('AccumulateTP','Free','Real',ipACCtp,nGross)
      Call GetMem('TEMP','Free','Real',iTEMP,nRedMO*iOrb(2))
*Jose****************************************************************
      If(lExtr(8)) then
         Call GetMem('qAOclMOOvl','Free','Real',iAOMOOvl,nHalf)
         Call GetMem('qAOclMOOvlE','Free','Real',iAOMOOvlE,nHalf)
         Call GetMem('AuxAOp','Free','Real',ipAOAUX,nBaseQ**2)
         Call GetMem('AuxAOpTri','Free','Real',ipAOAUXtri,nGross)
      Endif
*********************************************************************

      Return
      End
