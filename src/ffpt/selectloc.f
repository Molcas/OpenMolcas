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
* Copyright (C) Anders Ohrn                                            *
************************************************************************
      Subroutine SelectLoc(H0,nSize)
************************************************************
*
*   <DOC>
*     <Name>SelectLoc</Name>
*     <Syntax>Call SelectLoc(H0,nSize)</Syntax>
*     <Arguments>
*       \Argument{H0}{The one-electron hamiltonain with perturbations so far. On output the localaized perturbed one-electron hamiltonian.}{}{inout}
*       \Argument{nSize}{Size of the triangular H0 with the additional origo and nuclear contribution.}{}{in}
*     </Arguments>
*     <Purpose>To localize the perturbation LoProp-style. This way a perturbation can be applied selectively on parts of a molecule.</Purpose>
*     <Dependencies></Dependencies>
*     <Author>A.Ohrn</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*      Collect H0 as it is, clean the vacuum part so only perturbation
*      is there. Collect LoProp transformation and transform. Put zeros
*      according to user specification and transform back. The localized
*      perturbation is added to the one-electron hamiltonian and H0 is
*      returned.
*     </Description>
*    </DOC>
*
************************************************************
      Implicit Real*8 (a-h,o-z)

#include "input.fh"
#include "WrkSpc.fh"

      Parameter (ZERO=0.0D+0,ONE=1.0D+0)
      Real*8 H0(nSize)
      Character*8 Label
      Logical Debug,OneOrNot1,OneOrNot2,OneOrNot3,OneOrNot4,OneOrNot
      Logical CrazySet
      Data Debug /.false./

*
*-- Commence!
*
      Call qEnter('SELECTLOC')
      Write(6,*)
      Write(6,*)' The perturbation will be localized "LoProp style".'
      Write(6,*)
      Write(6,*)' -- Number of basis subsets:',nSets
      Do k=1,nSets
        Write(6,*)'    ',iSelection(1,k),iSelection(2,k)
      Enddo
      Write(6,*)' -- Atoms and bonds logical flags:'
      Do i=1,nSets
        Write(6,*)'      Set atom:  ',i,Atoms(i)
        Do j=i+1,nSets
          Write(6,*)'      Sets bond: ',i,j,Bonds(i,j)
        Enddo
      Enddo

*
*-- No symmetry allowed.
*
      Call Get_iScalar('nSym',nSym)
      If(nSym.ne.1) then
        Write(6,*)
      Write(6,*)' You have specified symmetry. The keyword "SELEctive"',
     & ' in FFPT is incompatible with this.'
        Write(6,*)' Aborting....'
        Call Abend()
      Endif

*
*-- Collect the overlap and some auxiliaries for LoProp.
*
      Call GetMem('Orbital_Type','Allo','Inte',ipType,nBas(1))
      Call GetMem('Center_Index','Allo','Inte',ipCenter,nBas(1))
      Call Get_iArray('Orbital Type',iWork(ipType),nBas(1))
      Call Get_iArray('Center Index',iWork(ipCenter),nBas(1))
      Do i=ipType,ipType+nBas(1)-1
        If(iWork(i).ne.1.and.iWork(i).ne.0) then
          Write(6,*)'Orbital type vector is corrupted!'
          Call Abend()
        Endif
      Enddo

      iOpt2=2
      iOpt1=1
      iOpt0=0
      Label='MltPl  0'
      iRc=-1
      iSymLbl=1
      Call iRdOne(iRc,iOpt1,Label,1,nInts,iSymLbl)
      Call GetMem('SMatTr','Allo','Real',ipSTr,nInts+4)
      Call RdOne(iRc,iOpt0,Label,1,Work(ipSTr),iSymLbl)
      If(iRc.ne.0) then
        Write(6,*)'Error reading overlap matrix in SELECTLOC!'
        Call QTrace
        Call Abend()
      Endif
*-- Lets be square.
      Call GetMem('SMatSq','Allo','Real',ipSSq,nBas(1)**2)
      Call Square(Work(ipSTr),Work(ipSSq),1,nBas(1),nBas(1))

*
*-- Call the localization utility and get the transformation matrix.
*
      Call GetMem('T','Allo','Real',ipT,nBas(1)**2)
      Call GetMem('Tinv','Allo','Real',ipTinv,nBas(1)**2)
      Call Localize_LoProp(Work(ipT),Work(ipTinv),nBas(1),Work(ipSSq)
     &                    ,iWork(ipCenter),iWork(ipType))
      If(DeBug) then
        Call RecPrt('Total transfMat',' ',Work(ipT),nBas(1),nBas(1))
      Endif

*
*-- Transform the perturbation to the LoProp basis. FFPT accumulates the
*   perturbation to H0, but we only want the perturbation V, hence first
*   a subtraction is necessary.
*
      Call GetMem('VacH0','Allo','Real',iHVac,nInts+4)
      If(LCumulate) Then
         Label='OneHam  '
      Else
         Label='OneHam 0'
      EndIf
      iRc=-1
      Call RdOne(iRc,iOpt2,Label,1,Work(iHVac),iSymLbl)
      If(iRc.ne.0) then
        Write(6,*)'Error reading H0 in SELECTLOC!'
        Call QTrace
        Call Abend()
      Endif
      Call GetMem('Pert','Allo','Real',iV,nInts)
      Do i=1,nInts
        Work(iV+i-1)=H0(i)-Work(iHVac+i-1)
      Enddo
      H01=H0(nInts+1)
      H02=H0(nInts+2)
      H03=H0(nInts+3)
      H04=H0(nInts+4)
*----But first translate the perturbation origo to the relevant centre
      Call TransNow(iV,ipSTr)
*----You may proceed.
      Call GetMem('PertSq','Allo','Real',iVS,nBas(1)**2)
      Call Square(Work(iV),Work(iVS),1,nBas(1),nBas(1))
      If(DeBug) then
        Call RecPrt('Pert:(Basis:ord)',' ',Work(iVS),nBas(1),nBas(1))
      Endif

      Call GetMem('TEMP','Allo','Real',iTEMP,nBas(1)**2)
      Call GetMem('PertL','Allo','Real',iVLoP,nBas(1)**2)
*----Go to basis where overlap matrix, S, is diagonal.
      Call DGEMM_('T','N',nBas(1),nBas(1),nBas(1),ONE,Work(ipT)
     &          ,nBas(1),Work(iVS),nBas(1),ZERO,Work(iTEMP),nBas(1))
      Call DGEMM_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(iTEMP),
     &    nBas(1),Work(ipT),nBas(1),ZERO,Work(iVLoP),nBas(1))
      If(DeBug) then
        Call RecPrt('Pert:(Basis:LoP)',' ',Work(iVLop),nBas(1),nBas(1))
      Endif

*
*-- Set elements to zero as designated in input. The routine below is
*   far from optimal, but we are not in need of great speed here so....
*

      kaunter=0
      Do i=1,nBas(1)
        Do j=1,nBas(1)
          Siff=0.0D+0
          Do k=1,nSets
            Do l=k,nSets
              If(k.ne.l) then
                If(.not.Bonds(k,l)) Goto 999
              Else
                If(.not.Atoms(k)) Goto 999
              Endif
              OneOrNot1=i.ge.iSelection(1,k).and.i.le.iSelection(2,k)
              OneOrNot2=j.ge.iSelection(1,l).and.j.le.iSelection(2,l)
              OneOrNot3=i.ge.iSelection(1,l).and.i.le.iSelection(2,l)
              OneOrNot4=j.ge.iSelection(1,k).and.j.le.iSelection(2,k)
              OneOrNot=(OneOrNot1.and.OneOrNot2).or.
     &                 (OneOrNot3.and.OneOrNot4)
              CrazySet=(OneOrNot1.and.OneOrNot2).and.
     &                 (OneOrNot3.and.OneOrNot4)
              If(CrazySet.and.Bonds(k,l).and..not.Atoms(k)) then
                Write(6,*)'Your set selection is not exclusive!'
              Endif
              If(OneOrNot) Siff=1.0D+0
*-- Here we enable to set the weight in the bond-domain to some
*   other number than one.
              If(OneOrNot.and.(Atoms(k).and.Bonds(k,l))) Siff=SiffBond
999           Continue
            Enddo
          Enddo
          Work(iVLop+kaunter)=Work(iVLop+kaunter)*Siff
          kaunter=kaunter+1
        Enddo
      Enddo

*
*-- Transform back. Due to the non-unitarian and non-orthogonal basis
*   the inverse is is contravariant (if the transformation was
*   covariant). See the Book by Lanczos.
*
      Call DGEMM_('T','N',nBas(1),nBas(1),nBas(1),ONE,Work(ipTinv)
     &          ,nBas(1),Work(iVLop),nBas(1),ZERO,Work(iTEMP),nBas(1))
      Call DGEMM_('N','N',nBas(1),nBas(1),nBas(1),ONE,Work(iTEMP),
     & nBas(1),Work(ipTinv),nBas(1),ZERO,Work(iVS),nBas(1))
      If(DeBug) then
        Call RecPrt('Pert.Zeroed',' ',Work(iVS),nBas(1),nBas(1))
      Endif

*
*-- Add this perturbation to the one-electron hamiltonian.
*
      kaunter=0
      Do i=1,nBas(1)
        Do j=1,i
          kaunter=kaunter+1
          ind=(i-1)*nBas(1)+j
          H0(kaunter)=Work(iHVac+kaunter-1)+Work(iVS+ind-1)
        Enddo
      Enddo
*----And don't forget the orgio and the nuclear repulsion.
      H0(nInts+1)=H01
      H0(nInts+2)=H02
      H0(nInts+3)=H03
      H0(nInts+4)=H04

*
*-- Deallocations en masse.
*
      Call GetMem('Orbital_Type','Free','Inte',ipType,nBas(1))
      Call GetMem('Center_Index','Free','Inte',ipCenter,nBas(1))
      Call GetMem('SMatSq','Free','Real',ipSSq,nBas(1)**2)
      Call GetMem('T','Free','Real',ipT,nBas(1)**2)
      Call GetMem('Tinv','Free','Real',ipTinv,nBas(1)**2)
      Call GetMem('VacH0','Free','Real',iHVac,nInts+4)
      Call GetMem('Pert','Free','Real',iV,nInts)
      Call GetMem('PertSq','Free','Real',iVS,nBas(1)**2)
      Call GetMem('TEMP','Free','Real',iTEMP,nBas(1)**2)
      Call GetMem('PertL','Free','Real',iVLoP,nBas(1)**2)
      Call GetMem('SMatTr','Free','Real',ipSTr,nInts+4)

*
*-- Exit
*
      Write(6,*)
      Write(6,*)'  ....Done!'
      Write(6,*)
      Call qExit('SELECTLOC')
      Return
      End
