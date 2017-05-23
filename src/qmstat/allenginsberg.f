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
      Subroutine AllenGinsberg(QMMethod,Eint,Poli,dNuc,Cha,Dip,Qua
     &                        ,MxBaux,iVEC,nDim,iExtr_Atm,lEig,iEig
     &                        ,iQ_Atoms,ip_ExpCento,E_Nuc_Part
     &                        ,lSlater,Eint_Nuc)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.fh"

      Dimension Eint(MxQCen,10),Poli(MxQCen,10)
      Dimension Cha(MxBaux,MxQCen),Dip(MxBaux,3,MxQCen)
      Dimension Qua(MxBaux,6,MxQCen),dNuc(MxAt)
      Dimension iExtr_Atm(MxAt),Eint_Nuc(MxAt)
      Dimension iCenSet(MxAt**2)
      Character QMMethod*5
      Logical lEig,Check1,Check2,lSlater

*
*-- Set up centre index set. The order of centres are decided by
*   the MpProp-program and are hence collected in the get_center
*   routine.
*
*
*-- Atom centres
*
      Do 901, iAt=1,MxAt
        If(iExtr_Atm(iAt).eq.-1) then
          GoTo 902
        Else
          iCenSet(iAt)=iExtr_Atm(iAt)
        Endif
901   Continue
      NExtrAt=iAt
902   Continue
      NExtrAt=iAt-1
*
*-- Bond centres
*
      kaunter=iQ_Atoms
      kaunt=NExtrAt
      Do 905, iAt=2,iQ_Atoms
        Do 906, jAt=1,iAt-1
          kaunter=kaunter+1
          Check1=.false.
          Check2=.false.
          Do 907, i1=1,NExtrAt
            If(iAt.eq.iCenSet(i1)) Check1=.true.
            If(jAt.eq.iCenSet(i1)) Check2=.true.
907       Continue
          If(Check1.and.Check2) then
            kaunt=kaunt+1
            iCenSet(kaunt)=kaunter
          Endif
906     Continue
905   Continue
*
*-- A minor check.
*
      NExpect=NExtrAt*(nExtrAt+1)/2
      NTotal=kaunt
      If(NTotal.ne.NExpect) then
        Write(6,*)
        Write(6,*)' Error in atom specification for partial'
     &//' perturbation extraction.'
        Call Quit(_RC_GENERAL_ERROR_)
      Endif

*
*-- Compute partial nuclear contribution.
*
      E_Nuc_Part=0.0d0
      Do 931, iAt=1,NExtrAt
        iCx=iCenSet(iAt)
        If(lSlater) then
          E_Nuc_Part=E_Nuc_Part
     &         -(Eint_Nuc(iCx)+Poli(iCx,1))*dNuc(iCx)
        Else
          E_Nuc_Part=E_Nuc_Part
     &         -(Eint(iCx,1)+Poli(iCx,1))*dNuc(iCx)
        Endif
931   Continue

*
*-- Set up matrix elements for the partial perturbations.
*   Compare with hel, helstate, polink and polins.
*
      Call GetMem('VelPart','Allo','Real',iVelP,nDim*(nDim+1)/2)
      Call GetMem('VpoPart','Allo','Real',iVpoP,nDim*(nDim+1)/2)
      call dcopy_(nDim*(nDim+1)/2,ZERO,iZERO,Work(iVelP),iONE)
      call dcopy_(nDim*(nDim+1)/2,ZERO,iZERO,Work(iVpoP),iONE)
      kk=0
      Do 911, i=1,nDim
        Do 912, j=1,i
          kk=kk+1
          Do 913, k=1,NTotal
            iCx=iCenSet(k)
            dMp=Cha(kk,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,1)*dMp
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,1)*dMp
            dMp=Dip(kk,1,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,2)*dMp
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,2)*dMp
            dMp=Dip(kk,2,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,3)*dMp
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,3)*dMp
            dMp=Dip(kk,3,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,4)*dMp
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,4)*dMp
            dMp=Qua(kk,1,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,5)*dMp
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,5)*dMp
            dMp=Qua(kk,3,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,7)*dMp
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,7)*dMp
            dMp=Qua(kk,6,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,10)*dMp
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,10)*dMp
            dMp=Qua(kk,2,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,6)*dMp*2.0d0
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,6)*dMp*2.0d0
            dMp=Qua(kk,4,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,8)*dMp*2.0d0
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,8)*dMp*2.0d0
            dMp=Qua(kk,5,iCx)
            Work(iVelP+kk-1)=Work(iVelP+kk-1)+Eint(iCx,9)*dMp*2.0d0
            Work(iVpoP+kk-1)=Work(iVpoP+kk-1)+Poli(iCx,9)*dMp*2.0d0
913       Continue
912     Continue
911   Continue

*
*-- Collect expectation value for the partial perturbation.
*
      Call Expectus(QMMethod,Work(iVelP),Work(iVelP),Work(iVpoP)
     &             ,Work(iVpoP),MxBaux,iVEC,nDim,lEig,iEig
     &             ,ip_ExpCento)

*
*-- Deallocate.
*
      Call GetMem('VelPart','Free','Real',iVelP,nDim*(nDim+1)/2)
      Call GetMem('VpoPart','Free','Real',iVpoP,nDim*(nDim+1)/2)

*
*-- Howl
*

      Return
      End
