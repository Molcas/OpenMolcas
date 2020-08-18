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
      Subroutine SOOUT(label,cnt_ico,phase_ico)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "info.fh"
#include "lundio.fh"
      Integer cnt_ico(0:7,*),phase_ico(0:7,*)
      Character Label(MaxBfn+MaxBfn_Aux)*(LENIN8)
      Character ChOper(0:7)*3
      Data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
*
      Call SOCtl_mod(ChOper,Label,Maxbfn+MaxBfn_Aux,0,0,Cnt_ico,
     &               Phase_ico)
*
      Return
      End
*
      Subroutine SOCtl_mod(ChOper,Mamn,nMamn,nDkroll,nCall,Cnt_ico,
     &                     Phase_ico)
      use Basis_Info
      Implicit Real*8 (a-h,o-z)
*
#include "itmax.fh"
#include "info.fh"

#include "WrkSpc.fh"
#include "real.fh"
#include "print.fh"
*
      Character ChOper(0:7)*3, ChTemp*8, Mamn(nMamn)*(LENIN8)
      Logical kECP, TstFnc
      Integer cnt_ico(0:7,*),phase_ico(0:7,*)
*
*     Generate list of symmetry adapted or petite list basis functions
*
*     Loop over Irreps
      iSO = 0
*
*     Loop over irreducible representations and symmetry operations,
*     respectively, for SO and Petite list, respectively.
*
      Do 200 iIrrep = 0, nIrrep-1
*
*        Loop over distinct shell types
*
         mdc = 0
         mc  = 1
         iShell = 0
         Do 201 iCnttp = 1, nCnttp
            kECP = dbsc(iCnttp)%ECP
            If (AuxCnttp(iCnttp).or.dbsc(iCnttp)%Frag) Go To 201
*
*           Loop over distinct centers
*
            Do 202 iCnt = 1, dbsc(iCnttp)%nCntr
               mdc = mdc + 1
*
*              Loop over shells associated with this center
*              Start with s type shells
*
               kComp = 0
               iSh = ipVal(iCnttp) - 1
               Do 203 iAng = 0, nVal_Shells(iCnttp)-1
                  iSh = iSh + 1
                  iShell = iShell + 1
                  nExpi=Shells(iSh)%nExp
                  If (nExpi.eq.0) Go To 2033
                  nBasisi=Shells(iSh)%nBasis
                  If (nBasisi.eq.0) Go To 2033
                  jComp = (iAng+1)*(iAng+2)/2
                  If(Shells(iSh)%Prjct ) jComp = 2*iAng + 1
                  Do 204 iComp = 1, jComp
                     lComp = kComp + iComp
*                    Get character of basis function
                     iChBs = iChBas(lComp)
                     If (Shells(iSh)%Transf) iChBs=iChBas(iSphCr(lComp))
*
*                    Skip if function not a basis of irreps.
*
                     If (.Not.TstFnc(iOper,nIrrep,iCoSet(0,0,mdc),
     &                   nIrrep/nStab(mdc),iChTbl,iIrrep,iChBs,
     &                   nStab(mdc))) Go To 204
                     IrrCmp(IndS(iShell)+iComp) =
     &                    iOr(IrrCmp(IndS(iShell)+iComp),2**iIrrep)
*
                     Do 205 iCntrc = 1, nBasisi
                        iSO = iSO + 1
                        If (iSO.gt.nMamn) Then
                           Call WarningMessage(2,'SOout: iSO.gt.nMamn')
                           Call Abend()
                        End If
                        ChTemp=LblCBs(lComp)
                        If (Shells(iSh)%Transf) ChTemp=LblSbs(lComp)
                        Do ico=0,nIrrep/nStab(mdc)-1
                        Cnt_ico(ico,iso)=mc+ico
                        Phase_ico(ico,iso)=
     &                        iPrmt(NrOpr(iCoSet(iCo,0,mdc),
     &                        iOper,nIrrep),iChbs)*
     &                        iChTbl(iIrrep,NrOpr(iCoSet(iCo,0,mdc),
     &                        iOper,nIrrep))
                        End Do
                        Mamn(iSO)=LblCnt(mdc)(1:LENIN)//ChTemp(1:8)
 205                 Continue
*
 204              Continue
 2033             continue
                  kComp = kComp + (iAng+1)*(iAng+2)/2
 203           Continue
               mc = mc + nIrrep/nStab(mdc)
 202        Continue
*
 201     Continue
 200  Continue
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_character(ChOper)
         Call Unused_integer(nDkroll)
         Call Unused_integer(nCall)
      End If
      End
