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
      Subroutine Extract(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,HMatOld
     &                  ,xyzQuQ,ip_ExpVal,ip_ExpCento,ENR,ENP)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"

      Dimension xyzMy(3),Hmat(*),HMatOld(*),xyzQuQ(6)
      Dimension iDt(3)

*
*---  Just pass on the numbers according to QM-method.
*
      If(QmType(1:4).eq.'RASS') then
        Call ExtractR(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,HMatOld
     &               ,xyzQuQ,lExtr,iExtr_Eig,iExtr_Atm,ip_ExpVal
     &               ,ip_ExpCento,ENR,ENP)
      ElseIf(QmType(1:3).eq.'SCF') then
        Call ExtractS(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,xyzQuQ
     &               ,lExtr,iExtr_Atm,ip_ExpVal,ip_ExpCento,ENR,ENP)
      Endif

      Return
      End


      Subroutine ExtractR(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nState,HMatOld
     &                   ,xyzQu,lExtr,iExtr_Eig,iExtr_Atm,ip_ExpVal
     &                   ,ip_ExpCento,ENR,ENP)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"

      Dimension xyzMy(3),Hmat(*),HMatOld(*),xyzQu(6)
      Dimension iDt(3),iExtr_Atm(*)
      Logical lExtr(*)

      Write(iLu,*)'<<<<<<<Configuration ',i9,'>>>>>>>'
      If(lExtr(1)) then
        Write(iLu,*)'Total Energy'
        Write(iLu,'(F15.8)')Etot
      Endif
      If(lExtr(2)) then
        Write(iLu,*)'QM-Dipole'
        Write(iLu,'(3(F12.5))')(xyzMy(k),k=1,3)
      Endif
      If(lExtr(3)) then
        Write(iLu,*)'QM-Quadrupole'
        Write(iLu,'(6(F12.5))')(xyzQu(k),k=1,6)
      Endif
      If(lExtr(4)) then
        Write(iLu,*)'Eigenvalues of RASSI-matrix'
        Do 1, i=1,iExtr_Eig
          ind=i*(i+1)/2
          Write(iLu,'(F15.8)')Hmat(ind)
1       Continue
      Endif
      If(lExtr(5)) then
        Write(iLu,*)'Corresponding eigenvectors'
        Do 2, j=0,iExtr_Eig-1
          Write(iLu,'(5(F15.8))')(Work(iC+j*nState+k),k=0,nState-1)
2       Continue
      Endif
      If(lExtr(6)) then
        Write(iLu,*)'Expectation values (H_0,V_el,V_pol,V_pp)'
        Write(iLu,*)'  Nuc cont:',ENR
        If(lExtr(4)) nDim=iExtr_Eig
        If(.not.lExtr(4)) nDim=nState
        Do 3, i=1,nDim
          Write(iLu,'(4(F15.8))')(Work(ip_ExpVal+4*(i-1)+k),k=0,3)
3       Continue
        Call GetMem('ExpVals','Free','Real',ip_ExpVal,4*nDim)
      Endif
      If(lExtr(7)) then
        Write(iLu,*)'Expectation values partial V_el, V_pol'
        Write(iLu,*)'  Nuc cont:',ENP
        If(lExtr(4)) nDim=iExtr_Eig
        If(.not.lExtr(4)) nDim=nState
        Do 4, j=1,nDim
          Write(iLu,'(2(F15.8))')(Work(ip_ExpCento+4*(j-1)+k),k=1,2)
4       Continue
        Call GetMem('ExpVals','Free','Real',ip_ExpCento,4*nDim)
      Endif

      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer_array(iDt)
        Call Unused_real_array(HMatOld)
        Call Unused_integer_array(iExtr_Atm)
      End If
      End


      Subroutine ExtractS(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nBas,xyzQu
     &                   ,lExtr,iExtr_Atm,ip_ExpVal,ip_ExpCento
     &                   ,ENR,ENP)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"

      Dimension xyzMy(3),Hmat(*),xyzQu(6)
      Dimension iDt(3),iExtr_Atm(*)
      Logical lExtr(*)

      Write(iLu,*)'<<<<<<<Configuration ',i9,'>>>>>>>'
      If(lExtr(1)) then
        Write(iLu,*)'Total Energy'
        Write(iLu,'(F15.8)')Etot
      Endif
      If(lExtr(2)) then
        Write(iLu,*)'QM-Dipole'
        Write(iLu,'(3(F12.5))')(xyzMy(k),k=1,3)
      Endif
      If(lExtr(3)) then
        Write(iLu,*)'QM-Quadrupole'
        Write(iLu,'(6(F12.5))')(xyzQu(k),k=1,6)
      Endif
      If(lExtr(6)) then
        Write(iLu,*)'Expectation values (T+H_nuc,V_el,V_pol,V_pp)'
        Write(iLu,*)'  Nuc cont:',ENR
        Write(iLu,'(4(F15.8))')(Work(ip_ExpVal+k),k=0,3)
        Call GetMem('ExpVals','Free','Real',ip_ExpVal,4)
      Endif
      If(lExtr(7)) then
        Write(iLu,*)'Expectation values partial V_el, V_pol'
        Write(iLu,*)'  Nuc cont:',ENP
        Write(iLu,'(2(F15.8))')(Work(ip_ExpCento+k),k=1,2)
        Call GetMem('ExpVals','Free','Real',ip_ExpCento,4)
      Endif

      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_real_array(Hmat)
        Call Unused_integer(iC)
        Call Unused_integer_array(iDt)
        Call Unused_integer(nBas)
        Call Unused_integer_array(iExtr_Atm)
      End If
      End
