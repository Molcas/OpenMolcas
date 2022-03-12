!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Extract(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,HMatOld  &
     &                  ,xyzQuQ,ip_ExpVal,ip_ExpCento,ENR,ENP)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "qminp.fh"

      Dimension xyzMy(3),Hmat(*),HMatOld(*),xyzQuQ(6)
      Dimension iDt(3)

!
!---  Just pass on the numbers according to QM-method.
!
      If(QmType(1:4).eq.'RASS') then
        Call ExtractR(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,HMatOld     &
     &               ,xyzQuQ,lExtr,iExtr_Eig,iExtr_Atm,ip_ExpVal        &
     &               ,ip_ExpCento,ENR,ENP)
      ElseIf(QmType(1:3).eq.'SCF') then
        Call ExtractS(iLu,i9,Etot,xyzMy,Hmat,iC,iDt,nMatBas,xyzQuQ      &
     &               ,lExtr,iExtr_Atm,ip_ExpVal,ip_ExpCento,ENR,ENP)
      Endif

      Return
      End
