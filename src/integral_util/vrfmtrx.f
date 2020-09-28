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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      SubRoutine VrfMtrx(Label,lOper,nComp,ip,Matrix)
************************************************************************
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, Sweden, January '91                  *
************************************************************************
      use Basis_Info, only: nBas
      use Temporary_Parameters, only: PrPrt
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
*     Local arrays
      Real*8 Matrix(*)
      Character Label*(*), Line*80, Status
      Integer ip(nComp), lOper(nComp)
*
      Call GetEnvf('MOLCAS_TEST_not_yet_here',Status)
      If (Status.eq.' ') Return
      Do 10 iComp = 1, nComp
         VrfSum=0.0D0
         ip1 = ip(iComp)
         iSmLbl = lOper(iComp)
         If (Prprt) iSmLbl = iAnd(1,iSmLbl)
         Do 30 iIrrep = 0, nIrrep - 1
            If (nBas(iIrrep).le.0) Go To 30
            Do 40 jIrrep = 0, iIrrep
               If (nBas(jIrrep).le.0) Go To 40
               If (iAnd(iSmLbl,2**iEor(iIrrep,jIrrep)).eq.0) Go To 40
               If (iIrrep.eq.jIrrep) Then
                  n2=nBas(iIrrep)*(nBas(iIrrep)+1)/2
                  VrfSum=VrfSum+DDot_(n2,Matrix(ip1),1,Matrix(ip1),1)
                  ip1 = ip1 + n2
               Else
                  n2 =  nBas(iIrrep)*nBas(jIrrep)
                  VrfSum=VrfSum+DDot_(n2,Matrix(ip1),1,Matrix(ip1),1)
                  ip1 = ip1 + n2
               End If
 40         Continue
 30      Continue
*        Add the nuclear contribution and the operator position.
         n2=4
         VrfSum=VrfSum+DDot_(n2,Matrix(ip1),1,Matrix(ip1),1)
         Write (Line,'(A,I5)') Label,iComp
         Call Add_info(Line,[VrfSum],1,8)
 10   Continue
*
      Return
      End
