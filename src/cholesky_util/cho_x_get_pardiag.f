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
* Copyright (C) Francesco Aquilante                                    *
************************************************************************

      SUBROUTINE CHO_X_GET_PARDIAG(jSym,ip_List_rs,iSO_ab)

************************************************************************
*      Returns an array of the "a" and "b" indeces that give rise to the
*      parent diagonal from which a given J-index has been originated
*      by the (molecular) Cholesky decomposition procedure
*
*        ip_List_rs : pointer to the portion of InfVec corresponding
*                     to loc=1 and JSym, thus
*                     ip_of_iWork(InfVec(1,1,jSym)
*
*        ia=iSO_ab(1,numcho(jSym))  contains the index of the basis "a"
*                                   within its symm. block.
*                                   (Note: there is an integer function
*                                   "cho_isao(ia)" that returns the
*                                   irrep to which "a" belongs)
*
*        "Location 2" of iSO_ab returns the same kind of info for "b"
*
*
*  Author: F. Aquilante
*
************************************************************************
#ifdef _MOLCAS_MPP_
      Use Para_Info, Only: MyRank, nProcs
#endif
      use ChoArr, only: iRS2F
      Implicit Real*8 (a-h,o-z)
      Integer ip_List_rs, iSO_ab(2,*)

#include "WrkSpc.fh"
#include "cholesky.fh"
#include "choptr.fh"

      iOff=0

#ifdef _MOLCAS_MPP_
      Call GetMem('List','Allo','Inte',ip_List,nProcs)
      Call IZero(iWork(ip_List),nProcs)
      iWork(ip_List+MyRank)=NumCho(jSym)
      Call Cho_GAIGOP(iWork(ip_List),nProcs,'+')
      nTot=0
      Do iRank=1,nProcs
         nTot = nTot + iWork(ip_List+iRank-1)
         If (iRank.eq.MyRank) iOff = nTot
      End Do
      Call GetMem('List','Free','Inte',ip_List,nProcs)
#endif

      Do jv=1,NumCho(jSym)

         jRab = iWork(ip_List_rs+jv-1) ! addr. 1st red set within jSym

         kRab = iiBstr(jSym,1) + jRab ! global addr. 1st red set

         iSO_ab(1,jv+iOff)=iRS2F(1,kRab) !global addr. of a in SO list
         iSO_ab(2,jv+iOff)=iRS2F(2,kRab) !same for b, always (a .ge. b)
      End Do

#ifdef _MOLCAS_MPP_
      Call GetMem('Aux_iOS_ab','Allo','Inte',ip_Aux,2*nTot)
      Call Izero(iWork(ip_Aux),2*nTot)
      Do jv=1,NumCho(jSym)
         kv=iWork(ip_List_rs+4*MaxVec+jv-1) ! InfVec(jv,5,jSym)
         kAux=ip_Aux+2*(kv-1)
         iWork(kAux) = iSO_ab(1,jv+iOff)
         kAux=kAux+1
         iWork(kAux) = iSO_ab(2,jv+iOff)
      End Do
      Call Icopy(2*nTot,iWork(ip_Aux),1,iSO_ab(1,1),1)
      Call GetMem('Aux_iOS_ab','Free','Inte',ip_Aux,2*nTot)
      Call Cho_GAIGOP(iSO_ab,2*nTot,'+')
#endif
      Return
      End
