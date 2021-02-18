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
#ifdef _MOLCAS_MPP_
#include "stdalloc.fh"
      Integer, Allocatable:: ip_List(:), ip_Aux(:,:)
#endif

      iOff=0

#ifdef _MOLCAS_MPP_
      Call mma_allocate(ip_List,[0,nProcs-1],Label='ip_List')
      ip_List(:)=0
      ip_List(MyRank)=NumCho(jSym)
      Call Cho_GAIGOP(ip_List,nProcs,'+')
      nTot=0
      Do iRank=1,nProcs
         nTot = nTot + ip_List(iRank-1)
         If (iRank.eq.MyRank) iOff = nTot
      End Do
      Call mma_deallocate(ip_List)
#endif

      Do jv=1,NumCho(jSym)

         jRab = iWork(ip_List_rs+jv-1) ! addr. 1st red set within jSym

         kRab = iiBstr(jSym,1) + jRab ! global addr. 1st red set

         iSO_ab(1,jv+iOff)=iRS2F(1,kRab) !global addr. of a in SO list
         iSO_ab(2,jv+iOff)=iRS2F(2,kRab) !same for b, always (a .ge. b)
      End Do

#ifdef _MOLCAS_MPP_
      Call mma_allocate(ip_Aux,2,nTot,Label='ip_Aux')
      ip_Aux(:,:)=0
      Do jv=1,NumCho(jSym)
         kv=iWork(ip_List_rs+4*MaxVec+jv-1) ! InfVec(jv,5,jSym)
         ip_Aux(1,kv) = iSO_ab(1,jv+iOff)
         ip_Aux(2,kv) = iSO_ab(2,jv+iOff)
      End Do
      iSO_ab(:,1:nTot)=ip_Aux(:,1:nTot)
      Call mma_deallocate(ip_Aux)
      Call Cho_GAIGOP(iSO_ab,2*nTot,'+')
#endif
      Return
      End
