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
      Subroutine Merge_Lists(Mode,nAt)
      Implicit Real*8 (a-h,o-z)
#include "stdalloc.fh"
      Character*1 Mode
      Logical Found
      Real*8, Dimension(:,:), Allocatable :: rList
      Integer, Dimension(:,:), Allocatable :: iList
*                                                                      *
************************************************************************
*                                                                      *
*     Get the sizes of the arrays on the runfile. Same size for both
*     files.
*
      Call qpg_iArray('Slapaf Info 1',Found,n1)
      Call qpg_dArray('Slapaf Info 2',Found,n2)
*                                                                      *
************************************************************************
*                                                                      *
*     Loop over both files and pick up the fields and pick up pointers,
*     and interations counts.
*
      ip_Dummy=-1
      i_R  = 0
      ipEner_R= ip_Dummy
      ipCx_R  = ip_Dummy
      ipGx_R  = ip_Dummy
      i_P  = 0
      ipEner_P= ip_Dummy
      ipCx_P  = ip_Dummy
      ipGx_P  = ip_Dummy
*
      Call mma_allocate(iList,n1,2,label='iList')
      Call mma_allocate(rList,n2,2,label='rList')
      Call NameRun('RUNREAC')
      Call Get_iArray('Slapaf Info 1',iList(1,1),n1)
      Call Get_dArray('Slapaf Info 2',rList(1,1),n2)
      i_R=1
      iter_R  =     iList(2,1)
      ipEner_R= 1 + iList(5,1)
      ipCx_R  = 1 + iList(6,1)
      ipGx_R  = 1 + iList(7,1)
      Call NameRun('RUNPROD')
      Call Get_iArray('Slapaf Info 1',iList(1,2),n1)
      Call Get_dArray('Slapaf Info 2',rList(1,2),n2)
      i_P=2
      iter_P  =     iList(2,2)
      ipEner_P= 1 + iList(5,2)
      ipCx_P  = 1 + iList(6,2)
      ipGx_P  = 1 + iList(7,2)
*                                                                      *
************************************************************************
*                                                                      *
*     Translate into generic variables, _1=from, _2=to.

      If (Mode.eq.'R') Then
         i_1 = i_P
         iter_1=iter_P
         ipEner_1=ipEner_P
         ipCx_1=ipCx_P
         ipGx_1=ipGx_P
*
         i_2 = i_R
         iter_2=iter_R
         ipEner_2=ipEner_R
         ipCx_2=ipCx_R
         ipGx_2=ipGx_R
      Else
         i_1 = i_R
         iter_1=iter_R
         ipEner_1=ipEner_R
         ipCx_1=ipCx_R
         ipGx_1=ipGx_R
*
         i_2 = i_P
         iter_2=iter_P
         ipEner_2=ipEner_P
         ipCx_2=ipCx_P
         ipGx_2=ipGx_P
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Start moving stuff around
*
*     Update iteration counter
*
      iter_3=iter_2+1
      iList(2,i_2)=iter_3
      iOff_1=(iter_1-1)*3*nAt
      iOff_2=(iter_2-1)*3*nAt
      iOff_3=(iter_3-1)*3*nAt
*
*     Move the last item(s) in the "to" file up one step.
*
      rList(ipEner_2+iter_3-1,i_2)=rList(ipEner_2+iter_2-1,i_2)
      call dcopy_(3*nAt,rList(ipCx_2+iOff_2,i_2),1,
     &                  rlist(ipCx_2+iOff_3,i_2),1)
      call dcopy_(3*nAt,rList(ipGx_2+iOff_2,i_2),1,
     &                  rlist(ipGx_2+iOff_3,i_2),1)
*
*     Copy the last item(s) in the "from" file into the
*     second last position of the "to" file.
*
      rList(ipEner_2+iter_2-1,i_2)=rList(ipEner_1+iter_1-1,i_1)
      call dcopy_(3*nAt,rList(ipCx_1+iOff_1,i_1),1,
     &                  rlist(ipCx_2+iOff_2,i_2),1)
      call dcopy_(3*nAt,rList(ipGx_1+iOff_1,i_1),1,
     &                  rlist(ipGx_2+iOff_2,i_2),1)
*                                                                      *
************************************************************************
*                                                                      *
*     Write out the stuff on the appropriate run file.
*
      If (Mode.eq.'R') Then
         Call NameRun('RUNREAC')
      Else
         Call NameRun('RUNPROD')
      End If
*
      Call Put_iArray('Slapaf Info 1',iList(1,i_2),n1)
      Call Put_dArray('Slapaf Info 2',rList(1,i_2),n2)
      Call qpg_iScalar('iOff_Iter',Found)
      If (Found) Then
         Call Get_iScalar('iOff_Iter',iOff_Iter)
         Call Put_iScalar('iOff_Iter',iOff_Iter+1)
      EndIf
*                                                                      *
************************************************************************
*                                                                      *
*     Deallocate the memory allocations.
*
      Call mma_deallocate(rList)
      Call mma_deallocate(iList)
*                                                                      *
************************************************************************
*                                                                      *
*     Open the default run file.
*
      Call NameRun('RUNFILE')
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
