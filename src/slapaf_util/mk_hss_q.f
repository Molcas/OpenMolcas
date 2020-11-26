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
      Subroutine Mk_Hss_Q()
      use Slapaf_Info, only: Cx, Coor, Shift
      Implicit Real*8 (a-h,o-z)
#include "info_slapaf.fh"
#include "real.fh"
#include "WrkSpc.fh"

*     Compute the Hessian in internal coordinates.
*
      If ((lNmHss.or.lRowH).and.iter.eq.NmIter) Then
         Call Put_dArray('Unique Coordinates',Cx,3*nsAtom)
         Call Put_Coord_New(Cx,nsAtom)
         If (lRowH) Then
            If (BSet.and.HSet) Call Hss_q()
            Call RowHessian(NmIter,mInt,nRowH,mRowH,Delta/2.5d0,
     &                      Work(ipdqInt))
         Else
            Call FormNumHess(iter,Work(ipdqInt),Shift,mInt,Delta,
     &                       Stop,Work(ipqInt),nsAtom,Cubic,iNeg,
     &                       Work(ipDipM),mTROld,Smmtrc,Degen,UserT,
     &                       UserP,nUserPT,nsRot,lTherm,lDoubleIso,
     &                       Curvilinear)
         End If
*
         call dcopy_(3*nsAtom,Cx,1,Coor,1)
         Call Get_dArray('BMxOld',Work(ipB),3*nsAtom*mInt)
         ipIn = ipqInt + (Iter-1)*mInt
         call dcopy_(mInt,Work(ipqInt),1,Work(ipIn),1)
         ipIn = ipdqInt + (Iter-1)*mInt
         call dcopy_(mInt,Work(ipdqInt),1,Work(ipIn),1)
      Else
         If (BSet.and.HSet) Call Hss_q()
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
