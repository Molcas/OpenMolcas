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
      Subroutine Do_Lebedev(L_Eff,mPt,ipR)
************************************************************************
*                                                                      *
*     Computes datas useful for the angular quadrature.                *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "nq_info.fh"
#include "real.fh"
#include "WrkSpc.fh"
      Parameter (nSet=11)
      Integer Lebedev_order(nSet), Lebedev_npoints(nSet)
      Data Lebedev_order/5,7,11,17,23,29,35,41,47,53,59/
      Data Lebedev_npoints/14,26,50,110,194,302,434,590,770,974,1202/
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate angular grid a la Lebedev
*
      Do iSet = 1, nSet
         If (Lebedev_order(iSet).eq.L_Eff) Then
            mPt=Lebedev_npoints(iSet)
            Call GetMem('AngRW','Allo','Real',ipR,4*mPt)
            Call GetMem('temp','Allo','Real',iptemp,4*mPt)
*
            iOffR=iptemp
            iOffW=iOffR+3*mPt
            Call Lebedev(Work(iOffR),Work(iOffW),nPt,mPt,L_Eff)
            If (nPt.ne.mPt) Then
               Call WarningMessage(2,'Lebedev_Grid: nPt.ne.mPt')
               Write (6,*) 'nPt=',nPt
               Write (6,*) 'mPt=',mPt
               Call Abend()
            End If
            Call DScal_(nPt,Four*Pi,Work(iOffW),1)
*
            Call DGEMM_('N','N',
     &                  3,nPt,3,
     &                  1.0d0,Work(ip_O),3,
     &                  Work(iOffR),3,
     &                  0.0d0,Work(ipR),4)
            call dcopy_(nPt,Work(iOffW  ),1,Work(ipR+3),4)
*
            Call GetMem('temp','Free','Real',iptemp,4*nPt)
*
            Go To 99
*
         End If
      End Do
      Write (6,'(A,I3)') 'Failed to find a Lebedev grid of order', L_EFF
      Write (6,'(A)') 'Available orders are:'
      Write (6,'(11(1X,I3))') (Lebedev_order(i),i=1,nSet)
      Call Abend()
*                                                                      *
************************************************************************
*                                                                      *
 99   Return
      End
