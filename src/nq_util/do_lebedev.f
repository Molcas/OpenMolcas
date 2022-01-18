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
      Subroutine Do_Lebedev(L_Eff,mPt,R)
************************************************************************
*                                                                      *
*     Computes datas useful for the angular quadrature.                *
*                                                                      *
************************************************************************
      use nq_Grid, only: Pax
      Implicit Real*8 (a-h,o-z)
#include "nq_info.fh"
#include "real.fh"
#include "stdalloc.fh"
      Parameter (nSet=11)
      Integer Lebedev_order(nSet), Lebedev_npoints(nSet)
      Data Lebedev_order/5,7,11,17,23,29,35,41,47,53,59/
      Data Lebedev_npoints/14,26,50,110,194,302,434,590,770,974,1202/
      Real*8, Allocatable:: TempR(:,:), TempW(:)
      Real*8, Allocatable:: R(:,:)
*                                                                      *
************************************************************************
*                                                                      *
*---- Generate angular grid a la Lebedev
*
      Do iSet = 1, nSet
         If (Lebedev_order(iSet).eq.L_Eff) Then
            mPt=Lebedev_npoints(iSet)
            Call mma_allocate(R,4,mPt,Label='R')
            Call mma_allocate(TempR,3,mPt,Label='TempR')
            Call mma_allocate(TempW,mPt,Label='TempW')
*
            Call Lebedev(TempR,TempW,nPt,mPt,L_Eff)
            If (nPt.ne.mPt) Then
               Call WarningMessage(2,'Lebedev_Grid: nPt.ne.mPt')
               Write (6,*) 'nPt=',nPt
               Write (6,*) 'mPt=',mPt
               Call Abend()
            End If
            Call DScal_(nPt,Four*Pi,TempW,1)
*
            Call DGEMM_('N','N',
     &                  3,nPt,3,
     &                  1.0d0,Pax,3,
     &                        TempR,3,
     &                  0.0d0,R,4)
            call dcopy_(nPt,TempW,1,R(4,1),4)
*
            Call mma_deallocate(TempW)
            Call mma_deallocate(TempR)
*
            Return
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
      End
