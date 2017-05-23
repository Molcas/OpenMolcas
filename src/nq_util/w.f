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
      Subroutine W(R,ilist_p,Weights,list_p,nlist_p,nGrid)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "itmax.fh"
#include "nq_info.fh"
#include "debug.fh"
      Real*8 R(3,nGrid), Weights(nGrid)
      Integer list_p(nlist_p)
*                                                                      *
************************************************************************
*                                                                      *
#include "nq_structure.fh"
      declare_ip_coor
      p(x)=(x*0.5D0)*(3.0D0-x**2)
*                                                                      *
************************************************************************
*                                                                      *
      P_i = Zero ! dummy initialize
*
*     iNQ is the index of the current atomic grid to which these grid
*     points belong.
*
      iNQ=list_p(ilist_p)
C     Write (*,*) 'ilist_p=',ilist_p
C     Write (*,*) 'nlist_p=',nlist_p
C     Write (*,*) 'nGrid=',nGrid
C     Write (*,*) 'iNQ=',iNQ
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, nGrid
*        Write (*,*) 'iGrid=',iGrid
*                                                                      *
************************************************************************
*                                                                      *
*----    Becke's partitioning
*
         Sum_P_k=Zero
         Do klist_p = 1, nlist_p
            kNQ=list_p(klist_p)
*           Write (*,*) 'kNQ=',kNQ
            kx=ip_Coor(kNQ)
            ky=kx+1
            kz=ky+1
            r_k=sqrt((R(1,iGrid)-Work(kx))**2
     &               +(R(2,iGrid)-Work(ky))**2
     &               +(R(3,iGrid)-Work(kz))**2)
            P_k=One
            Do llist_p = 1, nlist_p
               lNQ=list_p(llist_p)
*
               If (kNQ.ne.lNQ) Then
                  lx=ip_Coor(lNQ)
                  ly=lx+1
                  lz=ly+1
*
                  r_l=sqrt((R(1,iGrid)-Work(lx))**2
     &                     +(R(2,iGrid)-Work(ly))**2
     &                     +(R(3,iGrid)-Work(lz))**2)
                  R_kl=sqrt((Work(kx)-Work(lx))**2
     &                      +(Work(ky)-Work(ly))**2
     &                      +(Work(kz)-Work(lz))**2)
                  rMU_kl=(r_k-r_l)/R_kl
                  If (rMU_kl.le.0.5D0) Then
                     p1=p(rMU_kl)
                     p2=p(p1)
                     p3=p(p2)
                     s=Half*(One-p3)
                  Else
                     xdiff=rMU_kl-1.0D0
                     xdiff=(-1.5D0-0.5D0*xdiff)*xdiff**2
                     xdiff=(-1.5D0-0.5D0*xdiff)*xdiff**2
                     p3=   ( 1.5D0+0.5D0*xdiff)*xdiff**2
                     s=Half*p3
                  End If
                  P_k=P_k*s
               End If
            End Do
*
            If (kNQ.eq.iNQ) P_i=P_k
            Sum_P_k = Sum_P_k + P_k
         End Do
         Fact=Weights(iGrid)
         Weights(iGrid)=Fact*P_i/Sum_P_k
C        Write (*,*) 'Fact,P_A,Z,Weights=',Fact,P_i,Sum_P_k,
C    &               Weights(iGrid)
*                                                                      *
************************************************************************
*                                                                      *
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
