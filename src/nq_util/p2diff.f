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
      Subroutine P2Diff(mGrid,dF_dRho,ndF_dRho)
      use nq_Grid, only: Rho, GradRho
      Implicit Real*8 (A-H,O-Z)
      Dimension dF_dRho(ndF_dRho,mGrid)
#include "nq_info.fh"
#include "real.fh"
#include "nq_index.fh"
      Real*8 MinDns
      Logical GradRho_Allocated
*
      MinDns=0.0001d0

      GradRho_Allocated=Allocated(GradRho)
      If (.NOT.GradRho_Allocated) Then
         Do iGrid=1,mGrid
            If (Abs(Rho(1,iGrid)-Rho(2,iGrid)).le.MinDns  .or.
     &         (Rho(1,iGrid)+Rho(2,iGrid))**2-
     &          4.0d0*Rho(1,iGrid)*Rho(2,iGrid).le. MinDns) Then
               dF_dRho(ipRa,iGrid) = 0.5d0*(dF_dRho(ipRa,iGrid) +
     &                                   dF_dRho(ipRb,iGrid))
               dF_dRho(ipRb,iGrid) = 0.0d0
            Else
               rab1 = 1.0d0/(Rho(1,iGrid)-Rho(2,iGrid))
               dFdRa = dF_dRho(ipRa,iGrid)
               dFdRb = dF_dRho(ipRb,iGrid)
               dF_dRho(ipRb,iGrid) = rab1* (dFdRb-dFdRa)
               dF_dRho(ipRa,iGrid) = rab1*
     &                              (Rho(1,iGrid)*dFdRa-
     &                              Rho(2,iGrid)*dFdRb)
            End If
         End Do   ! iGrid
      Else


         Do iGrid=1,mGrid
            If(abs(Rho(1,iGrid)-Rho(2,iGrid)).le.MinDns.or.
     &        (Rho(1,iGrid)+Rho(2,iGrid))**2-
     &         4.0d0*Rho(1,iGrid)*Rho(2,iGrid).le.
     &         MinDns) Then
               dF_dRho_ = 0.5d0*(
     &              dF_dRho(ipRa,iGrid) + dF_dRho(ipRb,iGrid))
               dF_dRhox = 0.5d0*(
     &              dF_dRho(ipdRxa,iGrid) + dF_dRho(ipdRxb,iGrid))
               dF_dRhoy = 0.5d0*(
     &              dF_dRho(ipdRya,iGrid) + dF_dRho(ipdRyb,iGrid))
               dF_dRhoz = 0.5d0*(
     &              dF_dRho(ipdRza,iGrid) + dF_dRho(ipdRzb,iGrid))
*
               If (.True.) Then
               dF_dRho(ipRa,iGrid) = dF_dRho_
               dF_dRho(ipdRxa,iGrid) = dF_dRhox
               dF_dRho(ipdRya,iGrid) = dF_dRhoy
               dF_dRho(ipdRza,iGrid) = dF_dRhoz
               dF_dRho(ipRb,iGrid) = 0.0d0
               dF_dRho(ipdRxb,iGrid) = 0.0d0
               dF_dRho(ipdRyb,iGrid) = 0.0d0
               dF_dRho(ipdRzb,iGrid) = 0.0d0
               Else
               rho_tot=Rho(1,iGrid)+Rho(2,iGrid)
               dF_dRho(ipRa,iGrid) = 0.0d0
               dF_dRho(ipdRxa,iGrid) = 0.0d0
               dF_dRho(ipdRya,iGrid) = 0.0d0
               dF_dRho(ipdRza,iGrid) = 0.0d0
               dF_dRho(ipRb,iGrid) = 2.0d0*dF_dRho_/rho_tot -
     &         4.0d0*(GradRho(1,iGrid)*dF_dRhox+
     &                GradRho(2,iGrid)*dF_dRhoy+
     &                GradRho(3,iGrid)*dF_dRhoz)/rho_tot/rho_tot
               dF_dRho(ipdRxb,iGrid) = 2.0d0*dF_dRhox/rho_tot
               dF_dRho(ipdRyb,iGrid) = 2.0d0*dF_dRhoy/rho_tot
               dF_dRho(ipdRzb,iGrid) = 2.0d0*dF_dRhoz/rho_tot
               End If

            Else
               rho_tot=Rho(1,iGrid)+Rho(2,iGrid)
               rab1 = 1.0d0/(Rho(1,iGrid)-Rho(2,iGrid))
               dFdRa  = dF_dRho(ipRa,iGrid)
               dFdRb  = dF_dRho(ipRb,iGrid)
               dFdRax = dF_dRho(ipdRxa,iGrid)
               dFdRbx = dF_dRho(ipdRxb,iGrid)
               dFdRay = dF_dRho(ipdRya,iGrid)
               dFdRby = dF_dRho(ipdRyb,iGrid)
               dFdRaz = dF_dRho(ipdRza,iGrid)
               dFdRbz = dF_dRho(ipdRzb,iGrid)
*
              tmp1= (Rho(1,iGrid)*Rho(4,iGrid)-Rho(2,iGrid)*
     &                            GradRho(3,iGrid))*(dFdRbx-dFdRax)+
     &              (Rho(1,iGrid)*Rho(5,iGrid)-Rho(2,iGrid)*
     &                            GradRho(4,iGrid))*(dFdRby-dFdRay)+
     &              (Rho(1,iGrid)*Rho(6,iGrid)-Rho(2,iGrid)*
     &                            GradRho(5,iGrid))*(dFdRbz-dFdRaz)
*
               tmp2=
     &          (GradRho(4,iGrid)-GradRho(1,iGrid))*(dFdRbx-dFdRax)+
     &          (GradRho(5,iGrid)-GradRho(2,iGrid))*(dFdRby-dFdRay)+
     &          (GradRho(6,iGrid)-GradRho(3,iGrid))*(dFdRbz-dFdRaz)
*
               dF_dRho(ipRa,iGrid) = rab1*
     &                             (Rho(1,iGrid)*dFdRa-
     &                              Rho(2,iGrid)*dFdRb)
     &                              -rab1*rab1*tmp1
               dF_dRho(ipdRxa,iGrid) = rab1*
     &                             (Rho(1,iGrid)*dFdRax-
     &                              Rho(2,iGrid)*dFdRbx)
               dF_dRho(ipdRya,iGrid) = rab1*
     &                             (Rho(1,iGrid)*dFdRay-
     &                              Rho(2,iGrid)*dFdRby)
               dF_dRho(ipdRza,iGrid) = rab1*
     &                             (Rho(1,iGrid)*dFdRaz-
     &                              Rho(2,iGrid)*dFdRbz)
*
               dF_dRho(ipRb,iGrid) = rab1*
     &                             (dFdRb-dFdRa) +
     &                              rab1*rab1*tmp2
               dF_dRho(ipdRxb,iGrid) = rab1*
     &                             (dFdRbx-dFdRax)
               dF_dRho(ipdRyb,iGrid) = rab1*
     &                             (dFdRby-dFdRay)
               dF_dRho(ipdRzb,iGrid) = rab1*
     &                             (dFdRbz-dFdRaz)
            End If
         End Do   ! iGrid

      End If
*
      Return
      End
