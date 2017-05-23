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
      Subroutine W_Partitioning(R,ilist_p,Fact,Weights,list_p,nlist_p,
     &                          dW_dR,nGrad_Eff,iTab,Do_Grad,dW_Temp,
     &                          dPB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "itmax.fh"
#include "nq_info.fh"
#include "debug.fh"
      Real*8 R(3), Weights, dW_dR(nGrad_Eff), dW_Temp(3,nlist_p),
     &       dPB(3,nlist_p,nlist_p), sxyz(3)
*
      Integer list_p(nlist_p), iTab(4,nGrad_Eff)
      Logical Do_Grad
*                                                                      *
************************************************************************
*                                                                      *
#include "nq_structure.fh"
      declare_ip_coor
      p(x)=(x/Two)*(Three-x**2)
*                                                                      *
************************************************************************
*                                                                      *
ct    R(1)=Zero
ct    R(2)=Zero
ct    R(3)=Zero
ct    Fact=One
      P_i = Zero ! dummy initialize
*
*     iNQ is the index of the current atomic grid to which these grid
*     points belong.
*
      iNQ=list_p(ilist_p)
*                                                                      *
************************************************************************
*                                                                      *
      If (.Not.(Do_Grad.and.Grid_Type.eq.Moving_Grid)) Then
*                                                                      *
************************************************************************
*                                                                      *
*----    Becke's partitioning
*
         Sum_P_k=Zero
*        P_k=Zero
         Do klist_p = 1, nlist_p
            kNQ=list_p(klist_p)
            kx=ip_Coor(kNQ)
            ky=kx+1
            kz=ky+1
            r_k=sqrt((R(1)-Work(kx))**2
     &               +(R(2)-Work(ky))**2
     &               +(R(3)-Work(kz))**2)
            P_k=One
            Do llist_p = 1, nlist_p
               lNQ=list_p(llist_p)
               lx=ip_Coor(lNQ)
               ly=lx+1
               lz=ly+1
               If (kNQ.ne.lNQ) Then
                  r_l=sqrt((R(1)-Work(lx))**2
     &                     +(R(2)-Work(ly))**2
     &                     +(R(3)-Work(lz))**2)
                  R_kl=sqrt((Work(kx)-Work(lx))**2
     &                      +(Work(ky)-Work(ly))**2
     &                      +(Work(kz)-Work(lz))**2)
                  rMU_kl=(r_k-r_l)/R_kl
c                 p1=p(rMU_kl)
c                 p2=p(p1)
c                 p3=p(p2)
                  p3=p(p(p(rMU_kl)))
                  s=Half*(One-p3)
                  P_k=P_k*s
               End If
            End Do
            If (kNQ.eq.iNQ) P_i=P_k
            Sum_P_k = Sum_P_k + P_k
         End Do
         Weights=Fact*P_i/Sum_P_k
*                                                                      *
************************************************************************
*                                                                      *
      Else
*                                                                      *
************************************************************************
*                                                                      *
         Call FZero(dW_Temp,3*nlist_p)
*
*------- The current grid point is associated with center A and the
*        "atomic" displacement vector relative to center A is computed.
*
         iA=ilist_p
         sxyz(1) = R(1)-Work(ip_Coor(iNQ)  )
         sxyz(2) = R(2)-Work(ip_Coor(iNQ)+1)
         sxyz(3) = R(3)-Work(ip_Coor(iNQ)+2)
*
         Z=Zero
         Call FZero(dPB,3*nlist_p**2)
*
*        Compute all P_B and corresponding derivatives.
*
         P_A=Zero
         Do iB = 1, nlist_p
            kNQ=list_p(iB)
            kx=ip_Coor(kNQ)
            r_Bx=R(1)-Work(kx  )
            r_By=R(2)-Work(kx+1)
            r_Bz=R(3)-Work(kx+2)
            r_B=sqrt(r_Bx**2+r_By**2+r_Bz**2)
*
*           loop over C=/=B for all s(mu_BC), see Eq. B3
*
            P_B=One
            Do iC = 1, nlist_p
*
               If (iC.ne.iB) Then
                  lNQ=list_p(iC)
                  lx=ip_Coor(lNQ)
*
                  r_Cx=R(1)-Work(lx  )
                  r_Cy=R(2)-Work(lx+1)
                  r_Cz=R(3)-Work(lx+2)
                  r_C=sqrt(r_Cx**2+r_Cy**2+r_Cz**2)
                  R_BCx=Work(kx  )-Work(lx  )
                  R_BCy=Work(kx+1)-Work(lx+1)
                  R_BCz=Work(kx+2)-Work(lx+2)
                  R_BC=sqrt(R_BCx**2+R_BCy**2+R_BCz**2)
*
*                 Eq. B6
*
                  rMU_BC=(r_B-r_C)/R_BC
*
*                 Eqs. B5
*
                  p1=p(rMU_BC)
                  p2=p(p1)
                  p3=p(p2)
*
*                 Eq. B4
*
                  s_MU_BC=Half*(One-p3)
                  If (s_MU_BC.eq.Zero) s_MU_BC=1.0D-99
*
*                 Eq. B10
*
                  P_B=P_B*s_MU_BC
                  tMU_BC=-27D0*(One-p2**2)*(One-p1**2)*(One-rMU_BC**2)
     &                  /(16d0*s_MU_BC)
*
*---------------- Differentiate mu_BC with respect to the center, D.
*
                  Do iD = 1, nlist_p
*
                     If (iD.eq.iB) Then
*
*                       d mu_BC(r_A) / dB, Eq. B10
*
                        dmu_BC_dBx =  -r_Bx/(r_B*R_BC)
     &                             - (r_B-r_C)*R_BCx/R_BC**3
                        dmu_BC_dBy =  -r_By/(r_B*R_BC)
     &                             - (r_B-r_C)*R_BCy/R_BC**3
                        dmu_BC_dBz =  -r_Bz/(r_B*R_BC)
     &                             - (r_B-r_C)*R_BCz/R_BC**3
*
                        dPB(1,iD,iB)=dPB(1,iD,iB) + tMU_BC*dmu_BC_dBx
                        dPB(2,iD,iB)=dPB(2,iD,iB) + tMU_BC*dmu_BC_dBy
                        dPB(3,iD,iB)=dPB(3,iD,iB) + tMU_BC*dmu_BC_dBz
*
                     Else If (iD.eq.iC) Then
*
*                       d mu_BC(r_A) / dC, Eq, B10
*
                        dmu_BC_dCx =  r_Cx/(r_C*R_BC)
     &                             + (r_B-r_C)*R_BCx/R_BC**3
                        dmu_BC_dCy =  r_Cy/(r_C*R_BC)
     &                             + (r_B-r_C)*R_BCy/R_BC**3
                        dmu_BC_dCz =  r_Cz/(r_C*R_BC)
     &                             + (r_B-r_C)*R_BCz/R_BC**3
                        dPB(1,iD,iB)=dPB(1,iD,iB) + tMU_BC*dmu_BC_dCx
                        dPB(2,iD,iB)=dPB(2,iD,iB) + tMU_BC*dmu_BC_dCy
                        dPB(3,iD,iB)=dPB(3,iD,iB) + tMU_BC*dmu_BC_dCz
*
                     End If
*
*                    d mu_BC(r_A) / dr_A
*
                     dmu_BC_dAx = (r_Bx/r_B - r_Cx/r_C)/R_BC
                     dmu_BC_dAy = (r_By/r_B - r_Cy/r_C)/R_BC
                     dmu_BC_dAz = (r_Bz/r_B - r_Cz/r_C)/R_BC
*
                     If (iD.eq.iA) Then
*
*                       The direct term
*
                        dPB(1,iD,iB)=dPB(1,iD,iB) + tMU_BC*dmu_BC_dAx
                        dPB(2,iD,iB)=dPB(2,iD,iB) + tMU_BC*dmu_BC_dAy
                        dPB(3,iD,iB)=dPB(3,iD,iB) + tMU_BC*dmu_BC_dAz
                     End If
*
                  End Do  ! iD
*
               End If
            End Do        ! iC
*
*           Multiply derivatives with P_B as in Eq. B8
            Do iD = 1, nlist_p
               dPB(1,iD,iB)=P_B*dPB(1,iD,iB)
               dPB(2,iD,iB)=P_B*dPB(2,iD,iB)
               dPB(3,iD,iB)=P_B*dPB(3,iD,iB)
            End Do
*
            If (iB.eq.iA) P_A=P_B
*
*           Denominator Eq. B2
            Z = Z + P_B
         End Do           ! iB
*
*        Eq. B2
*
         Weights=Fact*P_A/Z
*                                                                      *
************************************************************************
*                                                                      *
*------- Assemble the gradient
*
         If (Debug) Call RecPrt('dPB',' ',dPB,3*nlist_p,nlist_p)
         Do iB = 1, nlist_p
*
            dZ_dBx=Zero
            dZ_dBy=Zero
            dZ_dBz=Zero
            Do iC = 1, nlist_p
               dZ_dBx=dZ_dBx+dPB(1,iB,iC)
               dZ_dBy=dZ_dBy+dPB(2,iB,iC)
               dZ_dBz=dZ_dBz+dPB(3,iB,iC)
            End Do
*
*           Eq. B7
*
            dW_Temp(1,iB) = dPB(1,iB,iA)/Z - (P_A*dZ_dBx)/Z**2
            dW_Temp(2,iB) = dPB(2,iB,iA)/Z - (P_A*dZ_dBy)/Z**2
            dW_Temp(3,iB) = dPB(3,iB,iA)/Z - (P_A*dZ_dBz)/Z**2
         End Do
         If (Debug) Call RecPrt('dW_Temp',' ',dW_Temp,3,nlist_p)
*                                                                      *
************************************************************************
*                                                                      *
*        Pick up the relevant gradients
*
         Call FZero(dW_dR,nGrad_Eff)
         Do iGrad = 1, nGrad_Eff
            iCar=iTab(1,iGrad)
            kNQ =iTab(3,iGrad)
            Fact0=Fact*DBLE(Abs(iTab(4,iGrad)))
            Do iB = 1, nlist_p
               lNQ=list_p(iB)
               If (kNQ.eq.lNQ) dW_dR(iGrad)=Fact0*dW_Temp(iCar,iB)
            End Do
         End Do
         If (Debug) Call RecPrt('dW_dR',' ',dW_dR,1,nGrad_Eff)
*                                                                      *
************************************************************************
*                                                                      *
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
