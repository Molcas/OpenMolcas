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
      Subroutine dWdR(R,ilist_p,Weights,list_p,nlist_p,
     &                dW_dR,nGrad_Eff,iTab,dW_Temp,dPB,nGrid)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "itmax.fh"
#include "nq_info.fh"
#include "debug.fh"
      Real*8 R(3,nGrid), Weights(nGrid), dW_dR(nGrad_Eff,nGrid),
     &       dW_Temp(3,nlist_p), dPB(3,nlist_p,nlist_p), sxyz(3),
     &       dOdxs(3), Osxyz(3)
      Integer list_p(nlist_p), iTab(4,nGrad_Eff)
*                                                                      *
************************************************************************
*                                                                      *
#include "nq_structure.fh"
      declare_ip_coor
      declare_ip_dodx
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
      iA=ilist_p
      O11=Work(ip_O)
      O12=Work(ip_O+1)
      O13=Work(ip_O+2)
      O21=Work(ip_O+3)
      O22=Work(ip_O+4)
      O23=Work(ip_O+5)
      O31=Work(ip_O+6)
      O32=Work(ip_O+7)
      O33=Work(ip_O+8)
*                                                                      *
************************************************************************
*                                                                      *
      Do iGrid = 1, nGrid
         Call FZero(dW_dR(1,iGrid),nGrad_Eff)
*                                                                      *
************************************************************************
*                                                                      *
*------- The current grid point is associated with center A and the
*        "atomic" displacement vector relative to center A is computed.
*
         Osxyz(1) = R(1,iGrid)-Work(ip_Coor(iNQ)  )
         Osxyz(2) = R(2,iGrid)-Work(ip_Coor(iNQ)+1)
         Osxyz(3) = R(3,iGrid)-Work(ip_Coor(iNQ)+2)
*
         sxyz(1)=O11*Osxyz(1)+O12*Osxyz(2)+O13*Osxyz(3)
         sxyz(2)=O21*Osxyz(1)+O22*Osxyz(2)+O23*Osxyz(3)
         sxyz(3)=O31*Osxyz(1)+O32*Osxyz(2)+O33*Osxyz(3)
*
         Z=Zero
         Call FZero(dPB,3*nlist_p**2)
*
*        Compute all P_B and corresponding derivatives.
*
C        P_A=Zero
         P_A=One ! Dummy set
         Do iiB = 1, nlist_p
            If (iiB.eq.1) Then
               iB = iA
            Else If (iiB.eq.iA) Then
               iB = 1
            Else
               iB = iiB
            End If
*
            kNQ=list_p(iB)
            kx=ip_Coor(kNQ)
            r_Bx=R(1,iGrid)-Work(kx  )
            r_By=R(2,iGrid)-Work(kx+1)
            r_Bz=R(3,iGrid)-Work(kx+2)
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
                  r_Cx=R(1,iGrid)-Work(lx  )
                  r_Cy=R(2,iGrid)-Work(lx+1)
                  r_Cz=R(3,iGrid)-Work(lx+2)
                  r_C=sqrt(r_Cx**2+r_Cy**2+r_Cz**2)
                  R_BCx=Work(kx  )-Work(lx  )
                  R_BCy=Work(kx+1)-Work(lx+1)
                  R_BCz=Work(kx+2)-Work(lx+2)
                  R_BC=sqrt(R_BCx**2+R_BCy**2+R_BCz**2)
*
*                 Eq. B6
*
                  rMU_BC=(r_B-r_C)/R_BC
                  If (rMU_BC.le.0.5D0) Then
                     p1=p(rMU_BC)
                     p2=p(p1)
                     p3=p(p2)
*
*                 Eq. B4
*
                     s_MU_BC=Half*(One-p3)
*
                     P_B=P_B*s_MU_BC
                     If (P_B.le.1.0D-20) Go To 99
                     tMU_BC=-27D0*(One-p2**2)
     &                           *(One-p1**2)
     &                           *(One-rMU_BC**2)
     &                     /(16d0*Max(s_MU_BC,1.0D-99))
                  Else
                     xdiff0=rMU_BC-1.0D0
                     xdiff1=(-1.5D0-0.5D0*xdiff0)*xdiff0**2
                     xdiff2=(-1.5D0-0.5D0*xdiff1)*xdiff1**2
                     p3=    ( 1.5D0+0.5D0*xdiff2)*xdiff2**2
                     s_MU_BC=Half*p3
*
                     P_B=P_B*s_MU_BC
                     If (P_B.le.1.0D-20) Go To 99
                     tMU_BC= 27D0*(2.0D0+xdiff2)*xdiff2
     &                           *(2.0D0+xdiff1)*xdiff1
     &                           *(2.0D0+xdiff0)*xdiff0
     &                     /(16d0*Max(s_MU_BC,1.0D-99))
                  End If
*
*---------------- Differentiate mu_BC with respect to the center, D.
*
                  Do iD = 1, nlist_p
C                    jNQ=list_p(iD)
*
                     If (iD.eq.iB) Then
*
*                       d mu_BC(r_A) / dB, Eq. B10
*
                        If (r_B.eq.Zero) Then
                           dmu_BC_dBx =
     &                                - (r_B-r_C)*R_BCx/R_BC**3
                           dmu_BC_dBy =
     &                                - (r_B-r_C)*R_BCy/R_BC**3
                           dmu_BC_dBz =
     &                                - (r_B-r_C)*R_BCz/R_BC**3
                        Else
                           dmu_BC_dBx =  -r_Bx/(r_B*R_BC)
     &                                - (r_B-r_C)*R_BCx/R_BC**3
                           dmu_BC_dBy =  -r_By/(r_B*R_BC)
     &                                - (r_B-r_C)*R_BCy/R_BC**3
                           dmu_BC_dBz =  -r_Bz/(r_B*R_BC)
     &                                - (r_B-r_C)*R_BCz/R_BC**3
                        End If
*
                        dPB(1,iB,iB)=dPB(1,iB,iB) + tMU_BC*dmu_BC_dBx
                        dPB(2,iB,iB)=dPB(2,iB,iB) + tMU_BC*dmu_BC_dBy
                        dPB(3,iB,iB)=dPB(3,iB,iB) + tMU_BC*dmu_BC_dBz
*
                     Else If (iD.eq.iC) Then
*
*                       d mu_BC(r_A) / dC, Eq, B10
*
                        If (r_C.eq.Zero) Then
                           dmu_BC_dCx =
     &                                + (r_B-r_C)*R_BCx/R_BC**3
                           dmu_BC_dCy =
     &                                + (r_B-r_C)*R_BCy/R_BC**3
                           dmu_BC_dCz =
     &                                + (r_B-r_C)*R_BCz/R_BC**3
                        Else
                           dmu_BC_dCx =  r_Cx/(r_C*R_BC)
     &                                + (r_B-r_C)*R_BCx/R_BC**3
                           dmu_BC_dCy =  r_Cy/(r_C*R_BC)
     &                                + (r_B-r_C)*R_BCy/R_BC**3
                           dmu_BC_dCz =  r_Cz/(r_C*R_BC)
     &                                + (r_B-r_C)*R_BCz/R_BC**3
                        End If
                        dPB(1,iC,iB)=dPB(1,iC,iB) + tMU_BC*dmu_BC_dCx
                        dPB(2,iC,iB)=dPB(2,iC,iB) + tMU_BC*dmu_BC_dCy
                        dPB(3,iC,iB)=dPB(3,iC,iB) + tMU_BC*dmu_BC_dCz
*
                     End If
*
*                    d mu_BC(r_A) / dr_A
*
                     If (r_B.eq.Zero) Then
                        dmu_BC_dAx = (         - r_Cx/r_C)/R_BC
                        dmu_BC_dAy = (         - r_Cy/r_C)/R_BC
                        dmu_BC_dAz = (         - r_Cz/r_C)/R_BC
                     Else If (r_C.eq.Zero) Then
                        dmu_BC_dAx = (r_Bx/r_B           )/R_BC
                        dmu_BC_dAy = (r_By/r_B           )/R_BC
                        dmu_BC_dAz = (r_Bz/r_B           )/R_BC
                     Else
                        dmu_BC_dAx = (r_Bx/r_B - r_Cx/r_C)/R_BC
                        dmu_BC_dAy = (r_By/r_B - r_Cy/r_C)/R_BC
                        dmu_BC_dAz = (r_Bz/r_B - r_Cz/r_C)/R_BC
                     End If
*
                     If (iD.eq.iA) Then
*
*                       The direct term
*
                        dPB(1,iA,iB)=dPB(1,iA,iB) + tMU_BC*dmu_BC_dAx
                        dPB(2,iA,iB)=dPB(2,iA,iB) + tMU_BC*dmu_BC_dAy
                        dPB(3,iA,iB)=dPB(3,iA,iB) + tMU_BC*dmu_BC_dAz
*
                     End If
*
                     jNQ=list_p(iD)
                     Do iCar = 1, 3
C                       Call xxDGeMul(Work(ip_dOdx(jNQ,iCar)),3,'N',
C    &                              sxyz,3,'N',
C    &                              dOdxs,3,
C    &                              3,3,1)
                        dOdx_11=Work(ip_dOdx(jNQ,iCar))
                        dOdx_21=Work(ip_dOdx(jNQ,iCar)+1)
                        dOdx_31=Work(ip_dOdx(jNQ,iCar)+2)
                        dOdx_12=Work(ip_dOdx(jNQ,iCar)+3)
                        dOdx_22=Work(ip_dOdx(jNQ,iCar)+4)
                        dOdx_32=Work(ip_dOdx(jNQ,iCar)+5)
                        dOdx_13=Work(ip_dOdx(jNQ,iCar)+6)
                        dOdx_23=Work(ip_dOdx(jNQ,iCar)+7)
                        dOdx_33=Work(ip_dOdx(jNQ,iCar)+8)
                        dOdxs(1)=dOdx_11*sxyz(1)
     &                          +dOdx_12*sxyz(2)
     &                          +dOdx_13*sxyz(3)
                        dOdxs(2)=dOdx_21*sxyz(1)
     &                          +dOdx_22*sxyz(2)
     &                          +dOdx_23*sxyz(3)
                        dOdxs(3)=dOdx_31*sxyz(1)
     &                          +dOdx_32*sxyz(2)
     &                          +dOdx_33*sxyz(3)
                        temp = tMU_BC*(dmu_BC_dAx*dOdxs(1)
     &                                +dmu_BC_dAy*dOdxs(2)
     &                                +dmu_BC_dAz*dOdxs(3))
*
                        dPB(iCar,iD,iB) = dPB(iCar,iD,iB) - temp
                     End Do
*
                  End Do  ! iD
*
               End If
            End Do        ! iC
*
 99         Continue
*
*           Multiply derivatives with P_B as in Eq. B8
*
            Do iD = 1, nlist_p
               dPB(1,iD,iB)=P_B*dPB(1,iD,iB)
               dPB(2,iD,iB)=P_B*dPB(2,iD,iB)
               dPB(3,iD,iB)=P_B*dPB(3,iD,iB)
            End Do
*
            If (iB.eq.iA) P_A=P_B
            If (P_A.le.1.0D-20) Go To 98
*
*           Denominator Eq. B2
            Z = Z + P_B
         End Do           ! iB
*
         If (P_A.eq.Zero) Then
            Fact=Zero
         Else
            Fact = Weights(iGrid)*Z/P_A
         End If
*                                                                      *
************************************************************************
*                                                                      *
*------- Assemble the gradient
*
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
*                                                                      *
************************************************************************
*                                                                      *
*        Pick up the relevant gradients
*
         Do iGrad = 1, nGrad_Eff
            iCar=iTab(1,iGrad)
            kNQ =iTab(3,iGrad)
            Fact0=Fact
            Do iB = 1, nlist_p
               lNQ=list_p(iB)
               If (kNQ.eq.lNQ) dW_dR(iGrad,iGrid)=Fact0*dW_Temp(iCar,iB)
            End Do
         End Do
*                                                                      *
************************************************************************
*                                                                      *
 98      Continue
      End Do ! iGrid
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
