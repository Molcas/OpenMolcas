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
* Copyright (C) 1996, Roland Lindh                                     *
************************************************************************
      Subroutine Trsn(xyz,nCent,Tau,Bt,lWrite,lWarn,Label,dBt,ldB)
************************************************************************
*                                                                      *
* Reference: Molecular Vibrations, E. Bright Wilson, Jr, J. C. Decicius*
*            and Paul C. Cross, Sec. 4-1, Eq. 20-24                    *
*                                                                      *
* R.Lindh May-June '96                                                 *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Bt(3,nCent), xyz(3,nCent), dBt(3,nCent,3,nCent),
     &       BRij(3,2), dBRij(3,2,3,2), BRjk(3,2), dBRjk(3,2,3,2),
     &       BRkl(3,2), dBRkl(3,2,3,2), Bf2(3,3),
     &       Bf3(3,3)
      Logical lWrite, lWarn, ldB
      Character*8 Label
      Dimension Dum(1)
*
      mCent=2
      Call Strtch(xyz(1,1),mCent,Rij1,BRij,.False.,Label,dBRij,ldB)
      Call Strtch(xyz(1,2),mCent,Rjk1,BRjk,.False.,Label,dBRjk,ldB)
      Call Strtch(xyz(1,3),mCent,Rkl1,BRkl,.False.,Label,dBRkl,ldB)
      mCent=3
      Call Bend(xyz(1,1),mCent,Fi2,Bf2,.False.,.False.,Label,Dum,
     &          .False.)
      SinFi2=Sin(Fi2)
      CosFi2=Cos(Fi2)
      Call Bend(xyz(1,2),mCent,Fi3,Bf3,.False.,.False.,Label,Dum,
     &          .False.)
      SinFi3=Sin(Fi3)
      CosFi3=Cos(Fi3)
*
      If (SinFi2*SinFi3.lt.1.0d-13) Then
         Tau=Zero
         dTau=Zero
         If (lWrite) Write (6,1) Label,-dTau,-Tau
         Return
      End If
*
*     Get the angle between the two planes, i.e. the
*     angle between the normal vectors.
*
*     r123 * r234 = CosTau
*
      CosTau = ( ( BRij(2,1)*BRjk(3,2) - BRij(3,1)*BRjk(2,2) ) *
     &           ( BRjk(2,1)*BRkl(3,2) - BRjk(3,1)*BRkl(2,2) ) +
     &           ( BRij(3,1)*BRjk(1,2) - BRij(1,1)*BRjk(3,2) ) *
     &           ( BRjk(3,1)*BRkl(1,2) - BRjk(1,1)*BRkl(3,2) ) +
     &           ( BRij(1,1)*BRjk(2,2) - BRij(2,1)*BRjk(1,2) ) *
     &           ( BRjk(1,1)*BRkl(2,2) - BRjk(2,1)*BRkl(1,2) ) )
     &         / (SinFi2*SinFi3)
*
*     For the vector product of the two vectors. This
*     will give a vector parallel to e23. The direction
*     relative to e23 defines the sign.
*
*     e123 X e234 = SinTau * e23
*
      SinTau = ( BRij(1,2) * (BRjk(2,1)*BRkl(3,2)-BRjk(3,1)*BRkl(2,2))
     &         + BRij(2,2) * (BRjk(3,1)*BRkl(1,2)-BRjk(1,1)*BRkl(3,2))
     &         + BRij(3,2) * (BRjk(1,1)*BRkl(2,2)-BRjk(2,1)*BRkl(1,2)) )
     &         / (SinFi2*SinFi3)
*
*     (-Pi < Tau <= Pi)
*
      Tau = ATan2(SinTau,CosTau)
      If (Abs(Tau).eq.Pi) Tau=Pi
*
      dTau = 180.0D+00*Tau/Pi
      dFi2 = 180.0D+00*Fi2/Pi
      dFi3 = 180.0D+00*Fi3/Pi
      If (lWarn) Then
         If (dTau.gt.177.5 .or. dTau.lt.-177.5) Then
            Call WarningMessage(1,
     &                  ' Warning: dihedral angle close to'
     &         //' end of range')
         End If
         If (dFi2.gt.177.5 .or. dFi2.lt.2.5) Then
            Call WarningMessage(1,
     &                  ' Warning: bond angle 2 close to'
     &         //' end of range')
         End If
         If (dFi3.gt.177.5 .or. dFi3.lt.2.5) Then
            Call WarningMessage(1,
     &                  ' Warning: bond angle 3 close to'
     &         //' end of range')
         End If
      End If
      If (LWRITE) Write (6,1) Label,-dTau,-Tau
1     FORMAT(1X,A,' : Dihedral= ',F10.4,
     & '   / Degree  ',F10.6,' / rad')
*
*---- Compute the WDC matrix.
*
      Do ix = 1, 3
         iy=ix+1
         If (iy.gt.3) iy=iy-3
         iz=iy+1
         If (iz.gt.3) iz=iz-3
         Bt(ix,1) = (BRij(iy,2)*BRjk(iz,2)-BRij(iz,2)*BRjk(iy,2))
     &           / (Rij1*SinFi2**2)
         Bt(ix,4) = (BRkl(iy,1)*BRjk(iz,1)-BRkl(iz,1)*BRjk(iy,1))
     &           / (Rkl1*SinFi3**2)
         Bt(ix,2) = -( (Rjk1-Rij1*CosFi2) * Bt(ix,1)
     &             +         Rkl1*CosFi3  * Bt(ix,4))/Rjk1
         Bt(ix,3) = - ( Bt(ix,1)+Bt(ix,2)+Bt(ix,4))
      End Do
*
      If (ldB) Then
*
*------- Compute the derivative of the WDC matrix.
*
         Do ix = 1, 3
            iy=ix+1
            If (iy.gt.3) iy=iy-3
            iz=iy+1
            If (iz.gt.3) iz=iz-3
            Do jx = 1, ix
               jy=jx+1
               If (jy.gt.3) jy=jy-3
               jz=jy+1
               If (jz.gt.3) jz=jz-3
*
            dBt(ix,1,jx,1) =(  dBRij(ix,1,jy,2)*BRjk(jz,2)
     &                       - dBRij(ix,1,jz,2)*BRjk(jy,2)
     &                       - Bt(jx,1)*(BRij(ix,1)*SinFi2**2
     &                       + Rij1*Two*SinFi2*CosFi2*Bf2(ix,1)) )
     &                       / (Rij1*SinFi2**2)
            dBt(ix,1,jx,2) =-( (-BRij(ix,1)*CosFi2
     &                         + Rij1*SinFi2*Bf2(ix,1) ) * Bt(jx,1)
     &                         + (Rjk1-Rij1*CosFi2) * dBt(ix,1,jx,1) )
     &                       / Rjk1
            dBt(jx,2,ix,1) = dBt(ix,1,jx,2)
            dBt(ix,1,jx,4) = Zero
            dBt(jx,4,ix,1) = dBt(ix,1,jx,4)
            dBt(ix,1,jx,3) = - (dBt(ix,1,jx,1) + dBt(ix,1,jx,2) )
            dBt(jx,3,ix,1) = dBt(ix,1,jx,3)
            dBt(ix,4,jx,4) =(  dBRkl(ix,2,jy,1)*BRjk(jz,1)
     &                       - dBRkl(ix,2,jz,1)*BRjk(jy,1)
     &                       - Bt(jx,4)*(BRkl(ix,2)*SinFi3**2
     &                       + Rkl1*Two*SinFi3*CosFi3*Bf3(ix,3)) )
     &                       / (Rkl1*SinFi3**2)
            dBt(ix,4,jx,3) =-( (-BRkl(ix,2)*CosFi3
     &                         + Rkl1*SinFi3*Bf3(ix,3) ) * Bt(jx,4)
     &                         + (Rjk1-Rkl1*CosFi3) * dBt(ix,4,jx,4) )
     &                       / Rjk1
            dBt(jx,3,ix,4) = dBt(ix,4,jx,3)
            dBt(ix,4,jx,2) = - ( dBt(ix,4,jx,4) + dBt(ix,4,jx,3) )
            dBt(jx,2,ix,4) = dBt(ix,4,jx,2)
            If (ix.ne.jx) Then
               dBt(jx,1,ix,1) = dBt(ix,1,jx,1)
               dBt(ix,4,jx,1) = Zero
               dBt(jx,4,ix,4) = dBt(ix,4,jx,4)
               dBt(jx,1,ix,4) = dBt(ix,4,jx,1)
               dBt(jx,1,ix,2) =-( (-BRij(jx,1)*CosFi2
     &                            + Rij1*SinFi2*Bf2(jx,1) ) * Bt(ix,1)
     &                            + (Rjk1-Rij1*CosFi2) * dBt(jx,1,ix,1))
     &                          / Rjk1
               dBt(ix,2,jx,1) = dBt(jx,1,ix,2)
               dBt(ix,3,jx,1) = - ( dBt(ix,1,jx,1) + dBt(ix,2,jx,1)
     &                          + dBt(ix,4,jx,1) )
               dBt(jx,1,ix,3) = dBt(ix,3,jx,1)
               dBt(jx,4,ix,3) =-( (-BRkl(jx,2)*CosFi3
     &                            + Rkl1*SinFi3*Bf3(jx,3) ) * Bt(ix,4)
     &                            + (Rjk1-Rkl1*CosFi3) * dBt(jx,4,ix,4))
     &                          / Rjk1
               dBt(ix,3,jx,4) = dBt(jx,4,ix,3)
               dBt(ix,2,jx,4) = - ( dBt(ix,4,jx,4) + dBt(ix,3,jx,4) )
               dBt(jx,4,ix,2) = dBt(ix,2,jx,4)
            End If
            dBt(ix,2,jx,3) = - ( ( BRjk(ix,1)
     &                           + Rkl1*SinFi3*Bf3(ix,1) ) * Bt(jx,4)
     &                         + ( Rjk1 - Rkl1*CosFi3 ) * dBt(ix,2,jx,4)
     &                         + ( BRij(ix,2)*CosFi2
     &                           - Rij1*SinFi2*Bf2(ix,2) ) * Bt(jx,1)
     &                         +  Rij1*CosFi2 * dBt(ix,2,jx,1)
     &                         + Bt(jx,3) * BRjk(ix,1) ) / Rjk1
            dBt(jx,3,ix,2) = dBt(ix,2,jx,3)
            dBt(ix,2,jx,2) = - ( dBt(ix,2,jx,1) + dBt(ix,2,jx,4)
     &                         + dBt(ix,2,jx,3) )
            dBt(ix,3,jx,3) = - ( dBt(ix,2,jx,3) + dBt(ix,1,jx,3)
     &                         + dBt(ix,4,jx,3) )
            If (ix.ne.jx) Then
               dBt(ix,3,jx,2) = - ( dBt(ix,2,jx,2) + dBt(ix,1,jx,2)
     &                            + dBt(ix,4,jx,2) )
               dBt(jx,2,ix,3) = dBt(ix,3,jx,2)
               dBt(jx,2,ix,2) = dBt(ix,2,jx,2)
               dBt(jx,3,ix,3) = dBt(ix,3,jx,3)
            End If
*
            End Do
         End Do
*
      End If
*
      Fac = -One
      Tau = Fac * Tau
      Call DScal_(3*nCent,Fac,Bt,1)
      If (ldB) Call DScal_((3*nCent)**2,Fac,dBt,1)
*
      Return
      End
