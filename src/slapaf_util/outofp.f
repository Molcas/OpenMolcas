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
      Subroutine OutofP(xyz,nCent,Teta,Bt,lWrite,lWarn,Label,dBt,ldB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Bt(3,nCent), xyz(3,nCent), dBt(3,nCent,3,nCent),
     &       R42(3),      R43(3),
     &       C14X(3,3),     BR14X(3,3), dBR14X(3,3,3,3)
      Logical lWrite, lWarn, ldB
      Character*8 Label
*#define _Test_Numerical_
#ifdef _Test_Numerical_
      Real*8 Bt_temp_p(3,4)
      Real*8 Bt_temp_m(3,4)
#endif
*
      Call qEnter('OutofP')
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     First some diagnostics
*
*     4-->1 (Bond)
      RX1=xyz(1,1)-xyz(1,4)
      RY1=xyz(2,1)-xyz(2,4)
      RZ1=xyz(3,1)-xyz(3,4)
      R41KV=RX1*RX1+RY1*RY1+RZ1*RZ1
      Q41=SQRT(R41KV)
      e41x = RX1 / Q41
      e41y = RY1 / Q41
      e41z = RZ1 / Q41
*     4-->2 (Bond in plane)
      RX2=xyz(1,2)-xyz(1,4)
      RY2=xyz(2,2)-xyz(2,4)
      RZ2=xyz(3,2)-xyz(3,4)
      R42KV=RX2*RX2+RY2*RY2+RZ2*RZ2
      Q42=SQRT(R42KV)
      e42x = RX2 / Q42
      e42y = RY2 / Q42
      e42z = RZ2 / Q42
*     4-->3 (Bond in plane)
      RX3=xyz(1,3)-xyz(1,4)
      RY3=xyz(2,3)-xyz(2,4)
      RZ3=xyz(3,3)-xyz(3,4)
      R43KV=RX3*RX3+RY3*RY3+RZ3*RZ3
      Q43=SQRT(R43KV)
      e43x = RX3 / Q43
      e43y = RY3 / Q43
      e43z = RZ3 / Q43
*
*     Get the angle between e43 and e42
*
      CosFi1 = e43x*e42x + e43y*e42y + e43z*e42z
*
      Fi1=ArCos(CosFi1)
      If (Abs(CosFi1).gt.One) Call RecPrt('xyz(1)',' ',xyz,3,4)
      dFi1 = 180.D0 * Fi1 / Pi
      If (lWarn.and.(dFi1.gt.177.5D0 .or. dFi1.lt.2.5D0)) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*
*     Dirty exit! This happens when an earlier structure is ill defined.
*
      If (Abs(Fi1-Pi).lt.1.0D-13) Then
         Teta=0.0D0
         Call FZero(Bt,3*nCent)
         Call qExit('OutofP')
         Return
      End If
*
*     Get the angle between e41 and e43
*
      CosFi2 = e41x*e43x + e41y*e43y + e41z*e43z
*
      Fi2=ArCos(CosFi2)
      If (Abs(CosFi2).gt.One) Call RecPrt('xyz(2)',' ',xyz,3,4)
      dFi2 = 180.D0 * Fi2 / Pi
      If (lWarn.and.(dFi2.gt.177.5D0 .or. dFi2.lt.2.5D0)) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*
*     Get the angle between e41 and e42
*
      CosFi3 = e41x*e42x + e41y*e42y + e41z*e42z
*
      Fi3=ArCos(CosFi3)
      If (Abs(CosFi3).gt.One) Call RecPrt('xyz(3)',' ',xyz,3,4)
      dFi3 = 180.D0 * Fi3 / Pi
      If (lWarn.and.(dFi3.gt.177.5D0 .or. dFi3.lt.2.5D0)) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call RecPrt('xyz',' ',xyz,3,nCent)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     The first two centers are trivially
*
      call dcopy_(3,xyz(1,1),1,C14X(1,1),1)
      call dcopy_(3,xyz(1,4),1,C14X(1,2),1)
*
*     The 3rd is
*
      R42(1)=xyz(1,2)-xyz(1,4)
      R42(2)=xyz(2,2)-xyz(2,4)
      R42(3)=xyz(3,2)-xyz(3,4)
      R43(1)=xyz(1,3)-xyz(1,4)
      R43(2)=xyz(2,3)-xyz(2,4)
      R43(3)=xyz(3,3)-xyz(3,4)
      C14X(1,3)=R42(2)*R43(3)-R42(3)*R43(2)
      C14X(2,3)=R42(3)*R43(1)-R42(1)*R43(3)
      C14X(3,3)=R42(1)*R43(2)-R42(2)*R43(1)
*
*     Exit if 2-3-4 are collinear
*     (equivalent to the above check, but this is more concrete)
*
      If ((C14X(1,3)**2+C14X(2,3)**2+C14X(3,3)**2).lt.1.0D-10) Then
         Teta=0.0D0
         Call FZero(Bt,3*nCent)
         Call qExit('OutofP')
         Return
      End If
      C14X(1,3)=C14X(1,3)+xyz(1,4)
      C14X(2,3)=C14X(2,3)+xyz(2,4)
      C14X(3,3)=C14X(3,3)+xyz(3,4)
*
      mCent=3
      Call Bend(C14X,mCent,Teta,BR14X,.False.,.False.,Label,dBR14X,
     &          ldB)
*
      Teta = Teta - Pi/Two
      dTeta= 180.D0 * Teta/ Pi
      If (lWarn.and.(dTeta.gt.87.5D0 .or. dTeta.lt.-87.5D0)) Then
         Write (6,*) 'Warning: Out of plane angle close to'
     &        //' end of range'
      End If
      If (LWRITE) Write (6,'(1X,A,A,F10.4,A,F10.4,A)') Label,
     &   ' : Out of plane angle=',dTeta,'/degree, ',Teta,'/rad'
*
*---- Compute the WDC matrix
*
      Do ix = 1, 3
         iy = Mod(ix+1, 4)+(ix+1)/4
         iz = Mod(iy+1, 4)+(iy+1)/4
*
         Bt(ix,1) = - BR14X(ix,1)
         Bt(ix,2) =   R43(iz)*BR14X(iy,3) - R43(iy)*BR14X(iz,3)
         Bt(ix,3) = - R42(iz)*BR14X(iy,3) + R42(iy)*BR14X(iz,3)
*
         Bt(ix,4) = - (Bt(ix,1)+Bt(ix,2)+Bt(ix,3))
*
      End Do
#ifdef _DEBUG_
      Call RecPrt('Outofp: R43',' ',R43,1,3)
      Call RecPrt('Outofp: R43',' ',R42,1,3)
      Call RecPrt('Outofp: BR14X',' ',BR14X,3,3)
      Call RecPrt('Outofp: B matrix',' ',Bt,3,nCent)
#endif
*
      If (ldB) Then
*
*------- Compute the derivative of the WDC matrix.
*
         call dcopy_(12**2, [.9D1],0,dBt,1)
         Do ix = 1, 3
            iy = Mod(ix+1, 4)+(ix+1)/4
            iz = Mod(iy+1, 4)+(iy+1)/4
            Do jx = 1, ix
               jy = Mod(jx+1, 4)+(jx+1)/4
               jz = Mod(jy+1, 4)+(jy+1)/4
*                                                                      *
************************************************************************
*                                                                      *
*------------- Do Block (1,1), (2,1), (3,1). Construct
*              block (1,2) and (1,3) by symmetry.
*
               dBt(ix,1,jx,1) = - dBR14X(ix,1,jx,1)
*
               dBt(ix,2,jx,1) = + R43(iz)*dBR14X(iy,3,jx,1)
     &                          - R43(iy)*dBR14X(iz,3,jx,1)
               dBt(jx,1,ix,2) = dBt(ix,2,jx,1)
*
               dBt(ix,3,jx,1) = - R42(iz)*dBR14X(iy,3,jx,1)
     &                          + R42(iy)*dBR14X(iz,3,jx,1)
               dBt(jx,1,ix,3) = dBt(ix,3,jx,1)
*
*------------- Do block (4,1) by translational invariance and
*              (1,4) by symmetry
*
               dBt(ix,4,jx,1) = - ( dBt(ix,1,jx,1) + dBt(ix,2,jx,1)
     &                            + dBt(ix,3,jx,1) )
               dBt(jx,1,ix,4) = dBt(ix,4,jx,1)
               If (ix.ne.jx) Then
                  dBt(jx,1,ix,1) = dBt(ix,1,jx,1)
                  dBt(jx,2,ix,1) = + R43(jz)*dBR14X(jy,3,ix,1)
     &                             - R43(jy)*dBR14X(jz,3,ix,1)
                  dBt(ix,1,jx,2) = dBt(jx,2,ix,1)
                  dBt(jx,3,ix,1) = - R42(jz)*dBR14X(jy,3,ix,1)
     &                             + R42(jy)*dBR14X(jz,3,ix,1)
                  dBt(ix,1,jx,3) = dBt(jx,3,ix,1)
                  dBt(jx,4,ix,1) = - ( dBt(jx,1,ix,1) + dBt(jx,2,ix,1)
     &                               + dBt(jx,3,ix,1) )
                  dBt(ix,1,jx,4) = dBt(jx,4,ix,1)
               End if
*                                                                      *
************************************************************************
*                                                                      *
*------------- Do block (2,2), and (3,2). Construct block
*              (2,3) by symmetry
*
               dBt(ix,2,jx,2) =- R43(iz)*(R43(jz)*dBR14X(iy,3,jy,3)
     &                                   -R43(jy)*dBR14X(iy,3,jz,3))
     &                         + R43(iy)*(R43(jz)*dBR14X(iz,3,jy,3)
     &                                   -R43(jy)*dBR14X(iz,3,jz,3))
               dBt(ix,3,jx,2) =+ R42(iz)*(R43(jz)*dBR14X(iy,3,jy,3)
     &                                   -R43(jy)*dBR14X(iy,3,jz,3))
     &                         - R42(iy)*(R43(jz)*dBR14X(iz,3,jy,3)
     &                                   -R43(jy)*dBR14X(iz,3,jz,3))
               If (ix.eq.jz)  dBt(ix,3,jx,2) =  dBt(ix,3,jx,2)
     &                                       + BR14X(jy,3)
               If (ix.eq.jy)  dBt(ix,3,jx,2) =  dBt(ix,3,jx,2)
     &                                       - BR14X(jz,3)
               dBt(jx,2,ix,3) = dBt(ix,3,jx,2)
*
*              Do block (4,2) by translational invariance and (2,4) by
*              symmetry
*
               dBt(ix,4,jx,2) = - ( dBt(ix,1,jx,2) + dBt(ix,2,jx,2)
     &                            + dBt(ix,3,jx,2) )
               dBt(jx,2,ix,4) = dBt(ix,4,jx,2)
               If (ix.ne.iy) Then
                  dBt(jx,2,ix,2) = dBt(ix,2,jx,2)
*
                  dBt(jx,3,ix,2) =+ R42(jz)*(R43(iz)*dBR14X(jy,3,iy,3)
     &                                      -R43(iy)*dBR14X(jy,3,iz,3))
     &                            - R42(jy)*(R43(iz)*dBR14X(jz,3,iy,3)
     &                                      -R43(iy)*dBR14X(jz,3,iz,3))
                  If (jx.eq.iz)  dBt(jx,3,ix,2) =  dBt(jx,3,ix,2)
     &                                          + BR14X(iy,3)
                  If (jx.eq.iy)  dBt(jx,3,ix,2) =  dBt(jx,3,ix,2)
     &                                          - BR14X(iz,3)
                  dBt(ix,2,jx,3) = dBt(jx,3,ix,2)
                  dBt(jx,4,ix,2) = - ( dBt(jx,1,ix,2) + dBt(jx,2,ix,2)
     &                               + dBt(jx,3,ix,2) )
                  dBt(ix,2,jx,4) = dBt(jx,4,ix,2)
               End If
*                                                                      *
************************************************************************
*                                                                      *
*------------- Do block (3,3)
*
               dBt(ix,3,jx,3) =- R42(iz)*(R42(jz)*dBR14X(iy,3,jy,3)
     &                                   -R42(jy)*dBR14X(iy,3,jz,3))
     &                         + R42(iy)*(R42(jz)*dBR14X(iz,3,jy,3)
     &                                   -R42(jy)*dBR14X(iz,3,jz,3))
*
*------------- Do (4,3) byh translational invariance and (3,4) by
*              symmetry
*
               dBt(ix,4,jx,3) = - ( dBt(ix,1,jx,3) + dBt(ix,2,jx,3)
     &                            + dBt(ix,3,jx,3) )
               dBt(jx,3,ix,4) = dBt(ix,4,jx,3)
               If (ix.ne.iy) Then
                  dBt(jx,3,ix,3) = dBt(ix,3,jx,3)
                  dBt(jx,4,ix,3) = - ( dBt(jx,1,ix,3) + dBt(jx,2,ix,3)
     &                               + dBt(jx,3,ix,3) )
                  dBt(ix,3,jx,4) = dBt(jx,4,ix,3)
               End If
*                                                                      *
************************************************************************
*                                                                      *
*              Finally do (4,4) by translational invariance and
*              symmetry
*
               dBt(ix,4,jx,4) = - ( dBt(ix,1,jx,4) + dBt(ix,2,jx,4)
     &                            + dBt(ix,3,jx,4) )
               If (ix.ne.jx) Then
                  dBt(jx,4,ix,4) = dBt(ix,4,jx,4)
               End If
*                                                                      *
************************************************************************
*                                                                      *
            End Do
         End Do
#ifdef _DEBUG_
         Call RecPrt('dBt','(4(3F7.2,2X))',dBt,12,12)
#endif
      End If
      Call DScal_(12,   -One, Bt,1)
C     Call DScal_(12**2,-One,dBt,1)
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _Test_Numerical_
*
*     test of Bt
*
      delta=1.0D-5
      Do iAtom=1,4
         Do iCar=1,3
            tmp=xyz(iCar,iAtom)
            xyz(iCar,iAtom)= tmp + delta
            Call OutofP0(xyz,4,Tetap,ldB)
            xyz(iCar,iAtom)= tmp - delta
            Call OutofP0(xyz,4,Tetam,ldB)
            xyz(iCar,iAtom)= tmp
*
            dbdx=(Tetap-Tetam)/(2.0D0*delta)
            If (Abs(dbdx-Bt(iCar,iAtom)).gt.delta) Then
               Write (6,*) dbdx, Bt(iCar,iAtom)
               call Abend()
            End If
         End Do
      End Do
*
*     test of dBt
*
      If (ldb) Then
      Do iAtom = 1, 4
         Do iCar = 1, 3
            tmp=xyz(iCar,iAtom)
            xyz(iCar,iAtom)= tmp + delta
            Call OutofP1(xyz,4,Tetap,Bt_temp_p,ldB)
            xyz(iCar,iAtom)= tmp - delta
            Call OutofP1(xyz,4,Tetam,Bt_temp_m,ldB)
            xyz(iCar,iAtom)= tmp
*
            Do jAtom=1, 4
               Do jCar=1, 3
                  ddbddx=(Bt_temp_p(jCar,jAtom)-Bt_temp_m(jCar,jAtom))
     &                   / (2.0D0*delta)
               End Do
               If (Abs(ddbddx-dBt(iCar,iAtom,jCar,jAtom)).gt.delta)
     &            Then
                  Write (6,*) ddbddx, dBt(iCar,iAtom,jCar,jAtom)
                  call Abend()
               End If
            End Do
*
         End Do ! iCar
      End Do    ! iAtom
      End If
#endif
      Call qExit('OutofP')
      Return
      End
#ifdef _Test_Numerical_
      Subroutine OutofP0(xyz,nCent,Teta,ldB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 xyz(3,nCent), C14X(3,3),  BR14X(3,3), dBR14X(3,3,3,3),
     &       R42(3),      R43(3)
      Logical ldB
      Character*8 Label
*
      Call qEnter('OutofP0')
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     First some diagnostics
*
*     4-->1 (Bond)
      RX1=xyz(1,1)-xyz(1,4)
      RY1=xyz(2,1)-xyz(2,4)
      RZ1=xyz(3,1)-xyz(3,4)
      R41KV=RX1*RX1+RY1*RY1+RZ1*RZ1
      Q41=SQRT(R41KV)
      e41x = RX1 / Q41
      e41y = RY1 / Q41
      e41z = RZ1 / Q41
*     4-->2 (Bond in plane)
      RX2=xyz(1,2)-xyz(1,4)
      RY2=xyz(2,2)-xyz(2,4)
      RZ2=xyz(3,2)-xyz(3,4)
      R42KV=RX2*RX2+RY2*RY2+RZ2*RZ2
      Q42=SQRT(R42KV)
      e42x = RX2 / Q42
      e42y = RY2 / Q42
      e42z = RZ2 / Q42
*     4-->3 (Bond in plane)
      RX3=xyz(1,3)-xyz(1,4)
      RY3=xyz(2,3)-xyz(2,4)
      RZ3=xyz(3,3)-xyz(3,4)
      R43KV=RX3*RX3+RY3*RY3+RZ3*RZ3
      Q43=SQRT(R43KV)
      e43x = RX3 / Q43
      e43y = RY3 / Q43
      e43z = RZ3 / Q43
*
*     Get the angle between e43 and e42
*
      CosFi1 = e43x*e42x + e43y*e42y + e43z*e42z
*
      Fi1=ArCos(CosFi1)
      If (Abs(CosFi1).gt.One) Call RecPrt('xyz(1)',' ',xyz,3,4)
      dFi1 = 180.D0 * Fi1 / Pi
      If (dFi1.gt.177.5D0 .or. dFi1.lt.2.5D0) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*
*     Dirty exit! This happens when an earlier structure is ill defined.
*
      If (Abs(Fi1-Pi).lt.1.0D-13) Then
         Teta=0.0D0
         Call qExit('OutofP0')
         Return
      End If
*
*     Get the angle between e41 and e43
*
      CosFi2 = e41x*e43x + e41y*e43y + e41z*e43z
*
      Fi2=ArCos(CosFi2)
      If (Abs(CosFi2).gt.One) Call RecPrt('xyz(2)',' ',xyz,3,4)
      dFi2 = 180.D0 * Fi2 / Pi
      If (dFi2.gt.177.5D0 .or. dFi2.lt.2.5D0) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*
*     Get the angle between e41 and e42
*
      CosFi3 = e41x*e42x + e41y*e42y + e41z*e42z
*
      Fi3=ArCos(CosFi3)
      If (Abs(CosFi3).gt.One) Call RecPrt('xyz(3)',' ',xyz,3,4)
      dFi3 = 180.D0 * Fi3 / Pi
      If (dFi3.gt.177.5D0 .or. dFi3.lt.2.5D0) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call RecPrt('xyz',' ',xyz,3,nCent)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     The first two centers are trivially
*
      call dcopy_(3,xyz(1,1),1,C14X(1,1),1)
      call dcopy_(3,xyz(1,4),1,C14X(1,2),1)
*
*     The 3rd is
*
      R42(1)=xyz(1,2)-xyz(1,4)
      R42(2)=xyz(2,2)-xyz(2,4)
      R42(3)=xyz(3,2)-xyz(3,4)
      R43(1)=xyz(1,3)-xyz(1,4)
      R43(2)=xyz(2,3)-xyz(2,4)
      R43(3)=xyz(3,3)-xyz(3,4)
      C14X(1,3)=R42(2)*R43(3)-R42(3)*R43(2)
      C14X(2,3)=R42(3)*R43(1)-R42(1)*R43(3)
      C14X(3,3)=R42(1)*R43(2)-R42(2)*R43(1)
*
*     Exit if 2-3-4 are collinear
*     (equivalent to the above check, but this is more concrete)
*
      If ((C14X(1,3)**2+C14X(2,3)**2+C14X(3,3)**2).lt.1.0D-10) Then
         Teta=0.0D0
         Call qExit('OutofP0')
         Return
      End If
      C14X(1,3)=C14X(1,3)+xyz(1,4)
      C14X(2,3)=C14X(2,3)+xyz(2,4)
      C14X(3,3)=C14X(3,3)+xyz(3,4)
*
      mCent=3
      Call Bend(C14X,mCent,Teta,BR14X,.False.,.False.,Label,dBR14X,
     &          ldB)
*
      Teta = Teta - Pi/Two
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('OutofP0')
      Return
      End
      Subroutine OutofP1(xyz,nCent,Teta,Bt,ldB)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 Bt(3,nCent), xyz(3,nCent),
     &       R42(3),      R43(3),
     &       C14X(3,3),     BR14X(3,3), dBR14X(3,3,3,3)
      Logical ldB
      Character*8 Label
*
      Call qEnter('OutofP1')
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*     First some diagnostics
*
*     4-->1 (Bond)
      RX1=xyz(1,1)-xyz(1,4)
      RY1=xyz(2,1)-xyz(2,4)
      RZ1=xyz(3,1)-xyz(3,4)
      R41KV=RX1*RX1+RY1*RY1+RZ1*RZ1
      Q41=SQRT(R41KV)
      e41x = RX1 / Q41
      e41y = RY1 / Q41
      e41z = RZ1 / Q41
*     4-->2 (Bond in plane)
      RX2=xyz(1,2)-xyz(1,4)
      RY2=xyz(2,2)-xyz(2,4)
      RZ2=xyz(3,2)-xyz(3,4)
      R42KV=RX2*RX2+RY2*RY2+RZ2*RZ2
      Q42=SQRT(R42KV)
      e42x = RX2 / Q42
      e42y = RY2 / Q42
      e42z = RZ2 / Q42
*     4-->3 (Bond in plane)
      RX3=xyz(1,3)-xyz(1,4)
      RY3=xyz(2,3)-xyz(2,4)
      RZ3=xyz(3,3)-xyz(3,4)
      R43KV=RX3*RX3+RY3*RY3+RZ3*RZ3
      Q43=SQRT(R43KV)
      e43x = RX3 / Q43
      e43y = RY3 / Q43
      e43z = RZ3 / Q43
*
*     Get the angle between e43 and e42
*
      CosFi1 = e43x*e42x + e43y*e42y + e43z*e42z
*
      Fi1=ArCos(CosFi1)
      If (Abs(CosFi1).gt.One) Call RecPrt('xyz(1)',' ',xyz,3,4)
      dFi1 = 180.D0 * Fi1 / Pi
      If (dFi1.gt.177.5D0 .or. dFi1.lt.2.5D0) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*
*     Dirty exit! This happens when an earlier structure is ill defined.
*
      If (Abs(Fi1-Pi).lt.1.0D-13) Then
         Teta=0.0D0
         Call FZero(Bt,3*nCent)
         Call qExit('OutofP1')
         Return
      End If
*
*     Get the angle between e41 and e43
*
      CosFi2 = e41x*e43x + e41y*e43y + e41z*e43z
*
      Fi2=ArCos(CosFi2)
      If (Abs(CosFi2).gt.One) Call RecPrt('xyz(2)',' ',xyz,3,4)
      dFi2 = 180.D0 * Fi2 / Pi
      If (dFi2.gt.177.5D0 .or. dFi2.lt.2.5D0) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*
*     Get the angle between e41 and e42
*
      CosFi3 = e41x*e42x + e41y*e42y + e41z*e42z
*
      Fi3=ArCos(CosFi3)
      If (Abs(CosFi3).gt.One) Call RecPrt('xyz(3)',' ',xyz,3,4)
      dFi3 = 180.D0 * Fi3 / Pi
      If (dFi3.gt.177.5D0 .or. dFi3.lt.2.5D0) Then
         Write (6,*) 'Warning: auxiliary Angle close to end of range'
      End If
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUG_
      Call RecPrt('xyz',' ',xyz,3,nCent)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     The first two centers are trivially
*
      call dcopy_(3,xyz(1,1),1,C14X(1,1),1)
      call dcopy_(3,xyz(1,4),1,C14X(1,2),1)
*
*     The 3rd is
*
      R42(1)=xyz(1,2)-xyz(1,4)
      R42(2)=xyz(2,2)-xyz(2,4)
      R42(3)=xyz(3,2)-xyz(3,4)
      R43(1)=xyz(1,3)-xyz(1,4)
      R43(2)=xyz(2,3)-xyz(2,4)
      R43(3)=xyz(3,3)-xyz(3,4)
      C14X(1,3)=R42(2)*R43(3)-R42(3)*R43(2)
      C14X(2,3)=R42(3)*R43(1)-R42(1)*R43(3)
      C14X(3,3)=R42(1)*R43(2)-R42(2)*R43(1)
*
*     Exit if 2-3-4 are collinear
*     (equivalent to the above check, but this is more concrete)
*
      If ((C14X(1,3)**2+C14X(2,3)**2+C14X(3,3)**2).lt.1.0D-10) Then
         Teta=0.0D0
         Call FZero(Bt,3*nCent)
         Call qExit('OutofP1')
         Return
      End If
      C14X(1,3)=C14X(1,3)+xyz(1,4)
      C14X(2,3)=C14X(2,3)+xyz(2,4)
      C14X(3,3)=C14X(3,3)+xyz(3,4)
*
      mCent=3
      Call Bend(C14X,mCent,Teta,BR14X,.False.,.False.,Label,dBR14X,
     &          ldB)
*
      Teta = Teta - Pi/Two
      dTeta= 180.D0 * Teta/ Pi
      If (dTeta.gt.87.5D0 .or. dTeta.lt.-87.5D0) Then
         Write (6,*) 'Warning: Out of plane angle close to'
     &        //' end of range'
      End If
*
*---- Compute the WDC matrix
*
      Do ix = 1, 3
         iy = Mod(ix+1, 4)+(ix+1)/4
         iz = Mod(iy+1, 4)+(iy+1)/4
*
         Bt(ix,1) = - BR14X(ix,1)
         Bt(ix,2) =   R43(iz)*BR14X(iy,3) - R43(iy)*BR14X(iz,3)
         Bt(ix,3) = - R42(iz)*BR14X(iy,3) + R42(iy)*BR14X(iz,3)
*
         Bt(ix,4) = - (Bt(ix,1)+Bt(ix,2)+Bt(ix,3))
*
      End Do
#ifdef _DEBUG_
      Call RecPrt('Outofp: R43',' ',R43,1,3)
      Call RecPrt('Outofp: R43',' ',R42,1,3)
      Call RecPrt('Outofp: BR14X',' ',BR14X,3,3)
      Call RecPrt('Outofp: B matrix',' ',Bt,3,nCent)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('OutofP1')
      Return
      End
#endif
