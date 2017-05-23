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
      Subroutine Assemble_PPGrd(Final,nZeta,la,lb,iZeta,Alpha,Beta,
     &                          A_laplb,A_lamlb,A_lalbp,A_lalbm,JfGrad)
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       A_laplb((la+2)*(la+3)/2,(lb+1)*(lb+2)/2),
     &       A_lamlb((la+0)*(la+1)/2,(lb+1)*(lb+2)/2),
     &       A_lalbp((la+1)*(la+2)/2,(lb+2)*(lb+3)/2),
     &       A_lalbm((la+1)*(la+2)/2,(lb+0)*(lb+1)/2)
      Logical JfGrad(3,2)
*                                                                      *
************************************************************************
*                                                                      *
*     Statement function for cartesian index
*
      Ind(iy,iz) = (iy+iz)*(iy+iz+1)/2 + iz + 1
*                                                                      *
************************************************************************
*                                                                      *
C     Call RecPrt('X',' ',A_laplb,(la+2)*(la+3)/2,(lb+1)*(lb+2)/2)
C     If (la.gt.0)
C    &Call RecPrt('X',' ',A_lamlb,(la+0)*(la+1)/2,(lb+1)*(lb+2)/2)
C     Call RecPrt('X',' ',A_lalbp,(la+1)*(la+2)/2,(lb+2)*(lb+3)/2)
C     If (lb.gt.0)
C    &Call RecPrt('X',' ',A_lalbm,(la+1)*(la+2)/2,(lb+0)*(lb+1)/2)
      Do ix = la, 0, -1
         Do iy = la-ix, 0, -1
            iz = la-ix-iy
*
            Do jx = lb, 0, -1
               Do jy = lb-jx, 0, -1
                  jz = lb-jx-jy
*
                  i6 = 0
*                                                                      *
************************************************************************
*                                                                      *
*                 Ax
*
                  If (.Not.JfGrad(1,1)) Go To 102
                  i6 = i6 + 1
                  If (ix.eq.0) Then
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Alpha *A_laplb(Ind(iy  ,iz  ),Ind(jy,jz))
                  Else
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz  ),i6)=
     &                 Two*Alpha *A_laplb(Ind(iy  ,iz  ),Ind(jy,jz))
     &                -DBLE(ix)*A_lamlb(Ind(iy  ,iz  ),Ind(jy,jz))
                  End If
 102              Continue
*                                                                      *
************************************************************************
*                                                                      *
*                 Bx
*
                  If (.Not.JfGrad(1,2)) Go To 105
                  i6 = i6 + 1
                  If (jx.eq.0) Then
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Beta  *A_lalbp(Ind(iy,iz),Ind(jy  ,jz  ))
                  Else
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Beta  *A_lalbp(Ind(iy,iz),Ind(jy  ,jz  ))
     &                -DBLE(jx)*A_lalbm(Ind(iy,iz),Ind(jy  ,jz  ))
                  End If
 105              Continue
*                                                                      *
************************************************************************
*                                                                      *
*                 Ay
*
                  If (.Not.JfGrad(2,1)) Go To 103
                  i6 = i6 + 1
                  If (iy.eq.0) Then
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Alpha *A_laplb(Ind(iy+1,iz  ),Ind(jy,jz))
                  Else
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Alpha *A_laplb(Ind(iy+1,iz  ),Ind(jy,jz))
     &                -DBLE(iy)*A_lamlb(Ind(iy-1,iz  ),Ind(jy,jz))
                  End If
 103              Continue
*                                                                      *
************************************************************************
*                                                                      *
*                 By
*
                  If (.Not.JfGrad(2,2)) Go To 106
                  i6 = i6 + 1
                  If (jy.eq.0) Then
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Beta  *A_lalbp(Ind(iy,iz),Ind(jy+1,jz  ))
                  Else
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Beta  *A_lalbp(Ind(iy,iz),Ind(jy+1,jz  ))
     &                -DBLE(jy)*A_lalbm(Ind(iy,iz),Ind(jy-1,jz  ))
                  End If
 106              Continue
*                                                                      *
************************************************************************
*                                                                      *
*                 Az
*
                  If (.Not.JfGrad(3,1)) Go To 104
                  i6 = i6 + 1
                  If (iz.eq.0) Then
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Alpha *A_laplb(Ind(iy  ,iz+1),Ind(jy,jz))
                  Else
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Alpha *A_laplb(Ind(iy  ,iz+1),Ind(jy,jz))
     &                -DBLE(iz)*A_lamlb(Ind(iy  ,iz-1),Ind(jy,jz))
                  End If
 104              Continue
*                                                                      *
************************************************************************
*                                                                      *
*                 Bz
*
                  If (.Not.JfGrad(3,2)) Go To 107
                  i6 = i6 + 1
                  If (jz.eq.0) Then
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Beta  *A_lalbp(Ind(iy,iz),Ind(jy  ,jz+1))
                  Else
                     Final(iZeta,Ind(iy,iz),Ind(jy,jz),i6)=
     &                 Two*Beta  *A_lalbp(Ind(iy,iz),Ind(jy  ,jz+1))
     &                -DBLE(jz)*A_lalbm(Ind(iy,iz),Ind(jy  ,jz-1))
                  End If
 107              Continue
*                                                                      *
************************************************************************
*                                                                      *
*
               End Do
            End Do
         End Do
      End Do
C     Call RecPrt('Final',' ',Final,nZeta*nElem(la)*nElem(lb),6)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
