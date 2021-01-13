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
      Subroutine Make_Labels(LblCbs,LblSbs,MxFnc,lMax)
      Implicit Real*8 (a-h,o-z)
#include "angtp.fh"
      Character*8 LblCBs(MxFnc), LblSBs(MxFnc)
      Character Sgn*3
!#define _DEBUGPRINT_
*
*---- Generate cartesian labels
*
      lxyz=0
      Do ixyz = 0, lMax
         Do ix = ixyz, 0, -1
            iyMax=ixyz-ix
            Do iy = iyMax, 0 , -1
               lxyz=lxyz+1
               iz=ixyz-ix-iy
*              Form labels for cartesian basis functions
               Write (LblCBs(lxyz),'(A,3I2.2)') AngTp(ixyz),ix,iy,iz
            End Do
         End Do
      End Do
      If (lMax.ge.0) LblCBs(1) = '01s     '
      If (lMax.ge.1) Then
         LblCBs(2) = '02px    '
         LblCBs(3) = '02py    '
         LblCBs(4) = '02pz    '
      End If
*
*     Do the same for the spherical gaussians.
*
      i = 0
      Do n = 0, lMax
         Do l = n, 0, -2
            Do m = -l, l
               If (m.lt.0) Then
                  Sgn='-  '
               Else If (m.gt.0) Then
                  Sgn='+  '
               Else
                  Sgn='   '
               End If
               i=i+1
               Write (LblSbs(i),'(I2.2,A,I2.2,A)')
     &            n+1,AngTp(l),Abs(m),Sgn
            End Do
         End Do
      End Do
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'lMax,MxFnc=',lMax,MxFnc
      Write (6,*)
      Write (6,*) 'LblCBs:'
      Do ixyz = 1, lxyz
         Write (6,'(A)')LblCBs(ixyz)
      End Do
      Write (6,*) 'LblSBs:'
      Do ixyz = 1, i
         Write (6,'(A)')LblCBs(ixyz)
      End Do
      Write (6,*)
#endif
*
      Return
      End
