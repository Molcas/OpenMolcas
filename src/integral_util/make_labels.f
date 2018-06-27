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
      Subroutine Make_Labels(LblCbs,LblSbs,MxFnc,iAngMx)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "angtp.fh"
      Character*8 LblCBs(MxFnc), LblSBs(MxFnc)
      Character Sgn*4
*
*-----Shlnm
*
      If (iAngMx.ge. 0) AngTp(0)='s'
      If (iAngMx.ge. 1) AngTp(1)='p'
      If (iAngMx.ge. 2) AngTp(2)='d'
      If (iAngMx.ge. 3) AngTp(3)='f'
      If (iAngMx.ge. 4) AngTp(4)='g'
      If (iAngMx.ge. 5) AngTp(5)='h'
      If (iAngMx.ge. 6) AngTp(6)='i'
      If (iAngMx.ge. 7) AngTp(7)='k'
      If (iAngMx.ge. 8) AngTp(8)='l'
      If (iAngMx.ge. 9) AngTp(9)='m'
      If (iAngMx.ge.10) AngTp(10)='n'
      If (iAngMx.ge.11) AngTp(11)='o'
      If (iAngMx.ge.12) AngTp(12)='q'
      If (iAngMx.ge.13) AngTp(13)='r'
      If (iAngMx.ge.14) AngTp(14)='t'
      If (iAngMx.ge.15) AngTp(15)='u'
*
*---- Generate cartesian labels
*
      lxyz=0
      Do ixyz = 0, Max(iAngMx,1)
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
      If (iAngMx.ge.0) LblCBs(1) = '1s      '
      If (iAngMx.ge.1) Then
         LblCBs(2) = '2px     '
         LblCBs(3) = '2py     '
         LblCBs(4) = '2pz     '
      End If
*
*     Do the same for the spherical gaussians.
*
      i = 1
      Do n = 0, iAngMx
         Do l = n, 0, -2
            Do m = -l, l
               If (m.lt.0) Then
                  Sgn='-   '
               Else If (m.gt.0) Then
                  Sgn='+   '
               Else
                  Sgn='    '
               End If
               Write (LblSbs(i),'(I1,A,I2.2,A)') n+1,AngTp(l),Abs(m),Sgn
               i=i+1
            End Do
         End Do
      End Do
*
      Return
      End
