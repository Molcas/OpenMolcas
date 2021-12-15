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
* Copyright (C) Yannick Carissan                                       *
*               Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine UpdateB(Col,nOrb2Loc,ipLbl,nComp,Gamma_rot,iMO_s,iMO_t,
     &                   Debug)
C
C     Author: T.B. Pedersen
C
C     Purpose: update MO dipole matrices for Boys localisation.
C              (Almost exact copy of UpdateP by Y. Carissan.)
C
      Implicit Real*8 (a-h,o-z)
      Real*8  Col(nOrb2Loc,2)
      Integer ipLbl(nComp)
      Logical Debug
#include "WrkSpc.fh"

      Character*18 Label

      cosg   = cos(Gamma_rot)
      sing   = sin(Gamma_rot)
      cos2g  = cosg*cosg
      sin2g  = sing*sing
      cosing = cosg*sing

      Do iComp = 1,nComp

         ip  = ipLbl(iComp)
         ip0 = ip - 1

         kOff_s = ip0 + nOrb2Loc*(iMO_s-1)
         kOff_t = ip0 + nOrb2Loc*(iMO_t-1)
         kOff_ss = kOff_s + iMO_s
         kOff_ts = kOff_s + iMO_t
         kOff_st = kOff_t + iMO_s
         kOff_tt = kOff_t + iMO_t
         Dss = Work(kOff_ss)
         Dst = Work(kOff_st)
         Dtt = Work(kOff_tt)
#if defined (_DEBUGPRINT_)
         Dts = Work(kOff_ts)
         Tst = Dst - Dts
         If (abs(Tst) .gt. 1.0d-14) Then
            Write(6,*) 'Broken symmetry in UpdateB!!'
            Write(6,*) 'MOs s and t: ',iMO_s,iMO_t
            Write(6,*) 'Component  : ',iComp
            Write(6,*) 'Dst  = ',Dst
            Write(6,*) 'Dts  = ',Dts
            Write(6,*) 'Diff = ',Tst
            Call SysAbendMsg('UpdateB','Broken symmetry!','[1]')
         End If
#endif

         Call dCopy_(nOrb2Loc,Work(kOff_s+1),1,Col(1,1),1)
         Call dCopy_(nOrb2Loc,Work(kOff_t+1),1,Col(1,2),1)

         Call dScal_(nOrb2Loc,cosg,Work(kOff_s+1),1)
         Call dAXPY_(nOrb2Loc,sing,Col(1,2),1,Work(kOff_s+1),1)
         Call dScal_(nOrb2Loc,cosg,Work(kOff_t+1),1)
         Call dAXPY_(nOrb2Loc,-sing,Col(1,1),1,Work(kOff_t+1),1)

         Work(kOff_s+iMO_s) = Dss*cos2g + Dtt*sin2g
     &                      + 2.0d0*Dst*cosing
         Work(kOff_s+iMO_t) = (Dtt-Dss)*cosing + Dst*(cos2g-sin2g)
         Work(kOff_t+iMO_s) = Work(kOff_s+iMO_t)
         Work(kOff_t+iMO_t) = Dtt*cos2g + Dss*sin2g
     &                      - 2.0d0*Dst*cosing

         Call dCopy_(nOrb2Loc,Work(kOff_s+1),1,Work(ip0+iMO_s),nOrb2Loc)
         Call dCopy_(nOrb2Loc,Work(kOff_t+1),1,Work(ip0+iMO_t),nOrb2Loc)

#if defined (_DEBUGPRINT_)
         Dst = Work(kOff_st)
         Dts = Work(kOff_ts)
         Tst = Dst - Dts
         If (abs(Tst) .gt. 1.0d-14) Then
            Write(6,*) 'Broken symmetry in UpdateB!!'
            Write(6,*) 'MOs s and t: ',iMO_s,iMO_t
            Write(6,*) 'Component  : ',iComp
            Write(6,*) 'Dst  = ',Dst
            Write(6,*) 'Dts  = ',Dts
            Write(6,*) 'Diff = ',Tst
            Call SysAbendMsg('UpdateB','Broken symmetry!','[2]')
         End If
#endif

      End Do

      If (Debug) Then
         Write(6,*) 'In UpdateB'
         Write(6,*) '----------'
         Do iComp = 1,nComp
            Write(Label,'(A,I2,A,I4)') 'MO Dip',iComp,'   col',iMO_s
            ip = ipLbl(iComp) + nOrb2Loc*(iMO_s-1)
            Call RecPrt(Label,' ',Work(ip),nOrb2Loc,1)
            Write(Label,'(A,I2,A,I4)') 'MO Dip',iComp,'   col',iMO_t
            ip = ipLbl(iComp) + nOrb2Loc*(iMO_t-1)
            Call RecPrt(Label,' ',Work(ip),nOrb2Loc,1)
         End Do
      End If

      End
