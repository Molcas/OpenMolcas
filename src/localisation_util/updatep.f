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
*               2005, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine UpdateP(PACol,Name,nBas_Start,
     &                   nOrb2Loc,nAtoms,iTab_ptr,PA,gamma_rot,
     &                   iMO_s,iMO_t,Debug)
c
c     Author: Yannick Carissan.
c
c     Modifications:
c        - October 6, 2005 (Thomas Bondo Pedersen):
c          Reduce operation count and use BLAS.
c
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "WrkSpc.fh"
#include "Molcas.fh"
      Integer iTab_Ptr(*),nBas_Start(*)
      Real*8 PACol(nOrb2Loc,2)
      Real*8 PA(nOrb2Loc,nOrb2Loc,nAtoms)
      Character*(LENIN8) Name(*),PALbl
      Logical Debug
c
      cosg   = cos(gamma_rot)
      sing   = sin(gamma_rot)
      cos2g  = cosg*cosg
      sin2g  = sing*sing
      cosing = cosg*sing
c
      Do iAt=1,nAtoms
#define _TEST4_
#ifdef _TEST4_
*       Call RecPrt('PA(1,1,iAt)',' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
c
c------ Copy out the PAss, PAtt, and PAst elements.
c
        PAss = PA(iMO_s,iMO_s,iAt)
        PAst = PA(iMO_s,iMO_t,iAt)
        PAtt = PA(iMO_t,iMO_t,iAt)
*       Write (6,*) 'updateP:',PAss,PAst,PAtt
#else
c
c------ The array iTab_ptr contains the value of the pointer on the
c       PA array for atom iAt
c
        ip  = iTab_ptr(iAt)
        ip0 = ip - 1
*       Call RecPrt('Work(ip)',' ',Work(ip),nOrb2Loc,nOrb2Loc)
c
c------ Copy out the PAss, PAtt, and PAst elements.
c
        kOff_s = ip0 + nOrb2Loc*(iMO_s-1)
        kOff_t = ip0 + nOrb2Loc*(iMO_t-1)
        kOff_ss = kOff_s + iMO_s
        kOff_ts = kOff_s + iMO_t
        kOff_st = kOff_t + iMO_s
        kOff_tt = kOff_t + iMO_t
        PAss = Work(kOff_ss)
        PAst = Work(kOff_st)
        PAtt = Work(kOff_tt)
*       Write (6,*) 'updateP:',PAss,PAst,PAtt
#endif
#if defined (_DEBUG_)
#ifdef _TEST4_
        PAts = PA(iMO_t,iMO_s,iAt)
#else
        PAts = Work(kOff_ts)
#endif
        Tst  = PAst - PAts
        If (abs(Tst) .gt. 1.0d-14) Then
           Write(6,*) 'Broken symmetry in UpdateP!!'
           Write(6,*) 'MOs s and t: ',iMO_s,iMO_t
           Write(6,*) 'PAst = ',PAst
           Write(6,*) 'PAts = ',PAts
           Write(6,*) 'Diff = ',Tst
           Call SysAbendMsg('UpdateP','Broken symmetry!',' ')
        End If
#endif
c
c------ Copy out columns s and t of PA.
c
#ifdef _TEST4_
        Call dCopy_(nOrb2Loc,PA(1,iMO_s,iAt),1,PACol(1,1),1)
        Call dCopy_(nOrb2Loc,PA(1,iMO_t,iAt),1,PACol(1,2),1)
#else
        Call dCopy_(nOrb2Loc,Work(kOff_s+1),1,PACol(1,1),1)
        Call dCopy_(nOrb2Loc,Work(kOff_t+1),1,PACol(1,2),1)
#endif
c
c------ Compute transformed columns.
c
#ifdef _TEST4_
        Call dScal_(nOrb2Loc,cosg,PA(1,iMO_s,iAt),1)
        Call dAXPY_(nOrb2Loc, sing,PACol(1,2),1,PA(1,iMO_s,iAt),1)
        Call dScal_(nOrb2Loc,cosg,PA(1,iMO_t,iAt),1)
        Call dAXPY_(nOrb2Loc,-sing,PACol(1,1),1,PA(1,iMO_t,iAt),1)
#else
        Call dScal_(nOrb2Loc,cosg,Work(kOff_s+1),1)
        Call dAXPY_(nOrb2Loc,sing,PACol(1,2),1,Work(kOff_s+1),1)
        Call dScal_(nOrb2Loc,cosg,Work(kOff_t+1),1)
        Call dAXPY_(nOrb2Loc,-sing,PACol(1,1),1,Work(kOff_t+1),1)
#endif
c
c------ Compute PAss, PAtt, PAst, and PAts (= PAst).
c
#ifdef _TEST4_
        PA(iMO_s,iMO_s,iAt)= PAss*cos2g + PAtt*sin2g
     &                     + Two*PAst*cosing
        PA(iMO_t,iMO_s,iAt)= (PAtt-PAss)*cosing + PAst*(cos2g-sin2g)
        PA(iMO_s,iMO_t,iAt)= PA(iMO_t,iMO_s,iAt)
        PA(iMO_t,iMO_t,iAt)= PAtt*cos2g + PAss*sin2g
     &                     - Two*PAst*cosing
#else
        Work(kOff_s+iMO_s) = PAss*cos2g + PAtt*sin2g
     &                     + Two*PAst*cosing
        Work(kOff_s+iMO_t) = (PAtt-PAss)*cosing + PAst*(cos2g-sin2g)
        Work(kOff_t+iMO_s) = Work(kOff_s+iMO_t)
        Work(kOff_t+iMO_t) = PAtt*cos2g + PAss*sin2g
     &                     - Two*PAst*cosing
#endif
c
c------ Copy columns to rows.
c
#ifdef _TEST4_
        Call dCopy_(nOrb2Loc,PA(1,iMO_s,iAt),1,PA(iMO_s,1,iAt),nOrb2Loc)
        Call dCopy_(nOrb2Loc,PA(1,iMO_t,iAt),1,PA(iMO_t,1,iAt),nOrb2Loc)
#else
        Call dCopy_(nOrb2Loc,Work(kOff_s+1),1,Work(ip0+iMO_s),nOrb2Loc)
        Call dCopy_(nOrb2Loc,Work(kOff_t+1),1,Work(ip0+iMO_t),nOrb2Loc)
#endif
c
      End Do
c
      If (Debug) Then
        Write(6,*) 'In UpdateP'
        Write(6,*) '----------'
        Do iAt=1,nAtoms
          PALbl='PA__'//Name(nBas_Start(iAt))(1:LENIN)
          ip=iTab_ptr(iAt)
#ifdef _TEST4_
          Call RecPrt(PALbl,' ',PA(1,1,iAt),nOrb2Loc,nOrb2Loc)
#else
          Call RecPrt(PALbl,' ',Work(ip),nOrb2Loc,nOrb2Loc)
#endif
        End Do
      End If
c
      Return
      End
