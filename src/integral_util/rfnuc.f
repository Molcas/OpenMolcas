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
* Copyright (C) 1990, Roland Lindh                                     *
*               1990, IBM                                              *
************************************************************************
      SubRoutine RFNuc(CoOP,rNucMm,ir)
************************************************************************
*                                                                      *
* Object: to compute the multipole moments for the nuclei.             *
*                                                                      *
* Called from: Input                                                   *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             November '90                                             *
************************************************************************
      use Basis_Info
      use Center_Info
#ifdef _OBSOLETE_
      use External_Centers, only: nOrd_XF, XF
#endif
      use Phase_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "print.fh"
#include "WrkSpc.fh"
      Real*8  rNucMm((ir+1)*(ir+2)/2), CoOp(3), A(3), RA(3)
#ifdef _OBSOLETE_
     &        ,rRMy(3)
      Integer iStb(0:7), jCoSet(0:7,0:7)
#endif
*                                                                      *
************************************************************************
*                                                                      *
      iRout = 124
      iPrint = nPrint(iRout)
      Call qEnter('RFNuc')
      If (iPrint.ge.99) Then
         Call RecPrt(' In RFNuc:CoOp',' ',CoOp,1,3)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     Compute the nuclear contribution to the multipole moments
*
*     Contributions due to the charges of nuclear charges
*
      iq = 0
      Do ix = ir, 0, -1
         Do iy = ir-ix, 0, -1
            iq = iq + 1
            iz = ir-ix-iy
            temp = Zero
C           Write (*,*) ' ix,iy,iz=',ix,iy,iz
*
            ndc = 0
            Do iCnttp = 1, nCnttp
               ZA = dbsc(iCnttp)%Charge
               If (ZA.eq.Zero) Go To 101
               If (iPrint.ge.99) Then
                  Write (6,*) ' Charge=',ZA
                  Call RecPrt(' Centers',' ',dbsc(iCnttp)%Coor,3,
     &                        dbsc(iCnttp)%nCntr)
               End If
               Do iCnt = 1, dbsc(iCnttp)%nCntr
                  A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)
                  mdc = ndc + iCnt
                  Do i = 0, nIrrep/dc(mdc)%nStab - 1
                     Call OA(dc(mdc)%iCoSet(i,0),A,RA)
C                    Call RecPrt(' RA',' ',RA,1,3)
C                    Call RecPrt(' CoOp',' ',CoOp,1,3)
#ifdef NAGFOR
                     If (iCnt.lt.-2) Write (6,*) 'Nag problem'
#endif

                     If (ix.eq.0) Then
                        CCoMx=One
                     Else
                        CCoMx=(RA(1)-CoOp(1))**ix
                     End If
                     If (iy.eq.0) Then
                        CCoMy=One
                     Else
                        CCoMy=(RA(2)-CoOp(2))**iy
                     End If
                     If (iz.eq.0) Then
                        CCoMz=One
                     Else
                        CCoMz=(RA(3)-CoOp(3))**iz
                     End If
C                    Write (*,*) CCoMx, CCoMy, CCoMz, temp
                     temp = temp + ZA * CCoMx * CCoMy * CCoMz
                  End Do
               End Do
 101           Continue
               ndc = ndc + dbsc(iCnttp)%nCntr
            End Do
            rNucMm(iq) = temp
         End Do
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _OBOLETE_
*     The remainder of this subroutine is obsolete and
*     is kept only for testing reasons. It is replaced by
*     the subroutine XFMoment which is more general.
      GoTo 99





      If ((.Not.Allocated(XF)).or.(nOrd_XF.lt.0)) Go To 99
*
*     Contributions due to the charges and dipoles of the
*     static external electric field.
*
*     Write (*,*) ' Adding contibutions from esef!'

      iq = 0
      Do ix = ir, 0, -1
         Do iy = ir-ix, 0, -1
            iq = iq + 1
            iz = ir-ix-iy
            temp = Zero
*           Write (*,*) ' ix,iy,iz=',ix,iy,iz
*
            Do iFd = 1, nXF
               DAx=Zero
               DAy=Zero
               DAz=Zero
               Qxx=Zero
               Qxy=Zero
               Qxz=Zero
               Qyy=Zero
               Qyz=Zero
               Qzz=Zero
               If (nOrd_XF.eq.0) Then
                  ZA = XF(4,iFd)
               Else If (nOrd_XF.eq.1) Then
                  ZA = XF(4,iFd)
                  DAx= XF(5,iFd)
                  DAy= XF(6,iFd)
                  DAz= XF(7,iFd)
               Else If (nOrd_XF.eq.2) Then
                  ZA = XF(4,iFd)
                  DAx= XF(5,iFd)
                  DAy= XF(6,iFd)
                  DAz= XF(7,iFd)
                  Qxx= XF(8,iFd)
                  Qxy= XF(9,iFd)
                  Qxz= XF(10,iFd)
                  Qyy= XF(11,iFd)
                  Qyz= XF(12,iFd)
                  Qzz= XF(13,iFd)
               Else
                  Call WarningMessage(2,
     &                      'RFNuc: Option not implemented yet!')
                  Call Abend()
               End If
               If (iPrint.ge.99) Then
                  Write (6,*) ' Charge=',ZA
                  Write (6,*) ' ixyz=',ixyz
                  Call RecPrt(' Centers',' ',XF(1,iXF),3,1)
               End If
*
               A(1:3) = XF(1:3,iXF)
*
*------------- Generate Stabilazor of C
*
               iChxyz=iChAtm(A)
               iDum=0
               Call Stblz(iChxyz,nStb,iStb,iDum,jCoSet)
*
*              Write (*,*) ' nStb=',nStb
               Do i = 0, nIrrep/nStb - 1
                  Call OA(jCoSet(i,0),A,RA)
                  rRMy(1)=DAx*DBLE(iPhase(1,jCoSet(i,0)))
                  rRMy(2)=DAy*DBLE(iPhase(2,jCoSet(i,0)))
                  rRMy(3)=DAz*DBLE(iPhase(3,jCoSet(i,0)))
                  QRAxx = QAxx
                  QRAyy = QAyy
                  QRAzz = QAzz
                  QRAxy=DBLE(iPhase(1,jCoSet(i,0))
     &                      *iPhase(2,jCoSet(i,0)))*QAxy
                  QRAxz=DBLE(iPhase(1,jCoSet(i,0))
     &                      *iPhase(3,jCoSet(i,0)))*QAxz
                  QRAyz=DBLE(iPhase(2,jCoSet(i,0))
     &                      *iPhase(3,jCoSet(i,0)))*QAyz

                  If (ix.eq.0) Then
                     CCoMx=One
                  Else
                     CCoMx=(RA(1)-CoOp(1))**ix
                  End If
                  If (iy.eq.0) Then
                     CCoMy=One
                  Else
                     CCoMy=(RA(2)-CoOp(2))**iy
                  End If
                  If (iz.eq.0) Then
                     CCoMz=One
                  Else
                     CCoMz=(RA(3)-CoOp(3))**iz
                  End If

*                 Write (*,*) CCoMx, CCoMy, CCoMz, temp
*
*---------------- The charge contibution
*
                  temp = temp + ZA * CCoMx * CCoMy * CCoMz
*
*---------------- Dipole contributions
*
                  If (ix.ge.1) Then
                     temp = temp + DBLE(ix)*rRmy(1)*CCoMy*CCoMz*
     &                      (RA(1)-CoOp(1))**(ix-1)
                  End If
                  If (iy.ge.1) Then
                     temp = temp + DBLE(iy)*rRmy(2)*CCoMx*CCoMz*
     &                      (RA(2)-CoOp(2))**(iy-1)
                  End If
                  If (iz.ge.1) Then
                     temp = temp + DBLE(iz)*rRmy(3)*CCoMx*CCoMy*
     &                      (RA(3)-CoOp(3))**(iz-1)
                  End If

               End Do
            End Do
c            Write (*,*) ' Temp=',temp
            rNucMm(iq) = rNucMm(iq) + temp

         End Do
      End Do
*
 99   Continue
#endif
      If (iPrint.ge.99) Call RecPrt(' Nuclear Multipole Moments',
     &                              ' ',rNucMm,ip,1)
*     Call GetMem(' Exit RFNuc','CHECK','REAL',iDum,iDum)
*                                                                      *
************************************************************************
*                                                                      *
      Call qExit('RFNuc')
      Return
      End
