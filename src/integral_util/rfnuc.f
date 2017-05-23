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
               If (Charge(iCnttp).eq.Zero) Go To 101
               ZA = Charge(iCnttp)
               ixyz = ipCntr(iCnttp)
               If (iPrint.ge.99) Then
                  Write (6,*) ' Charge=',ZA
                  Write (6,*) ' ixyz=',ixyz
                  Call RecPrt(' Centers',' ',Work(ixyz),3,nCntr(iCnttp))
               End If
               Do iCnt = 1, nCntr(iCnttp)
                  A(1) = Work(ixyz  )
                  A(2) = Work(ixyz+1)
                  A(3) = Work(ixyz+2)
                  mdc = ndc + iCnt
                  Do i = 0, nIrrep/nStab(mdc) - 1
                     RA(1)=A(1)*DBLE(iPhase(1,iCoset(i,0,mdc)))
                     RA(2)=A(2)*DBLE(iPhase(2,iCoset(i,0,mdc)))
                     RA(3)=A(3)*DBLE(iPhase(3,iCoset(i,0,mdc)))
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
                  ixyz = ixyz + 3
               End Do
 101           Continue
               ndc = ndc + nCntr(iCnttp)
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





      If ((.Not.lXF).or.(nOrd_XF.lt.0)) Go To 99
*
*     Contributions due to the charges and dipoles of the
*     static external electric field.
*
*     Write (*,*) ' Adding contibutions from esef!'

      nInp=(nOrd_XF+1)*(nOrd_XF+2)*(nOrd_XF+3)/6
      Inc=nInp+3
      If(iXPolType.gt.0) Inc = Inc + 6

      iq = 0
      Do ix = ir, 0, -1
         Do iy = ir-ix, 0, -1
            iq = iq + 1
            iz = ir-ix-iy
            temp = Zero
*           Write (*,*) ' ix,iy,iz=',ix,iy,iz
*
            ip = ipXF-1
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
                  ZA = Work(ip+(iFd-1)*Inc+4)
               Else If (nOrd_XF.eq.1) Then
                  ZA = Work(ip+(iFd-1)*Inc+4)
                  DAx= Work(ip+(iFd-1)*Inc+5)
                  DAy= Work(ip+(iFd-1)*Inc+6)
                  DAz= Work(ip+(iFd-1)*Inc+7)
               Else If (nOrd_XF.eq.2) Then
                  ZA = Work(ip+(iFd-1)*Inc+4)
                  DAx= Work(ip+(iFd-1)*Inc+5)
                  DAy= Work(ip+(iFd-1)*Inc+6)
                  DAz= Work(ip+(iFd-1)*Inc+7)
                  Qxx= Work(ip+(iFd-1)*Inc+8)
                  Qxy= Work(ip+(iFd-1)*Inc+9)
                  Qxz= Work(ip+(iFd-1)*Inc+10)
                  Qyy= Work(ip+(iFd-1)*Inc+11)
                  Qyz= Work(ip+(iFd-1)*Inc+12)
                  Qzz= Work(ip+(iFd-1)*Inc+13)
               Else
                  Call WarningMessage(2,
     &                      'RFNuc: Option not implemented yet!')
                  Call Abend()
               End If
               ixyz = ip+(iFd-1)*Inc+1
               If (iPrint.ge.99) Then
                  Write (6,*) ' Charge=',ZA
                  Write (6,*) ' ixyz=',ixyz
                  Call RecPrt(' Centers',' ',Work(ixyz),3,1)
               End If
*
               A(1) = Work(ixyz  )
               A(2) = Work(ixyz+1)
               A(3) = Work(ixyz+2)
*
*------------- Generate Stabilazor of C
*
               If (nIrrep.eq.8) Then
                  nOper=3
               Else If (nIrrep.eq.4) Then
                  nOper=2
               Else If (nIrrep.eq.2) Then
                  nOper=1
               Else
                  nOper=0
               End If
               iChxyz=iChAtm(A,iOper,nOper,iChBas(2))
               iDum=0
               Call Stblz(iChxyz,iOper,nIrrep,nStb,iStb,iDum,jCoSet)
*
*              Write (*,*) ' nStb=',nStb
               Do i = 0, nIrrep/nStb - 1
                  RA(1)=A(1)*DBLE(iPhase(1,jCoSet(i,0)))
                  RA(2)=A(2)*DBLE(iPhase(2,jCoSet(i,0)))
                  RA(3)=A(3)*DBLE(iPhase(3,jCoSet(i,0)))
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
