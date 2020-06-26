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
* Copyright (C) 2004, Par Soderhjelm                                   *
************************************************************************
      SubRoutine EFXF(coord,XF,nXF,nOrd_XF,iXPolType,dEF,
     &     XMolnr,nXMolnr,iGrid,scal14)

************************************************************************
*                                                                      *
*     Object:  Calculate electric field in one point                   *
*              from XFIELD multipoles                                  *
*              Note: Ignores symmetry totally!                         *
*                                                                      *
*     Authors: P. Soderhjelm                                           *
*              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
*                                                                      *
*              November 2004                                           *
************************************************************************

      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 coord(3),XF(*),dEF(3)
*     Integer XMolnr(nXMolnr,nXF)
      Real*8 XMolnr(nXMolnr,nXF)

      Logical LExcl

*
*     Statement function for Cartesian index
*
      nElem(ixyz) = (ixyz+1)*(ixyz+2)/2


      If(nOrd_XF.lt.0) Return
*
*     Calculate number of entries per XFIELD point
      Inc = 3
      Do iOrdOp = 0, nOrd_XF
         Inc = Inc + nElem(iOrdOp)
      End Do
      If(iXPolType.gt.0) Inc = Inc + 6

*     Loop over XF points
      Do iFd = 1, nXF
         scal=One
         If((iXPolType.gt.0).and.(iGrid.le.nXF)) Then
            LExcl=.False.
            If(iFd.eq.iGrid) LExcl=.True.
            Do i=1,nXMolnr
               If (INT(XMolnr(1,iGrid)).eq.INT(XMolnr(i,iFd)))
     &            LExcl=.True.
               If (INT(XMolnr(1,iGrid)).eq.-INT(XMolnr(i,iFd)))
     &            scal=scal14
            EndDo
            If(LExcl) Then
c               Write(6,*)'EXCLUDE ',iFd,' from field at ',iGrid
               Goto 1
            ElseIf(scal.lt.One) Then
c               Write(6,*)'SCALE ',iFd,' from field at ',iGrid,
c     &              ' with', scal
            EndIf
         EndIf
         ZA=Zero
         DAx=Zero
         DAy=Zero
         DAz=Zero
         QAxx=Zero
         QAxy=Zero
         QAxz=Zero
         QAyy=Zero
         QAyz=Zero
         QAzz=Zero
         If (nOrd_XF.eq.0) Then
            ZA = XF((iFd-1)*Inc+4)*scal
         Else If (nOrd_XF.eq.1) Then
            ZA = XF((iFd-1)*Inc+4)*scal
            DAx= XF((iFd-1)*Inc+5)*scal
            DAy= XF((iFd-1)*Inc+6)*scal
            DAz= XF((iFd-1)*Inc+7)*scal
         Else If (nOrd_XF.eq.2) Then
            ZA = XF((iFd-1)*Inc+4)*scal
            DAx= XF((iFd-1)*Inc+5)*scal
            DAy= XF((iFd-1)*Inc+6)*scal
            DAz= XF((iFd-1)*Inc+7)*scal
            QAxx= XF((iFd-1)*Inc+8)*scal
            QAxy= XF((iFd-1)*Inc+9)*scal
            QAxz= XF((iFd-1)*Inc+10)*scal
            QAyy= XF((iFd-1)*Inc+11)*scal
            QAyz= XF((iFd-1)*Inc+12)*scal
            QAzz= XF((iFd-1)*Inc+13)*scal
         Else
            Call WarningMessage(2,'Efxf: Option not implemented yet!')
            Call Abend()
         EndIf
         x = XF((iFd-1)*Inc+1)-coord(1)
         y = XF((iFd-1)*Inc+2)-coord(2)
         z = XF((iFd-1)*Inc+3)-coord(3)
         r12 = Sqrt(x**2 + y**2 + z**2 )


*     Z field
         dEF(1)=dEF(1)-ZA*x/r12**3
         dEF(2)=dEF(2)-ZA*y/r12**3
         dEF(3)=dEF(3)-ZA*z/r12**3

         If(nOrd_XF.lt.1) Goto 1

*     D field
         dEF(1)=dEF(1)+Three*(DAx* x+DAy* y+DAz *z)*x/r12**5-DAx/r12**3
         dEF(2)=dEF(2)+Three*(DAx* x+DAy* y+DAz *z)*y/r12**5-DAy/r12**3
         dEF(3)=dEF(3)+Three*(DAx* x+DAy* y+DAz *z)*z/r12**5-DAz/r12**3

         If(nOrd_XF.lt.2) Goto 1

*     Q field
         QAsum=(QAxx*x*x+QAyy*y*y+QAzz*z*z+2.0D0*
     &        (QAxy*x*y+QAxz*x*z+QAyz*y*z))

         dEF(1)=dEF(1)+0.5D0*(-15.0D0/r12**7*x*QAsum
     &        +3.0D0/r12**5*
     &   (3.0D0*QAxx*x
     &   +2.0D0*QAxy*y
     &   +2.0D0*QAxz*z
     &   +1.0D0*QAyy*x
     &   +1.0D0*QAzz*x))

         dEF(2)=dEF(2)+0.5D0*(-15.0D0/r12**7*y*QAsum
     &        +3.0D0/r12**5*
     &   (1.0D0*QAxx*y
     &   +2.0D0*QAxy*x
     &   +3.0D0*QAyy*y
     &   +2.0D0*QAyz*z
     &   +1.0D0*QAzz*y))

         dEF(3)=dEF(3)+0.5D0*(-15.0D0/r12**7*z*QAsum
     &        +3.0D0/r12**5*
     &   (1.0D0*QAxx*z
     &   +2.0D0*QAxz*x
     &   +1.0D0*QAyy*z
     &   +2.0D0*QAyz*y
     &   +3.0D0*QAzz*z))

c     These formulas gives the corresponding energy terms in drvn0
c         eZD=ZA*(DRBx*x+DRBy*y+DRBz*z)/r12**3
c         eDD=(DAx*DRBx+DAy*DRBy+DAz*DRBz)/r12**3
c     &        -Three*(DAx* x+DAy* y+DAz *z)
c     &        *(DRBx*x+DRBy*y+DRBz*z)/r12**5c
c         eQD=-0.5D0*(-15.0D0/r12**7*
c     & (DRBx*x+DRBy*y+DRBz*z)*QAsum
c     &                   +3.0D0/r12**5*
c     &   (3.0D0*DRBx*QAxx*x
c     &   +1.0D0*DRBy*QAxx*y
c     &   +1.0D0*DRBz*QAxx*z
c     &   +2.0D0*DRBx*QAxy*y
c     &   +2.0D0*DRBy*QAxy*x
c     &   +2.0D0*DRBx*QAxz*z
c     &   +2.0D0*DRBz*QAxz*x
c     &   +1.0D0*DRBx*QAyy*x
c     &   +3.0D0*DRBy*QAyy*y
c     &   +1.0D0*DRBz*QAyy*z
c     &   +2.0D0*DRBy*QAyz*z
c     &   +2.0D0*DRBz*QAyz*y
c     &   +1.0D0*DRBx*QAzz*x
c     &   +1.0D0*DRBy*QAzz*y
c     &   +3.0D0*DRBz*QAzz*z))


 1    Continue
      EndDo   !iFd

      Return
      End
