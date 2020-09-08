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
      Subroutine NucInd(coor,kdc,ifgrd,ifhss,indgrd,indhss,
     &                  jfgrd,jfhss,jndgrd,jndhss,tr,ifg)
      use Real_Spherical
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"
#include "disp2.fh"

      Real*8 Coor(3,4)
      Integer IndGrd(0:2,0:1,0:(nIrrep-1)),
     &        IndHss(0:1,0:2,0:1,0:2,0:(nIrrep-1))
      Logical IfHss(0:1,0:2,0:1,0:2),IfGrd(0:2,0:1), TstFnc, TF,
     &        IfG(0:3),Tr(0:3)
      Integer JndGrd(0:2,0:3,0:(nIrrep-1)),
     &        JndHss(0:3,0:2,0:3,0:2,0:(nIrrep-1))
*
      Logical JfHss(0:3,0:2,0:3,0:2),JfGrd(0:2,0:3),EQ
*                                                                      *
************************************************************************
*                                                                      *
*     Statement functions
*
      IX(i1,i2)=i1*(i1-1)/2+i2
      TF(mdc,iIrrep,iComp) = TstFnc(dc(mdc)%iCoSet,
     &                              iIrrep,iComp,dc(mdc)%nStab)
*                                                                      *
************************************************************************
*                                                                      *
c     iRout = 150
c     iPrint = nPrint(iRout)
c     Call qEnter('NAHSS')
*
      Call ICopy(nIrrep*16*9,[0],0,JndHss,1)
      Call iCopy(nIrrep*4*3,[0],0,JndGrd,1)
      Call LCopy(144,[.False.],0,jfHss,1)
      Call LCopy(4,[.False.],0,Tr,1)
      Call LCopy(12,[.False.],0,jfGrd,1)
*
*     COPY CNTLR MATRIXES
*
      Do  iAtom = 0, 1
         Do iCar  = 0, 2
            JfGrd(iCar,iAtom) = Ifgrd(iCar,iAtom)
            Do iIrrep=0,nIrrep-1
               JndGrd(iCar,iAtom,iIrrep)=
     &              IndGrd(iCar,iAtom,iIrrep)
            End Do
            Do  jAtom = 0, 1
               Do  jCar = 0, 2
                  JfHss(iAtom,iCar,jAtom,jCar) =
     &                 IfHss(iAtom,iCar,jAtom,jCar)
                  Do iIrrep=0,nIrrep-1
                     JndHss(iAtom,iCar,jAtom,jCar,iIrrep) =
     &                    IndHss(iAtom,iCar,jAtom,jCar,iIrrep)
                  End Do ! iirrep
               End Do !jCar
            End Do ! jAtom
         End Do !iCar
      End Do !iAtom

*
*
*-----------Derivatives with respect to the operator is computed via the
*     translational invariance.
*
      nnIrrep=nIrrep
      If (sIrrep) nnIrrep=1
      Do  iIrrep=0,nnIrrep-1
         nDisp = IndDsp(kdc,iIrrep)
         Do  iCar = 0, 2
            iComp = 2**iCar
            If (TF(kdc,iIrrep,iComp)) Then
               nDisp = nDisp + 1
*
*--------------------Reset flags for the basis set centers so that we
*     will explicitly compute the derivatives with
*     respect to those centers. Activate flag for the
*     third center so that its derivative will be comp-
*     uted by the translational invariance.
*
               JndGrd(iCar,0,iIrrep) = Abs(JndGrd(iCar,0,iIrrep))
               JndGrd(iCar,1,iIrrep) = Abs(JndGrd(iCar,1,iIrrep))
               JndGrd(iCar,2,iIrrep) = -nDisp
               JfGrd(iCar,0) = .True.
               JfGrd(iCar,1) = .True.
               JfGrd(iCar,2) = .False.
            Else
               JndGrd(iCar,2,iIrrep) = 0
            End If
        End DO
      End DO
*
*     The third center is calculated by translation invariance
*     This requires the 2nd derivatives on the other centers.
*

      Do iCar=0,2
         Do jAtom=0,2
            if (jAtom.eq.2) Then
               iStop=iCar
            Else
               iStop=2
            End If
            Do jCar=0,iStop
               Do iIrrep=0,nIrrep-1
                  If ((JndGrd(iCar,    2,iIrrep).ne.0).and.
     &                (JndGrd(jCar,jAtom,iIrrep).ne.0)) Then
                     JndHss(2,iCar,jAtom,jCar,iIrrep)=
     &                    -IX(Max(Abs(JndGrd(iCar,2,iIrrep)),
     &                    Abs(JndGrd(jCar,jAtom,iIrrep))),
     &                    Min(Abs(JndGrd(iCar,2,iIrrep)),
     &                    Abs(JndGrd(jCar,jAtom,iIrrep))))

                     Tr(2)=.true.
                     If (jAtom.eq.2) Then
                        Maxi=Max(iCar,jCar)
                        Mini=Min(iCar,jCar)
                        jfHss(0,Maxi,0,Mini)=.true.
                        jfHss(1,Maxi,1,Mini)=.true.
                        jfHss(1,iCar,0,jCar)=.true.
                        jfHss(1,jCar,0,iCar)=.true.
                     Else
                        Maxi=Max(iCar,jCar)
                        Mini=Min(iCar,jCar)
                        jfHss(jAtom,Maxi,jAtom,Mini)=.true.
                        jfHss(1,iCar,0,jCar)=.true.
                        jfHss(1,jCar,0,iCar)=.true.
                     End If ! jAtom == 2
                  End If ! if indgrd
               End Do  ! iirrep
            End Do ! jCar
         End Do ! jAtom
      End Do ! iCar
*
      IfG(0)=.true.
      IfG(1)=.true.
      IfG(2)=.false.
      IfG(3)=.false.
      Do iCent=0,1
         If (EQ(Coor(1,iCent+1),Coor(1,3) ) ) Then
            IfG(iCent)=.false.
            Do iCar=0,2
               jfGrd(iCar,iCent)=.false.
               Do kCar=0,2
                  Do KCent=0,3
                     jfHss(iCent,iCar,kCent,kCar)=.false.
                     jfHss(kCent,kCar,iCent,iCar)=.false.
                     Do iIrrep=0,nIrrep-1
                        jndHss(iCent,iCar,kCent,kCar,iIrrep)=0
                        jndHss(kCent,kCar,iCent,iCar,iIrrep)=0
                     End Do !iIrrep
                  End Do ! kcent
               End Do !kCar
               Do iIrrep=0,nIrrep-1
                  jndGrd(iCar,iCent,iIrrep)=0
               End Do !iIrrep
            End Do !ICat
         End If ! uf eq
      End Do !icent
*
      Return
      End
