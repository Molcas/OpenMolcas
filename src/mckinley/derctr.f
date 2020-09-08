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
* Copyright (C) 1995, Anders Bernhardsson                              *
************************************************************************
      SubRoutine DerCtr(mdci,mdcj,mdck,mdcl,ldot,
     &                  JfGrd,IndGrd,JfHss,IndHss,JfG,mbatch)
*                                                                      *
************************************************************************
*                                                                      *
      use Center_Info
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "info.fh"
#include "disp.fh"
#include "disp2.fh"
      Logical JfHss(4,3,4,3),IfHss(4,3,4,3),JfGrd(3,4),IfGrd(3,4),
     &        TF,TstFnc,IfG(4),JfG(4),ldot
      Integer IndHss(4,3,4,3,0:7),JndHss(4,3,4,3,0:7),
     &        IndGrd(3,4,0:7),JndGrd(3,4,0:7),iCo(4)
*define _OLD_CODE_
#ifdef _OLD_CODE_
      Integer iCom(0:7,0:7),iStabM(0:7), idcrr(0:7)
      Logical chck
#endif
*
      Ind(i1,i2)=i1*(i1-1)/2+i2
      TF(mdc,iIrrep,iComp) = TstFnc(iOper,nIrrep,dc(mdc)%iCoSet,
     &                       iChTbl,iIrrep,iComp,dc(mdc)%nStab)
*
       nnIrrep=nIrrep
       Call lCopy(12,[.false.],0,ifgrd,1)
       If (sIrrep) nnIrrep=1
*
      Do 3000 iIrrep=0,nnIrrep-1
         nDisp = IndDsp(mdci,iIrrep)
         Do 3001 iCar = 0, 2
            iComp = 2**iCar
            If (TF(mdci,iIrrep,iComp)) Then
               nDisp = nDisp + 1
               IndGrd(iCar+1,1,iIrrep) = nDisp
               IfGrd(iCar+1,1) = .True.
            Else
               IndGrd(iCar+1,1,iIrrep) = 0
            End If
 3001    Continue
 3000 Continue
      Do 3100 iIrrep=0,nnIrrep-1
         nDisp = IndDsp(mdcj,iIrrep)
         Do 3101 iCar = 0, 2
            iComp = 2**iCar
            If (TF(mdcj,iIrrep,iComp)) Then
               nDisp = nDisp + 1
               IndGrd(iCar+1,2,iIrrep) = nDisp
               IfGrd(iCar+1,2) = .True.
            Else
               IndGrd(iCar+1,2,iIrrep) = 0
            End If
 3101    Continue
 3100 Continue
      Do 3200 iIrrep=0,nnIrrep-1
         nDisp = IndDsp(mdck,iIrrep)
         Do 3201 iCar = 0, 2
            iComp = 2**iCar
            If (TF(mdck,iIrrep,iComp)) Then
               nDisp = nDisp + 1
               IndGrd(iCar+1,3,iIrrep) = nDisp
               IfGrd(iCar+1,3) = .True.
            Else
               IndGrd(iCar+1,3,iIrrep) = 0
            End If
 3201    Continue
 3200  Continue
       Do 3300 iIrrep=0,nnIrrep-1
          nDisp = IndDsp(mdcl,iIrrep)
          Do 3301 iCar = 0, 2
             iComp = 2**iCar
             If (TF(mdcl,iIrrep,iComp)) Then
                nDisp = nDisp + 1
                IndGrd(iCar+1,4,iIrrep) = nDisp
                IfGrd(iCar+1,4) = .True.
             Else
                IndGrd(iCar+1,4,iIrrep) = 0
             End If
 3301    Continue
 3300 Continue
      Do iIrrep=0,nnIrrep-1
         Do 3333 iCar = 1, 3
            Do 4444 iSh = 1, 4
               JndGrd(iCar,iSh,iIrrep) =
     &            IndGrd(iCar,iSh,iIrrep)
               JfGrd(iCar,iSh) =
     &            IfGrd(iCar,iSh)
 4444        Continue
 3333     Continue
      End Do
      iCo(1)=mdci
      iCo(2)=mdcj
      iCo(3)=mdck
      iCo(4)=mdcl
      Call iCopy(144*nirrep,[0],0,IndHss,1)
      Call iCopy(144*nirrep,[0],0,jndHss,1)
      Call lCopy(144,[.false.],0,IfHss,1)
      Call lCopy(144,[.false.],0,JfHss,1)
      if (.not.ldot) Return
*
      Do iAtom=1,4
         Do jAtom=1,iAtom
*
*        This segment of the code is not really needed.
*        If turned on it should not do much of a difference.
*
#ifdef _OLD_CODE_
            Call DCR(LmbdR,iOper,nIrrep,dc(iCo(iAtom))%iStab,
     &               dc(iCo(iAtom))%nStab,dc(iCo(jAtom))%iStab,
     &               dc(iCo(jAtom))%nStab,iDCRR,nDCRR)
*
*-----------Find the stabilizer for A and B
*
            Call Inter(dc(iCo(iAtom))%iStab,dc(iCo(iAtom))%nStab,
     &                 dc(iCo(jAtom))%iStab,dc(iCo(jAtom))%nStab,
     &                 iStabM,nStabM)
*
*          Generate all possible (left) CoSet
*          To the stabil. of A and B
*
            Do iIrrep = 0, nIrrep-1
               Do jOper = 0, nStabM-1
                  iCoM(iIrrep,jOper) =
     &                 iEor(iOper(iIrrep),iStabM(jOper))
               End Do
            End Do
*
*           Order the Coset so we will have the unique ones first
*           Check uniqueness
*
            nMax = 1
            Do 435 j = 1, nIrrep-1
               Do 436 i = 0, nMax - 1
                  Do 437 ielem = 0, nStabM-1
                    If (iCoM(i,1).eq.iCoM(j,ielem))
     &                  Go To 435
 437              Continue
 436           Continue
*
*----------Move unique CoSet
*
               nMax = nMax + 1
               Do 438 ielem = 0, nStabM-1
                  iTmp = iCoM(nMax-1,ielem)
                  iCoM(nMax-1,ielem) = iCoM(j,ielem)
                  iCoM(j,ielem) = iTmp
 438           Continue
               If (nMax.eq.nIrrep/nStabM) Go To 439
 435        Continue
 439        Continue
*
*           Check if the derivative is needed in the present symmetry
*
            nCoM=nIrrep/nStabM
*
            Do iCar=1,3
               if (iAtom.eq.jAtom) Then
                   istop=iCar
               Else
                   iStop=3
               End If
               Do jCar=1,istop
                  iComp=iEOr(2**(iCar-1),2**(jCar-1))
                  Chck=TstFnc(iOper,nIrrep,iCoM,iChTbl,0,iComp,nStabM)
                  If (Chck)
     &                 IfHss(iAtom,iCar,jAtom,jCar)=.true.
               End Do
            End Do
#endif
*
*               Calculate the index for the derivative
*
            Do  iIrrep=0,nnIrrep-1
                Do iCar=1,3
                   if (iAtom.eq.jAtom) Then
                      istop=iCar
                   Else
                      iStop=3
                   End If
                   Do jCar=1,istop
            IfHss(iAtom,iCar,jAtom,jCar)=.true.
                        If ((jndGrd(iCar,iAtom,iIrrep).gt.0).and.
     &                     (jndGrd(jCar,jAtom,iIrrep).gt.0))
     &                  Then
                            IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=
     &                      Ind(Max(JndGrd(iCar,iAtom,iIrrep),
     &                          jndGrd(jCar,jAtom,iIrrep)),
     &                      Min(jndGrd(iCar,iAtom,iIrrep),
     &                      jndGrd(jCar,jAtom,iIrrep)))
                         Else
                            IndHss(iAtom,iCar,jAtom,jCar,iIrrep)=0
                         End If
                   End Do
                End Do
            End Do
          End Do
      End Do
*
*              Scramble the control array for the hessian
*
      Call LCopy(4,[.true.],0,Ifg,1)
      Do iAtom=1,4
         JfG(iAtom)=IfG(iAtom)
         Do jAtom=1,iAtom
            Do iCar=1,3
              If (iAtom.eq.jAtom) Then
                 iStop=iCar
              Else
                 iStop=3
              End If
              Do jCar=1,iStop
                 If (iAtom.ge.jAtom) Then
                    JfHss(iAtom,iCar,jAtom,jCar)=
     &              IfHss(iAtom,iCar,jAtom,jCar)
                 Else If (iAtom.lt.jAtom) Then
                    JfHss(iAtom,iCar,jAtom,jCar)=
     &              IfHss(jAtom,jCar,iAtom,iCar)
                 End If
              End Do
           End Do
         End Do
       End Do
       If (sIrrep) Then
        Do ii=1,4
         Do ic1=1,3
          If (indgrd(ic1,ii,0).eq.0)
     &       jfgrd(ic1,ii)=.false.
          Do ij=1,4
           Do ic2=1,3
            If (Indhss(ii,ic1,ij,ic2,0).eq.0)
     &           JfHss(ii,ic1,ij,ic2)=.false.
           End Do
          End Do
         End Do
        End Do
       End If
*----------------------------------------------------------------
*
*  End Hess
*
*----------------------------------------------------------------------*
       Return
c Avoid unused argument warnings
       If (.False.) Call Unused_integer(mbatch)
       End
