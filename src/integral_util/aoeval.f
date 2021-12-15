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
* Copyright (C) 1995, Roland Lindh                                     *
*               2004, Takashi Tsuchiya                                 *
************************************************************************
      Subroutine AOEval(iAng,nCoor,Coor,xyz,RA,Transf,CffSph,nElem,nCmp,
     &                  Angular,nTerm,nForm,Thr_Rad,nRad,mExp,nExp,
     &                  Alpha,Radial,nBas,CffCnt,AOValue,mAO,
     &                  px,py,pz,ipx,ipy,ipz)
************************************************************************
* Object: to compute the values of the AOs on a grid. The calculation  *
*         is sectioned in an angular part and a radial part. Due to the*
*         gradients we have 3 angular parts and 2 radial part.         *
*                                                                      *
*      Author:Roland Lindh, Dept. of Theoretical Chemistry, University *
*             of Lund, SWEDEN. November '95                            *
*      Modified by: Takashi Tsuchiya, Dept. of Theoretical Chemistry,  *
*                   University of Lund, SWEDEN. February 2004          *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "print.fh"
      Real*8 xyz(nCoor,3,0:iAng+nRad-1), Coor(3,nCoor), RA(3),
     &       CffSph(nElem,nCmp), Alpha(nExp), Radial(nCoor,nRad,nBas),
     &       CffCnt(mExp,nBas), AOValue(mAO,nCoor,nBas,nCmp)
      Integer Angular(nTerm,5,nForm)
      Logical Transf
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Character*80 Label
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Index for derivatives in AOValue
*
*     mAO =   1:    F
*             2:   dFdx
*             3:   dFdy
*             4:   dFdz
*             5:  d2Fdxdx
*             6:  d2Fdxdy
*             7:  d2Fdxdz
*             8:  d2Fdydy
*             9:  d2Fdydz
*            10:  d2Fdzdz
*            11:  d3Fdxdxdx
*            12:  d3Fdxdxdy
*            13:  d3Fdxdxdz
*            14:  d3Fdxdydy
*            15:  d3Fdxdydz
*            16:  d3Fdxdzdz
*            17:  d3Fdydydy
*            18:  d3Fdydydz
*            19:  d3Fdydzdz
*            20:  d3Fdzdzdz
*            ... and so on;
*                in the same order for the higher order derivatives.
*                                                                      *
************************************************************************
*                                                                      *
      Ind(iy,iz) = (iy+iz)*(iy+iz+1)/2 + iz + 1
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      iRout=132
      iPrint=nPrint(iRout)
C     iPrint=99
      If (iPrint.ge.49) Then
         Write(6,*) '********** AOEval ***********'
         Write (6,*) 'In AOEval'
         Call RecPrt('Coor',' ',Coor,3,nCoor)
         Call RecPrt('CffCnt',' ',CffCnt,mExp,nBas)
         Call RecPrt('CffSph',' ',CffSph,nElem,nCmp)
         Write (6,*) 'RA=',RA
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*---- Initialize AOValue
*
      Call FZero(AOValue,mAO*nCoor*nBas*nCmp)
*
*---- Set the order of derivation
*
      nDrv=nRad-1
C     Write(6,*) '----- nDrv = ', nDrv
*                                                                      *
************************************************************************
*                                                                      *
*---- Evaluate the radial part
*----    Normal radial part and
*----    premultiplied with (minus two times the exponent)**iDrv
*
#ifdef _DEBUGPRINT_
      If (nRad.LE.0 .AND. nRad.GE.5) Then
         Write (6,*) 'AOEval: illegal value of nRad!'
         Call Abend()
      End If
#endif
*
      Thre=-99D0
      If (Thr_Rad.gt.0.0D0) Thre=Log(Thr_Rad)
      Exp_Min=1.0D10
      Do iExp = 1, nExp
         Exp_Min=Min(Exp_Min,Alpha(iExp))
      End Do
      Call FZero(Radial,nCoor*nRad*nBas)
      Do iCoor = 1, nCoor
         R2=(Coor(1,iCoor)-RA(1))**2
     &     +(Coor(2,iCoor)-RA(2))**2
     &     +(Coor(3,iCoor)-RA(3))**2
C        If (-Exp_Min*R2.lt.Thre) Go To 9898
         Do iExp = 1, nExp
            If (-Alpha(iExp)*R2.lt.Thre) Go To 9898
            Tmp=Exp(-Alpha(iExp)*R2)
            If (nRad.eq.1) Then
               Do iBas = 1, nBas
                  XCff=CffCnt(iExp,iBas)
                  Radial(iCoor,1,iBas)=Radial(iCoor,1,iBas) + XCff*Tmp
               End Do
            Else If (nRad.eq.2) Then
               Tmp2=-Two*Alpha(iExp)*Tmp
               Do iBas = 1, nBas
                  XCff=CffCnt(iExp,iBas)
                  Radial(iCoor,1,iBas)=Radial(iCoor,1,iBas) + XCff*Tmp
                  Radial(iCoor,2,iBas)=Radial(iCoor,2,iBas) + XCff*Tmp2
               End Do
            Else If (nRad.eq.3) Then
               Tmp2=-Two*Alpha(iExp)*Tmp
               Tmp3=-Two*Alpha(iExp)*Tmp2
               Do iBas = 1, nBas
                  XCff=CffCnt(iExp,iBas)
                  Radial(iCoor,1,iBas)=Radial(iCoor,1,iBas) + XCff*Tmp
                  Radial(iCoor,2,iBas)=Radial(iCoor,2,iBas) + XCff*Tmp2
                  Radial(iCoor,3,iBas)=Radial(iCoor,3,iBas) + XCff*Tmp3
               End Do
            Else If (nRad.eq.4) Then
               Tmp2=-Two*Alpha(iExp)*Tmp
               Tmp3=-Two*Alpha(iExp)*Tmp2
               Tmp4=-Two*Alpha(iExp)*Tmp3
               Do iBas = 1, nBas
                  XCff=CffCnt(iExp,iBas)
                  Radial(iCoor,1,iBas)=Radial(iCoor,1,iBas) + XCff*Tmp
                  Radial(iCoor,2,iBas)=Radial(iCoor,2,iBas) + XCff*Tmp2
                  Radial(iCoor,3,iBas)=Radial(iCoor,3,iBas) + XCff*Tmp3
                  Radial(iCoor,4,iBas)=Radial(iCoor,4,iBas) + XCff*Tmp4
               End Do
            Else
               Do iBas = 1, nBas
                  XCff=CffCnt(iExp,iBas)
                  Radial(iCoor,1,iBas)=Radial(iCoor,1,iBas) + XCff*Tmp
               End Do
               Do iDrv = 1, nDrv
                  Tmp=-Two*Alpha(iExp)*Tmp
                  Do iBas = 1, nBas
                     XCff=CffCnt(iExp,iBas)
                     Radial(iCoor,iDrv+1,iBas)=Radial(iCoor,iDrv+1,iBas)
     &                                        +  XCff*Tmp
                  End Do
               End Do
            End If
         End Do
 9898    Continue
      End Do
*
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
            Write (6,*) mExp,nExp
            Write (Label,'(A)')'Radial(nCoor*nRad,nBas)'
            Call RecPrt(Label,'(10G20.10)',Radial,nCoor*nRad,nBas)
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Evaluate the angular part
*
      If (iAng+nRad-1.ge.1) Then
         Do iCoor = 1, nCoor
            xyz(iCoor,1,0)=One
            xyz(iCoor,2,0)=One
            xyz(iCoor,3,0)=One
            xyz(iCoor,1,1)=px*(Coor(1,iCoor)-RA(1))
            xyz(iCoor,2,1)=py*(Coor(2,iCoor)-RA(2))
            xyz(iCoor,3,1)=pz*(Coor(3,iCoor)-RA(3))
         End Do
         Do i = 2, iAng+nRad-1
            Do iCoor = 1, nCoor
               xyz(iCoor,1,i)=xyz(iCoor,1,i-1)*xyz(iCoor,1,1)
               xyz(iCoor,2,i)=xyz(iCoor,2,i-1)*xyz(iCoor,2,1)
               xyz(iCoor,3,i)=xyz(iCoor,3,i-1)*xyz(iCoor,3,1)
            End Do
         End Do
      Else
         call dcopy_(3*nCoor,[One],0,xyz(1,1,0),1)
      End If
#ifdef _DEBUGPRINT_
      If (iPrint.ge.99) Then
         Do i = 0, iAng+nRad-1
            Write (Label,'(A,I2,A)')'xyz(nCoor,nCar,',i,')'
            Call RecPrt(Label,'(3E12.6)',xyz(1,1,i),nCoor,3)
         End Do
      End If
#endif
*                                                                      *
************************************************************************
*                                                                      *
*-----Calculate the angular components of the derivatives
*
      Call ICopy(5*nForm*nTerm,[0],0,Angular,1)
*
      Do ix = iAng, 0, -1
*
         Do iy = iAng-ix, 0, -1
            iz = iAng-ix-iy
            ip=Ind(iy,iz)
*
*---------- Initial values for 0th-order derivative
*
            Angular(1,1,1)=ix
            Angular(1,2,1)=iy
            Angular(1,3,1)=iz
            Angular(1,4,1)=0
            Angular(1,5,1)=1
*
*
*---------- Calculate derivatives of the angular components
*
*           jf: Formula to be derivated
*           if: New formula by way of derivation
*
            if = 1
C           Call PrntA(nterm,nform,Angular,if,0)
            Do jDrv = 0,nDrv-1
               Do jx = jDrv,0,-1
                  Do jy = jDrv-jx,0,-1
                     jz = jDrv-jx-jy
                     jf = Ind(jy,jz)
*
                     Do kDrv = 0,jDrv-1
                        jf=jf+(kDrv+1)*(kDrv+2)/2
                     End Do
*
#ifdef _DEBUGPRINT_
                     If (iPrint.ge.99) Then
                        Write (6,*) ' jx,jy,jz,jf=',jx,jy,jz,jf
                     End If
#endif
*
                     If (jy .EQ. 0 .AND. jz .EQ. 0) then
                       if=if+1
                       Call dFdxyz(nTerm,nForm,Angular,jf,if,1,ipx,jDrv)
C                      Call PrntA(nterm,nform,Angular,if,jdrv+1)
                       if=if+1
                       Call dFdxyz(nTerm,nForm,Angular,jf,if,2,ipy,jDrv)
C                      Call PrntA(nterm,nform,Angular,if,jdrv+1)
                       if=if+1
                       Call dFdxyz(nTerm,nForm,Angular,jf,if,3,ipz,jDrv)
C                      Call PrntA(nterm,nform,Angular,if,jdrv+1)
                     Else if (jz .EQ. 0) then
                       if=if+1
                       Call dFdxyz(nTerm,nForm,Angular,jf,if,2,ipy,jDrv)
C                      Call PrntA(nterm,nform,Angular,if,jdrv+1)
                       if=if+1
                       Call dFdxyz(nTerm,nForm,Angular,jf,if,3,ipz,jDrv)
C                      Call PrntA(nterm,nform,Angular,if,jdrv+1)
                     Else
                       if=if+1
                       Call dFdxyz(nTerm,nForm,Angular,jf,if,3,ipz,jDrv)
C                      Call PrntA(nterm,nform,Angular,if,jdrv+1)
                     Endif
                  End Do
               End Do
            End Do
*
*
*---------- Distribute contributions to the real spherical harmonics
*----------            and combine with angular part
*----------       to yield the values and the gradients
*
            If (Transf) Then
              Do iCmp = 1, nCmp
                Cff=CffSph(ip,iCmp)
                If (Cff.ne.Zero) Then
                  kForm=0
                  Do iDrv = 0,nDrv
                    mForm=(iDrv+1)*(iDrv+2)/2
                    Do iForm = kForm+1, kForm+mForm
                      mTerm=2**(iDrv)
                      Do iTerm = 1, mTerm
                        iCoef=Angular(iTerm,5,iForm)
                        If (iCoef.ne.0) Then
                          Coef=DBLE(iCoef)
                          iRad =Angular(iTerm,4,iForm)+1
                          Do iBas = 1, nBas
                            Do iCoor = 1, nCoor
                              AOValue(iForm,iCoor,iBas,iCmp)
     &                         = AOValue(iForm,iCoor,iBas,iCmp)
     &                         + xyz(iCoor,1,Angular(iTerm,1,iForm))
     &                         * xyz(iCoor,2,Angular(iTerm,2,iForm))
     &                         * xyz(iCoor,3,Angular(iTerm,3,iForm))
     &                         * Coef
     &                         * Cff
     &                         * Radial(iCoor,iRad,iBas)
                            End Do
                          End Do
                        End If
                      End Do
                    End Do
                    kForm=kForm+mForm
                  End Do
                End If
              End Do
            Else
              kForm=0
              Do iDrv = 0, nDrv
                mForm=(iDrv+1)*(iDrv+2)/2
                Do iForm = kForm+1, kForm+mForm
                  mTerm=2**(iDrv)
                  Do iTerm = 1, mTerm
                    iCoef=Angular(iTerm,5,iForm)
                    If (iCoef.ne.0) Then
                      Coef=DBLE(iCoef)
                      iRad =Angular(iTerm,4,iForm)+1
                      Do iBas = 1, nBas
                        Do iCoor = 1, nCoor
                          AOValue(iForm,iCoor,iBas,ip)
     &                     = AOValue(iForm,iCoor,iBas,ip)
     &                     + xyz(iCoor,1,Angular(iTerm,1,iForm))
     &                     * xyz(iCoor,2,Angular(iTerm,2,iForm))
     &                     * xyz(iCoor,3,Angular(iTerm,3,iForm))
     &                     * Coef
     &                     * Radial(iCoor,iRad,iBas)
                        End Do
                      End Do
                    End If
                  End Do
                End Do
                kForm=kForm+mForm
              End Do
            End If
*
*
         End Do
*
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      If (iPrint.ge.49) Then
         Write(6,*) 'mAO,nCoor,nBas,nCmp',mAO,nCoor,nBas,nCmp
         Call RecPrt('AOValue','(10G20.10)',AOValue,mAO,nCoor*nBas*nCmp)
      End If
#endif

      End
*
      Subroutine dFdxyz(mterm,mform,N,jp,ip,ixyz,ipf,jdrv)
*
      Implicit real*8 (a-h,o-z)
      Integer N(mterm,5,mform)
*
*
*     ipf: Phase factor in integer
*
*
*
*
      nterm=2**jdrv
      iterm=0
*
      Do jterm=1,nterm
*
*
*        downward operation
*        (derivation of angular part)
*
*
         iterm=iterm+1
         Do i=1,5
            If (i .EQ. ixyz) then
               N(iterm,ixyz,ip)=N(jterm,ixyz,jp)-1
            Else
               N(iterm,i,ip)=N(jterm,i,jp)
            End if
         End do
         N(iterm,5,ip)=N(iterm,5,ip)*N(jterm,ixyz,jp)
         N(iterm,5,ip)=N(iterm,5,ip)*ipf
*
*
*        upward operation
*        (derivation of radial  part)
*
*
         iterm=iterm+1
         Do i=1,5
            If (i .EQ. ixyz) then
               N(iterm,ixyz,ip)=N(jterm,ixyz,jp)+1
            Else
               N(iterm,i,ip)=N(jterm,i,jp)
            End if
         End do
         N(iterm,4,ip)=N(iterm,4,ip)+1
         N(iterm,5,ip)=N(iterm,5,ip)*ipf
      End do
      End
