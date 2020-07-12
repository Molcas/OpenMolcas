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
* Copyright (C) 1994, Anders Bernhardsson                              *
*               1994, Roland Lindh                                     *
************************************************************************
      SubRoutine CmbnT2(Rnxyz,nZeta,la,lb,Zeta,rKappa,Final,Alpha,Beta,
     &                  Hess,nHess,DAO,IfHss,IndHss,indgrd,iu,iv,nOp)
************************************************************************
*                                                                      *
* Object: compute the 2nd derivative  of the overlap matrix.           *
*                                                                      *
* Called from: OvrHss                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              DDot_   (ESSL)                                          *
*              QExit                                                   *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
c#include "print.fh"
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
      Real*8 Final(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2,6),
     &       Zeta(nZeta), rKappa(nZeta), Beta(nZeta),
     &       Rnxyz(nZeta,3,0:la+2,0:lb+2), Alpha(nZeta),
     &       Hess(nHess),
     &       DAO(nZeta,(la+1)*(la+2)/2,(lb+1)*(lb+2)/2)
      Logical IfHss(0:1,0:2,0:1,0:2)
      Integer IndHss(0:1,0:2,0:1,0:2,0:nIrrep-1),istab(0:1),
     &          nOp(2), ia(3),
     &          ib(3),indgrd(0:2,0:1,0:nirrep-1)
*
*     Statement function for Cartesian index
*
      Ind(ixyz,ix,iz) = (ixyz-ix)*(ixyz-ix+1)/2 + iz + 1
*
*     Index in the lower triang. local hessian i1 row i2 column
*
*                   itot*0  added to avoid compiler warning
      I(itot,i1,i2)=itot*0+i1*(i1-1)/2+i2
*
c     iRout = 134
      iStab(0)=iu
      iStab(1)=iv
c     iPrint = nPrint(iRout)
      iQ = 1
c     Call qEnter('CmbnT2')
*     Call GetMem(' Enter CmbnT2','LIST','REAL',iDum,iDum)
*
      exp32 = -Three/Two
      Do 25 iZeta = 1, nZeta
         rKappa(iZeta) = half*rKappa(iZeta) * Zeta(iZeta)**exp32
 25   Continue
c     If (iPrint.ge.99) Then
c        Call RecPrt(' In CmbnT2: Zeta  ',' ',Zeta  ,1,nZeta)
c        Call RecPrt(' In CmbnT2: rKappa',' ',rKappa,1,nZeta)
c        Call RecPrt(' In CmbnT2: Alpha ',' ',Alpha ,1,nZeta)
c        Call RecPrt(' In CmbnT2: Beta  ',' ',Beta  ,1,nZeta)
c     End If
      Do 10 iax = 0, la
         ia(1)=iax
         iyaMax=la-ia(1)
      Do 11 ibx = 0, lb
         ib(1)=ibx
         iybMax=lb-ib(1)
         Do 20 iay = 0, iyaMax
            ia(2)=iay
            ia(3) = la-ia(2)-ia(1)
            ipa= Ind(la,ia(1),ia(3))
         Do 21 iby = 0, iybMax
            ib(2)=iby
          ib(3) = lb-ib(2)-ib(1)
          ipb= Ind(lb,ib(1),ib(3))
*
*
*           Combine overlap integrals
*
*           Integrals like dI/dx1dx1
*
          Do 5 iCoor=1,3
            jCoor=Mod(iCoor,3)+1
            kCoor=Mod(jCoor,3)+1
            If (IfHss(0,iCoor-1,0,iCoor-1)) Then
                Do 30 iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=rKappa(iZeta)*
     &                 ((Two*Alpha(iZeta))**2 *
     &                 ((Two*Beta(iZeta))**2*
     &                 (Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor)+2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor)+2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)+2))
     &                  -Three*Two*Beta(iZeta)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  ))-
     &                  Two*Alpha(iZeta)*
     &                 ((Two*Beta(iZeta))**2*
     &                 (Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)+2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)+2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))
     &                  -Three*Two*Beta(iZeta)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )))
                  if (lb.gt.0) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-
     &                  rKappa(iZeta)*
     &                 ((Two*Alpha(iZeta))**2 *
     &                  Four*Beta(iZeta)*Dble(lb)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
     &                  -Two*Alpha(iZeta)*(Four*Beta(iZeta))*Dble(lb)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  ))
                  End If
                  if (ib(icoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+
     &                  rKappa(iZeta)*
     &                 ((Two*Alpha(iZeta))**2 *
     &                  Dble(ib(icoor)*(ib(icoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor)-2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
     &                  -Two*Alpha(iZeta)*
     &                  Dble(ib(icoor)*(ib(icoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)-2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  ))
                  End If
                  if (ib(jcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+
     &                  rKappa(iZeta)*
     &                 ((Two*Alpha(iZeta))**2 *
     &                  Dble(ib(jcoor)*(ib(jcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor)-2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
     &                  -Two*Alpha(iZeta)*
     &                  Dble(ib(jcoor)*(ib(jcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  ))
                  End If
                  if (ib(kcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+
     &                  rKappa(iZeta)*
     &                 ((Two*Alpha(iZeta))**2 *
     &                  Dble(ib(kcoor)*(ib(kcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+2,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),  ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)-2)
     &                  -Two*Alpha(iZeta)*
     &                  Dble(ib(kcoor)*(ib(kcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2))
                  End If
                  If (ia(iCoor).gt.0) Then
                        Final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) =
     &                   Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))
     &                   - rKappa(iZeta)*
     &                   Four *  Alpha(iZeta)* Dble(ia(iCoor)) *
     &                  ((Two*Beta(iZeta))**2*
     &                  (Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)+2)*
     &                   Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                   Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )+
     &                   Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                   Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)+2)*
     &                   Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )+
     &                   Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                   Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                   Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))-
     &                  (Three*Two*Beta(iZeta)*
     &                   Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)  )*
     &                   Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                   Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )))
                        If (lb.gt.0) Then
                         Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=
     &                   Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))+
     &                       rKappa(iZeta)*
     &                       Four*Beta(iZeta)*Dble(lb)*
     &                       Four *  Alpha(iZeta)* Dble(ia(iCoor)) *
     &                      (Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))*
     &                       Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))*
     &                       Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                        End If
                        If (ib(iCoor).gt.1) Then
                         Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=
     &                   Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-
     &                       rKappa(iZeta)*
     &                       Dble(ib(icoor)*(ib(icoor)-1))*
     &                       Four *  Alpha(iZeta)* Dble(ia(iCoor)) *
     &                      (Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor)-2)*
     &                       Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))*
     &                       Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                        End If
                       If (ib(jCoor).gt.1) Then
                         Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=
     &                   Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-
     &                       rKappa(iZeta)*
     &                       Dble(ib(jcoor)*(ib(jcoor)-1))*
     &                       Four *  Alpha(iZeta)* Dble(ia(iCoor)) *
     &                      (Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))*
     &                       Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)*
     &                       Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)))
                        End If
                       If (ib(kCoor).gt.1) Then
                         Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))=
     &                   Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))-
     &                       rKappa(iZeta)*
     &                       Dble(ib(kcoor)*(ib(kcoor)-1))*
     &                       Four *  Alpha(iZeta)* Dble(ia(iCoor)) *
     &                      (Rnxyz(iZeta,iCoor,ia(iCoor),ib(iCoor))*
     &                       Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))*
     &                       Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2))
                        End If
                     End If
                     If (ia(iCoor).gt.1) Then
                        Final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))
     &                     +rKappa(iZeta)*Dble(ia(iCoor)*(ia(iCoor)-1))*
     &                   ((Two*Beta(iZeta))**2*
     &                   (Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor)+2)*
     &                    Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                    Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )+
     &                    Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor)  )*
     &                    Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)+2)*
     &                    Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )+
     &                    Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor)  )*
     &                    Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                    Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)+2))-
     &                   (Three*Two*Beta(iZeta)*
     &                    Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor)  )*
     &                    Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                    Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  )))
                     If (lb.gt.0) Then
                        Final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))
     &                     -rKappa(iZeta)*Dble(ia(iCoor)*(ia(iCoor)-1))*
     &                      Four*Beta(iZeta)*Dble(lb)*
     &                   (Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))*
     &                    Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                    Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  ))
                     End If
                     If (ib(iCoor).gt.1) Then
                       Final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))
     &                      +rKappa(iZeta)*Dble(ia(iCoor)*(ia(iCoor)-1)*
     &                     ib(iCoor)*(ib(iCoor)-1))*
     &                   (Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor)-2)*
     &                    Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)  )*
     &                    Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  ))
                     End If
                     If (ib(jCoor).gt.1) Then
                       Final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))
     &                      +rKappa(iZeta)*Dble(ia(iCoor)*(ia(iCoor)-1)*
     &                     ib(jCoor)*(ib(jCoor)-1))*
     &                   (Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))*
     &                    Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor)-2)*
     &                    Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)  ))
                     End If
                     If (ib(kCoor).gt.1) Then
                       Final(iZeta,ipa,ipb,I(6,iCoor,iCoor)) =
     &                      Final(iZeta,ipa,ipb,I(6,iCoor,iCoor))
     &                      +rKappa(iZeta)*Dble(ia(iCoor)*(ia(iCoor)-1)*
     &                     ib(kCoor)*(ib(kCoor)-1))*
     &                   (Rnxyz(iZeta,iCoor,ia(iCoor)-2,ib(iCoor))*
     &                    Rnxyz(iZeta,jCoor,ia(jCoor),ib(jCoor))*
     &                    Rnxyz(iZeta,kCoor,ia(kCoor),ib(kCoor)-2))
                     End If
                   End If
 30               Continue
                  End If
 5       Continue
*
*
*           Integrals like dI/dxdz
*
         Do 57 kCoor=1,3
            iCoor=Mod(kCoor,3)+1
            jCoor=Mod(iCoor,3)+1
            iMax=Max(iCoor,jCoor)
            jCoor=Min(iCoor,jCoor)
            iCoor=iMax
            If (IfHss(0,iCoor-1,0,jCoor-1)) Then
               Do 35 iZeta = 1, nZeta
                  Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=rKappa(iZeta)*
     &                 ((Two*Alpha(iZeta))**2 *
     &                 ((Two*Beta(iZeta))**2*
     &                 (Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)+2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)+2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)+2))
     &                  -Three*Two*Beta(iZeta)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))))
                  if (lb.gt.0) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                 (Two*Alpha(iZeta))**2 *
     &                  Four*Beta(iZeta)*Dble(lb)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(icoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                 (Two*Alpha(iZeta))**2 *
     &                  Dble(ib(icoor)*(ib(icoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)-2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(jcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                 (Two*Alpha(iZeta))**2 *
     &                  Dble(ib(jcoor)*(ib(jcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)-2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(kcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                 (Two*Alpha(iZeta))**2 *
     &                  Dble(ib(kcoor)*(ib(kcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)-2)
                  End If
                  if (ia(icoor).gt.0) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Dble(ia(icoor))*Two*Alpha(iZeta)*
     &                 ((Two*Beta(iZeta))**2*
     &                 (Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)+2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)+2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)+2))
     &                  -Three*Two*Beta(iZeta)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  ))
                  if (lb.gt.0) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                  Two*Alpha(iZeta)*Dble(ia(iCoor)) *
     &                  Four*Beta(iZeta)*Dble(lb)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(icoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Two*Alpha(iZeta)*Dble(ia(iCoor) *
     &                  ib(icoor)*(ib(icoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)-2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(jcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Two*Alpha(iZeta)*Dble(ia(iCoor) *
     &                  ib(jcoor)*(ib(jcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor)-2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(kcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Two*Alpha(iZeta) *Dble(ia(iCoor)*
     &                  ib(kcoor)*(ib(kcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)+1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)-2)
                  End If
                  End If
                  if (ia(jcoor).gt.0) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Dble(ia(jcoor))*Two*Alpha(iZeta)*
     &                 ((Two*Beta(iZeta))**2*
     &                 (Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)+2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)+2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)+2))
     &                  -Three*Two*Beta(iZeta)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,  ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  ))
                  if (lb.gt.0) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                  Two*Alpha(iZeta)*Dble(ia(jCoor)) *
     &                  Four*Beta(iZeta)*Dble(lb)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(icoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Two*Alpha(iZeta)*Dble(ia(jCoor)*
     &                  ib(icoor)*(ib(icoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor)-2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(jcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Two*Alpha(iZeta)*Dble(ia(jCoor)*
     &                  ib(jcoor)*(ib(jcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)-2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(kcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Two*Alpha(iZeta)*Dble(ia(jCoor)) *
     &                  Dble(ib(kcoor)*(ib(kcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)+1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)-2)
                  End If
                  End If
                  if ((ia(iCoor).gt.0).and.(ia(jCoor).gt.0)) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                  Dble(ia(iCoor)*ia(jCoor))*
     &                 ((Two*Beta(iZeta))**2*
     &                 (Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)+2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)+2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  )+
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)+2))
     &                  -Three*Two*Beta(iZeta)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)  )*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)  )*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)  ))
                  if (lb.gt.0) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))-
     &                  rKappa(iZeta)*
     &                  Dble(ia(iCoor)*ia(jCoor)) *
     &                  Four*Beta(iZeta)*Dble(lb)*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(icoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                  Dble(ia(iCoor)*ia(jCoor)*
     &                  ib(icoor)*(ib(icoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor)-2)*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(jcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                  Dble(ia(iCoor)*ia(jCoor)*
     &                  ib(jcoor)*(ib(jcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor)-2)*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor))
                  End If
                  if (ib(kcoor).gt.1) Then
                    Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))=
     &              Final(iZeta,ipa,ipb,I(6,iCoor,jCoor))+
     &                  rKappa(iZeta)*
     &                  Dble(ia(iCoor)*ia(jCoor)*
     &                  ib(kcoor)*(ib(kcoor)-1))*
     &                  Rnxyz(iZeta,iCoor,ia(iCoor)-1,ib(iCoor))*
     &                  Rnxyz(iZeta,jCoor,ia(jCoor)-1,ib(jCoor))*
     &                  Rnxyz(iZeta,kCoor,ia(kCoor),  ib(kCoor)-2)
                  End If
                  End If

 35            Continue
            End If
 57         Continue
 21         Continue
 20         Continue
 11      Continue
 10      Continue
*
*     Trace the Hessian integrals
*
      nDAO = nZeta * (la+1)*(la+2)/2 * (lb+1)*(lb+2)/2
c     If (iPrint.ge.99) Then
c        Call RecPrt(' S(1)',' ',Final,nDAO,21)
c        Call RecPrt('   D ',' ',DAO,nDAO,1)
c     End If
      Do 90 iIrrep=0,nIrrep-1
      Do 100 iCnt=0,1
        Do 105  iCar=1,3
          Do 110 jCnt=0,1
            if (iCnt.eq.jCnt) Then
              iStop=iCar
            Else
              iStop=3
            End If
            Do 115 jCar=1,3
            If (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep).ne.0) Then
*
*              Accumulate contribution to the Hessian
*
*
*            Get the characteristics of the diff operator
*
               iCh=iEOr(2**(iCar-1)*iCnt,2**(jCar-1)*jCnt)
*
*             Get the character of the operator in the present irrep
*
               ps=DBLE(iChTbl(iIrrep,nOp(2))**
     &             (iCnt+jCnt))
*
*              Get the transf. character of the diff. operator
*
               ps = ps*DBLE(iPrmt(nOp(2),iCh))
*
*        If the over triangular diff. are needed multiply by two instead
*
               if ((iCnt.ne.jCnt).and.(iCar.eq.jCar).and.
     &               (Abs(indgrd(iCar-1,iCnt,iIrrep)).eq.
     &                iAbs(indgrd(jCar-1,jCnt,iIrrep)))) Then
                         ps=ps*Two
               End If
               iHess = Abs(IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep))
               Fact = DBLE(iStab(iCnt)*iStab(jCnt))/DBLE(nIrrep**2)
               Fact = Fact * ps
               if (IndHss(iCnt,iCar-1,jCnt,jCar-1,iIrrep).gt.0) Then
                oj=DDot_(nDAO,DAO,1,Final(1,1,1,I(6,iCar,jCar)),1)
                Hess(iHess) = Hess(iHess) + Fact*oj
               Else
                Fact=Fact*DBLE((-1)**(icnt+jcnt))
                oj=DDot_(nDAO,DAO,1,Final(1,1,1,I(6,Max(iCar,jCar),
     &                 Min(iCar,jCar))),1)
                Hess(iHess) = Hess(iHess) + Fact*oj
               End If
            End If
 115       Continue
 110     Continue
 105    Continue
 100  Continue
 90   Continue
*
c     Call qExit('CmbnT2')
      Return
      End
