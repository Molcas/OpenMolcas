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
      Subroutine Translation(ifg,jfgrd,jfhss,tr,jndgrd,jndhss,coorm,
     &                       nirrep,indgrd,indhss)
      Implicit Integer(a-z)
      Logical jfhss(4,3,4,3),jfgrd(3,4),tr(4),ifg(4),eq,alike

      Integer jndgrd(3,4,0:nirrep-1),jndhss(4,3,4,3,0:nirrep-1)
      Integer indgrd(3,4,0:nirrep-1),indhss(4,3,4,3,0:nirrep-1)
      Real*8 Coorm(3,4)
      Alike=.false.
      If (IfG(1).and.IfG(2).and.IfG(3).and.IfG(4)) Then
        Do 3333 iCent=1,3
          If (.not.Alike) Then
            Do 3334  jCent=iCent+1,4
               If (EQ(CoorM(1,iCent),CoorM(1,jCent))) Then
                  Do kCent=1,4
                     iMax=Max(jCent,kCent)
                     iMin=Min(jCent,kCent)
                     Do kCar=1,3
                        Do lCar=1,3
                           Do mIrrep=0,nIrrep-1
                              JndHss(iMax,kCar,iMin,lCar,mIrrep)=0
                           End Do
                           JfHss(iMax,kCar,iMin,lCar)=.false.
                        End Do
                     End Do
                  End Do
                  Do jCar=1,3
                     Do mIrrep=0,nIrrep-1
                        JndGrd(jCar,jCent,mIrrep)=0
                     End Do
                     jfGrd(jCar,jCent)=.false.
                  End Do
                  IfG(jCent)=.false.
                  If (.not.Alike) Then
                     IfG(iCent)=.false.
                     Tr(iCent)=.true.
                     Alike=.true.
                     Do kCent=1,4
                        If ((.not.(EQ(CoorM(1,iCent),
     &                       CoorM(1,kCent)))) .or.
     &                       (kCent.eq.iCent))
     &                       then
                           iMax=Max(iCent,kCent)
                           iMin=Min(iCent,kCent)
                           Do kCar=1,3
                             If (iMax.eq.iMin) Then
                                iStop=kCar
                             Else
                                iStop=3
                             End If
                             Do lCar=1,iStop
                               Do mIrrep=0,nIrrep-1
                                 JndHss(iMax,kCar,iMin,lCar,mIrrep)=
     &                            -IndHss(iMax,kCar,iMin,lCar,mIrrep)
                               End Do
                               JfHss(iMax,kCar,iMin,lCar)=.false.

*     Set the derivatives that are needed for the translation
*     invarians calculations.
*
                                    Do iiCent=1,4
                                    If (.not.EQ(CoorM(1,iiCent),
     &                                   CoorM(1,iCent)))
     &                                 Then
                                       Do jjCent=1,iiCent
                                         If (.not.EQ(CoorM(1,jjCent),
     &                                        CoorM(1,iCent)))
     &                                   Then
                                           Do kkCar=1,3
                                           if (iiCent.eq.jjCent) Then
                                              iStop=kkCar
                                           Else
                                               iStop=3
                                           End If
                                           Do llCar=1,iStop
                                               JfHss(iiCent,kkCar,
     &                                          jjCent,llCar)=.true.
                                           End Do
                                        End Do
                                    End If
                                   End Do !  icent
                                   Do kkCar=1,3
                                      JfGrd(kkCar,iiCent)=.true.
                                   End Do
                                   End If
                                 End Do
                                End Do
                              End Do
                           End If
                        End Do
                        Do jCar=1,3
                           Do mIrrep=0,nIrrep-1
                              JndGrd(jCar,iCent,mIrrep)=
     &                             -IndGrd(jCar,iCent,mIrrep)
                           End Do
                           JfGrd(jCar,iCent)=.false.
                        End  Do
                     End If
                  End If
 3334          Continue
            End If
 3333    Continue
      End If
*
*     If all centers are different delete the fourth center
*
      If (.not.Alike) Then
         IfG(4)=.false.
         Tr(4)=.true.
         Call lCopy(12,[.true.],0,JFGRD,1)
*     Call lCopy(144,[.true.],0,JFHss,1)
         Do iiC=1,4
            Do jjC=1,iiC
               Do iiCar=1,3
                  iStop=3
                  If (iic.eq.jjc) iStop=iiCar
                  Do jjCar=1,iStop
                     JfHss(iiC,iiCar,jjc,jjCar)=.true.
                  End Do
               End Do
            End Do
         End Do

         Do kCar=1,3
            Do lCent=1,4
               iStop=3
               If (lCent.eq.4) iStop=kCar
               Do lCar=1,iStop
                  Do mIrrep=0,nirrep-1
                     JndHss(4,kCar,lCent,lCar,mIrrep)=
     &                    -IndHss(4,kCar,lCent,lCar,mIrrep)
                  End Do
                  JfHss(4,kCar,lCent,lCar)=.false.
               End Do
            End Do
         End Do
*
         Do lCar=1,3
            Do mIrrep=0,nirrep-1
               Jndgrd(lCar,4,mIrrep)=
     &              -IndGrd(lCar,4,mIrrep)
            End Do
            jfgrd(lCar,4)=.false.
         End Do
*
      End If
      return
      end
