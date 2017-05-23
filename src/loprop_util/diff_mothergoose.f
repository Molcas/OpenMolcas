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
      Subroutine Diff_MotherGoose(Diffuse,nAt,nB,ipMP,ipC,nij,ip_EC
     &                           ,ip_ANr,ip_Ttot
     &                           ,ip_Ttot_Inv,lMax,iTP,dLimmo
     &                           ,Thrs1,Thrs2,nThrs,iPrint
     &                           ,ThrsMul,LuYou)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"

      Dimension dLimmo(2),Pot_Expo(nij*2),Pot_Point(nij)
      Dimension Pot_Fac(nij*4)

      Logical Diffuse(3),Diffed(nij*2)

*
*-- Sag hej till publiken.
*
      Write(6,'(A)')'  Enter Slater charge distribution section.'
*
*-- Take different route for the different methods for getting
*   diffuse distributions.
*
      If(Diffuse(2)) then
        Write(6,'(A)')'    ---Run a non-linear fit,'
     &              //' (Levenberg-Marquart).'
        Write(6,'(A)')'        Thresholds'
        Write(6,991)'           Delta                   :',Thrs1
        Write(6,991)'           Lambda                  :',Thrs2
        Write(6,991)'           Factor                  :',ThrsMul
        Write(6,992)'           Min. decreasing steps   :',nThrs
        Write(6,'(A)')'        Local limit factors'
        Write(6,993)'           Low:',dLimmo(1),'     High:',dLimmo(2)
        Call Diff_Numerical(nAt,nB,ipMP,ipC,nij,Work(ip_EC)
     &                     ,iWork(ip_ANr),ip_Ttot
     &                     ,ip_Ttot_Inv,lMax,iTP,dLimmo
     &                     ,Thrs1,Thrs2,nThrs,iPrint
     &                     ,ThrsMul,Pot_Expo,Pot_Point,Pot_Fac
     &                     ,Diffed)
      Elseif(Diffuse(3)) then
        Write(6,*)
        Write(6,*)'Not programmed yet, bitte sehr.'
        Call Abend()
      Endif
991   Format(A,E12.5)
992   Format(A,I2)
993   Format(2(A,F10.5))

*
*-- Print, analyze uzw, the result of the diffuse stuff.
*
      Call WeGotThis(nAt,nB,ipMP,ipC,nij,Work(ip_EC)
     &              ,iWork(ip_ANr),ip_Ttot
     &              ,ip_Ttot_Inv,lMax,iTP,iPrint
     &              ,Pot_Expo,Pot_Point,Pot_Fac,Diffed)

*
*-- Generate file with information for other programs.
*
      lMaxF=1
      Call YouGetThis(nAt,Work(ip_EC),Pot_Expo,Pot_Point,Pot_Fac
     &               ,Diffed,ipMP,lMax,lMaxF,nij,LuYou)


      Return
      End
