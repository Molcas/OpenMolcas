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
************************************************************************
*  RdInput
*
*> @brief
*>   Reads the input of the quater program
*> @author Y. Carissan
*>
*> @param[out] U1 vector needed for the rotation
*> @param[out] U2 vector needed for the rotation
*> @param[out] V1 vector needed for the rotation
*> @param[out] V2 vector needed for the rotation
*>
*> @details
*> Reads the input of the quater program.
************************************************************************
      Subroutine RdInput(U1,U2,V1,V2)
      Implicit none
#include "WrkSpc.fh"
#include "debug.fh"
#include "geoms.fh"
#include "options.fh"
      Integer iLU
      integer ig,iat
      Real*8 U1(3),U2(3),V1(3),V2(3)
      Character*180 Get_Ln
      character*6 cName
      External Get_Ln
      Character Line*180
      Character key*20
      Logical AxisSet,NewAxisSet,Geo1Set,Geo2Set
      Logical XYZ1Set,XYZ2Set

      AxisSet=.false.
      NewAxisSet=.false.
      Geo1Set=.false.
      Geo2Set=.false.
      XYZ1Set=.false.
      XYZ2Set=.false.
      rotate=.true.
      translate=.true.
      ngeoms=1
      Call SpoolInp(iLU)
c
c skip &quater &end
c
      Line=Get_ln(iLU)
666   continue
c      Do
        Line=Get_ln(iLU)
        Call Put_ln(Line)
        Call Get_S(1,key,1)
        if (debug) Write(6,*) 'KEY:',key
        if (key(1:4).eq.'AXIS') then
          Line=Get_ln(iLU)
          Call Put_ln(Line)
          Call Get_F(1,U1,3)
          Line=Get_ln(iLU)
          Call Put_ln(Line)
          Call Get_F(1,U2,3)
          AxisSet=.true.
        else if (key(1:4).eq.'NEWA') then
          Line=Get_ln(iLU)
          Call Put_ln(Line)
          Call Get_F(1,V1,3)
          Line=Get_ln(iLU)
          Call Put_ln(Line)
          Call Get_F(1,V2,3)
          NewAxisSet=.true.
        else if (key(1:4).eq.'DEBU') then
          debug=.true.
        else if (key(1:4).eq.'GEO1') then
          call readgeo(iLU,1)
          Geo1Set=.true.
        else if (key(1:4).eq.'GEO2') then
          call readgeo(iLU,2)
          Geo2Set=.true.
        else if (key(1:4).eq.'XYZ1') then
          Line=Get_ln(iLU)
          Call Put_ln(Line)
          Call Get_I(1,XYZ1,3)
          XYZ1Set=.true.
        else if (key(1:4).eq.'XYZ2') then
          Line=Get_ln(iLU)
          Call Put_ln(Line)
          Call Get_I(1,XYZ2,3)
          XYZ2Set=.true.
!        else if (key(1:4).eq.'NGEO') then
!          Line=Get_ln(iLU)
!          Call Put_ln(Line)
!          Call Get_I(1,ngeoms,1)
        else if (key(1:4).eq.'NOTR') then
          translate=.false.
        else if (key(1:4).eq.'NORO') then
          rotate=.false.
        else if (key(1:4).eq.'END ') then
          GoTo 999
        else
          Call SysAbendMsg("RdInput",
     &          "Keyword not relevant : ",key)
        endif
      goto 666

999   continue
      if (.not.AxisSet) then
        if (.not.XYZ1Set) then
          Call SysAbendMsg("RdInput",
     &      "Reference Axis not set",
     &      "AXIS or XYZ1 Keyword mandatory")
        end if
      else
        if (XYZ1Set) then
          Call SysAbendMsg("RdInput",
     &      "Reference Axis not set properly",
     &      "AXIS and XYZ1 Keywords are exclusive")
        end if
      end if

      if (.not.NewAxisSet) then
        if (.not.XYZ2Set) then
          Call SysAbendMsg("RdInput",
     &      "New Axis not set :","NEWAXIS or XYZ2 Keyword mandatory")
        end if
      else
        if (XYZ2Set) then
          Call SysAbendMsg("RdInput",
     &      "New Axis not set properly",
     &      "NEWAXIS and XYZ2 Keywords are exclusive")
        end if
      end if

      if (XYZ1Set.and..not.GEO1Set) then
          Call SysAbendMsg("RdInput",
     &      "XYZ1 keyword requires GEO1 definition",
     &      "")
      end if

      if (XYZ2Set.and..not.GEO2Set) then
          Call SysAbendMsg("RdInput",
     &      "XYZ2 keyword requires GEO2 definition",
     &      "")
      end if

      if (translate.and..not.(XYZ1Set.and.XYZ2Set)) then
          Call SysAbendMsg("RdInput",
     &      "Translation cannot be done if both",
     &      "XYZ1 and XYZ2 are not set")
      end if

      if (XYZ1Set) Call SetVect(nat(1),Work(ipgeo(1)),XYZ1,V1,V2)
      if (XYZ2Set) Call SetVect(nat(2),Work(ipgeo(2)),XYZ2,U1,U2)

      Call Getmem("TEST","List","Real",ig,ig)

      do ig=3,ngeoms+2
        if (ig.lt.10) then
          Write(cName,'(a5,i1)') "GEOM0",ig
          title(ig)=cName
        else if (ig.lt.100) then
          Write(cName,'(a4,i2)') "GEOM",ig
          title(ig)=cName
        end if
        nat(ig)=nat(2)
        do iat=1,nat(ig)
          geolbl(iat,ig)=geolbl(iat,2)
        end do
        Call GetMem(cName,"Allo","Real",ipgeo(ig),3*nat(ig))
      end do

      End
