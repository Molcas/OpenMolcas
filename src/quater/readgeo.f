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
      subroutine readgeo(iLU,ig)
************************************************************
*
*   <DOC>
*     <Name>readgeo</Name>
*     <Syntax>readgeo(iLU,ig)</Syntax>
*     <Arguments>
*       \Argument{iLU}{logic unit number}{Integer}{in}
*       \Argument{ig}{geometry index}{Integer}{in}
*     </Arguments>
*     <Purpose>Read a geometry and stores it into memory</Purpose>
*     <Dependencies>Get Ln and utilities</Dependencies>
*     <Author>Y. Carissan</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects>none</Side_Effects>
*     <Description>
*        Raed a geometry in the XYZ format and stores it
*        in memory at address Work(ipgeo(ig)).
*     </Description>
*    </DOC>
*
************************************************************
      implicit none
#include "WrkSpc.fh"
#include "debug.fh"
#include "geoms.fh"
      Character*180 Get_Ln
      External Get_Ln
      integer iLU,ig,natoms
      integer iat
      Character Line*180
      character lbl*20
      character cName*6

      if ((ig.lt.1).or.(ig.gt.2))
     &    Call SysAbendMsg("ReadGeo",
     &       "Wrong ig ","Shoot the programmer")

c read natoms
      Line=Get_ln(iLU)
      Call Put_ln(Line)
      Call Get_I(1,natoms,1)

      If (natoms.gt.500) Call SysAbendMsg("ReadGeo",
     & "Too many atoms in geom","")

      if (debug) Write(6,*) 'In READGEO : Nat=',natoms

      nat(ig)=natoms
c
c allocate memory for the geometry
c
      if (ig.lt.10) then
        Write(cName,'(a5,i1)') "GEOM0",ig
        title(ig)=cName
      else if (ig.lt.100) then
        Write(cName,'(a4,i2)') "GEOM",ig
        title(ig)=cName
      end if
      Call GetMem(cName,"Allo","Real",ipgeo(ig),3*nat(ig))
c
c read title
c
      Line=Get_ln(iLU)
      title(ig)=Trim(line)
c
c read label and coords
c
      iat=0
666     continue
        Line=Get_ln(iLU)
        Call Put_ln(Line)
        Call Get_S(1,lbl,1)
        if (lbl.eq.'END ') then
          Goto 999
        else
          iat=iat+1
          if (iat.gt.nat(ig)) Call SysAbendMsg("ReadGeo",
     &        "More atoms read than declared", "")

          geoLbl(iat,ig)=lbl
          Call Get_F(2,Work(ipgeo(ig)+(iat-1)*3),3)
        end if
      goto 666

      Call SysAbendMsg("ReadGeo","you should never get here",
     &     "shoot the programmer")

999   continue
      call PrintGeom(-1,nat(ig),title(ig),Work(ipgeo(ig)),geolbl(1,ig))
      return
      end
