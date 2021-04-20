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
        subroutine Ext_W4 (V2,M1,
     c                     nc,dima,dimb,dimab,dimapp,dimbpp,dimabpp,
     c                     addapp,addbpp,aGrp,bGrp,aSGrp,bSGrp)
c
c        this routine is a control routine to:
c        Extract M1(m,a"b") <- V2(m,a'b')
c        for all combinations of aGrp,bGrp,aSGrp,bSGrp
c
        implicit none
        integer nc,dima,dimb,dimab,dimapp,dimbpp,dimabpp
        integer addapp,addbpp,aGrp,bGrp,aSGrp,bSGrp
        real*8 M1(1)
        real*8 V2(1)
c
        if (aGrp.eq.bGrp) then
c          case V2(m,a'b')
          if (aSGrp.eq.bSGrp) then
c            subcase M1(m,a"b") <- V2(m,a'b')
            call Ext_W4hlp1 (V2,M1,nc,dima,dimab,dimapp,dimabpp,addapp)
          else
c            subcase M1(m,a",b") <- V2(m,a'b')
            call Ext_W4hlp2 (V2,M1,
     c                       nc,dima,dimab,dimapp,dimbpp,addapp,addbpp)
          end if
        else
c          case M1(m,a",b") <- V2(m,a',b')
          call Ext_W4hlp3 (V2,M1,
     c                     nc,dima,dimb,dimapp,dimbpp,addapp,addbpp)
        end if
c
        return
        end
