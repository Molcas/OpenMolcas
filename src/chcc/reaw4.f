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
        subroutine ReaW4 (W,Wa,aSGrp,bSGrp,cSGrp,dSGrp,LunAux)
c
c        this routine do:
c        Add W(a",b",c",d") <<- (a",b"|c",d") for any a",b",c",d"
c        via reading records W4 (p">=q")>=(r">=s")
c
c        Wa is an auxiliary file
c
c        Structure of records:
c        1) only records for pSGrp>=qSGrp & rSGrp>=sSGrp & pqSGrp>=rsSGrp
c           exists
c        2) in individual records p">=q",r">=s" are stored
c           however p">=q")>=(r">=a") reduction is not yet used
c
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aSGrp,bSGrp,cSGrp,dSGrp,LunAux
        real*8 W(1)
        real*8 Wa(1)
c
c        help variables
        integer abPerm,cdPerm,abcdPerm
        integer abSGrp,cdSGrp
        integer abLen,cdLen,abcdLen
        integer dima,dimb,dimc,dimd
        integer pSGrp,qSGrp,rSGrp,sSGrp,dimp,dimq,dimr,dims
        character*10 LunName
        integer i
c
c
c        -------- part I - define basic parameters
c
c1        Def dima,dimb,dimc,dimc
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpbe(bSGrp)
        dimc=DimSGrpa(cSGrp)
        dimd=DimSGrpbe(dSGrp)
c
c
c        In steps 2.1 - 2.3 also dimp-dims and pSGrp-sSGrp
c2.1        Def abPerm, abSGrp
        if (aSGrp.ge.bSGrp) then
          abPerm=0
          abSGrp=(aSGrp*(aSGrp-1)/2)+bSGrp
          pSGrp=aSGrp
          qSGrp=bSGrp
          dimp=dima
          dimq=dimb
        else
          abPerm=1
          abSGrp=(bSGrp*(bSGrp-1)/2)+aSGrp
          pSGrp=bSGrp
          qSGrp=aSGrp
          dimp=dimb
          dimq=dima
        end if
c
c2.2        Def cdPerm, cdSGrp
        if (cSGrp.ge.dSGrp) then
          cdPerm=0
          cdSGrp=(cSGrp*(cSGrp-1)/2)+dSGrp
          rSGrp=cSGrp
          sSGrp=dSGrp
          dimr=dimc
          dims=dimd
        else
          cdPerm=1
          cdSGrp=(dSGrp*(dSGrp-1)/2)+cSGrp
          rSGrp=dSGrp
          sSGrp=cSGrp
          dimr=dimd
          dims=dimc
        end if
c
c2.3        Def abcdPerm
        if (abSGrp.ge.cdSGrp) then
          abcdPerm=0
        else
          abcdPerm=1
          i=pSGrp
          pSGrp=rSGrp
          rSGrp=i
          i=qSGrp
          qSGrp=sSGrp
          sSGrp=i
          i=dimp
          dimp=dimr
          dimr=i
          i=dimq
          dimq=dims
          dims=i
        end if
c
c
c3.1        Def abLen
        if (aSGrp.eq.bSGrp) then
          abLen=dima*(dima+1)/2
        else
          abLen=dima*dimb
        end if
c
c3.2        Def cdLen
        if (cSGrp.eq.dSGrp) then
          cdLen=dimc*(dimc+1)/2
        else
          cdLen=dimc*dimd
        end if
c
c3.3        Def abcdLen
        abcdLen=abLen*cdLen
c
c
c        -------- part II - read proper integrals (pq|rs) from disc
c
c4        Def proper LunName
        call MkNameV4 (pSGrp,qSGrp,rSGrp,sSGrp,'W4',LunName)
c
c
c5        Read Wa(pqrs) - velka udalost
C       open (unit=LunAux,file=LunName,form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,LunName)
        call rea_chcc (LunAux,abcdLen,Wa(1))
        close (LunAux)
c
c
c        -------- part III - Upgrade W(a",b",c",d") from Wx(pqrs)
c        symbolics used: [xy] means xy or yx
c
        if (abcdPerm.eq.0) then
c          case ([ab]|[cd])
c
          if ((abPerm.eq.0).and.(cdPerm.eq.0)) then
c            subcase (ab|cd)
            call DefW4abcd (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else if ((abPerm.eq.0).and.(cdPerm.eq.1)) then
c            subcase (ab|dc)
            call DefW4abdc (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else if ((abPerm.eq.1).and.(cdPerm.eq.0)) then
c            subcase (ba|cd)
            call DefW4bacd (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else
c            subcase (ba|cd)
            call DefW4badc (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          end if
c
        else
c        case ([cd]|[ab])
c
          if ((abPerm.eq.0).and.(cdPerm.eq.0)) then
c            subcase (cd|ab)
            call DefW4cdab (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else if ((abPerm.eq.0).and.(cdPerm.eq.1)) then
c            subcase (cd|ab)
            call DefW4dcab (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else if ((abPerm.eq.1).and.(cdPerm.eq.0)) then
c            subcase (cd|ba)
            call DefW4cdba (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          else
c            subcase (dc|ba)
            call DefW4dcba (W,Wa,dima,dimb,dimc,dimd,abLen,cdLen,
     c                      aSGrp,bSGrp,cSGrp,dSGrp)
c
          end if
c
        end if
c
c
        return
        end
