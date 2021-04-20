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
        subroutine ReaW3 (Ww,Wx,aSGrp,beSGrp,bSGrp,LunAux)
c
c        this routine do:
c        define Ww(a",be",b",i) <- (a",be"|b",i)
c
c        integrals (a",be"|b",i) are stored in files V3xxyyzz for xx>=yy
c
c        Storing of integrals in V3files
c        do i=1,no
c          record of V3(a"be",b",_i): a">=be",b" for given i
c        end do
c
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer aSGrp,beSGrp,bSGrp,LunAux
        real*8 Ww(1)
        real*8 Wx(1)
c
c        help variables
c
        integer dima,dimb,dimbe
        character*8 LunName

c0        def dimensions
c
        dima=DimSGrpa(aSGrp)
        dimb=DimSGrpa(bSGrp)
        dimbe=DimSGrpbe(beSGrp)
c
c
        if (aSGrp.gt.beSGrp) then
c1        case aSGrp>beSGrp, integrals (a",be"|b",i) to be read
c1.1          make Name
          call MkNameV3 (aSGrp,beSGrp,bSGrp,'W3',LunName)
c1.2          read (a",be"|b",i) from V3:(a",be"|b",i)
          call ReaW3hlp1 (Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)
c
        else if (aSGrp.eq.beSGrp) then
c2        case aSGrp=beSGrp, integrals (a">=be"|b",i)to be read and expand
c2.1          make Name
          call MkNameV3 (aSGrp,beSGrp,bSGrp,'W3',LunName)
c2.2          read (a",be"|b",i) from V3:(a">=be"|b",i)
          call ReaW3hlp2 (Ww,Wx,dima,dimb,no,LunName,LunAux)
c
        else
c3        case aSGrp<beSGrp, integrals (be",a"|b",i) to be read and mapped
c3.1          make Name
          call MkNameV3 (beSGrp,aSGrp,bSGrp,'W3',LunName)
c3.2          read (a",be"|b",i) from V3:(be",a"|b",i)
          call ReaW3hlp3 (Ww,Wx,dima,dimbe,dimb,no,LunName,LunAux)

        end if
c
        return
        end
