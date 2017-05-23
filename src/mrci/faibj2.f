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
      subroutine faibj2(IFTA,IFTB,ICOUP1,ICOUP,
     & INDA,INDB,MYSYM,INTSYM,
     & NYSYM,NSIJ,MYL,NYL,FACS,IPOA,IPOB,
     & INMY,INNY,INDX,iTYP)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
      DIMENSION INTSYM(*),INDX(*)
      DIMENSION IPOA(9),IPOB(9)

      IFTA=0
      IFTB=0
c      GO TO (109,110,111,112,113),ITYP
      IF(ITYP.eq.1) then
      INDA=IRC(2)+ICOUP1
      INDB=IRC(2)+ICOUP
      IFTA=1
      IFTB=1
      endif
      if(ITYP.eq.2) then
      INDA=IRC(3)+ICOUP1
      INDB=IRC(3)+ICOUP
      endif
      if(ITYP.eq.3) then
      INDA=IRC(2)+ICOUP1
      INDB=IRC(3)+ICOUP
      IFTA=1
      endif
      if(ITYP.eq.4) then
      INDA=IRC(3)+ICOUP1
      INDB=IRC(2)+ICOUP
      IFTB=1
      endif
      if(ITYP.eq.5) then
      INDA=IRC(1)+ICOUP1
      INDB=IRC(1)+ICOUP
      endif
cvv : unroll inline function to make GCC compiler works proper..
      MYSYM=JSUNP(INTSYM,INDA)
      NYSYM=MUL(MYSYM,NSIJ)
      MYL=MUL(MYSYM,LSYM)
      NYL=MUL(NYSYM,LSYM)
      FACS=1.0D00
      CALL IPO(IPOA,NVIR,MUL,NSYM,MYL,IFTA)
      CALL IPO(IPOB,NVIR,MUL,NSYM,NYL,IFTB)
      INMY=INDX(INDA)+1
      INNY=INDX(INDB)+1
      return
      end
