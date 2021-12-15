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
* Copyright (C) 2004,2005, Alexander Wolf                              *
*               2004,2005, Markus Reiher                               *
************************************************************************
        subroutine finalize (length,operator,stimes,ttimes,t)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 13.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer length,stimes(maxsnumber),ttimes(maxsnumber),k,pos,idum
      character*(maxlength) operator
      character*(4) t(maxsnumber)
      integer dkh_char2int
1000  continue
      pos=0
      pos=index(operator(1:length),'[PVP]')
      if (pos.gt.0) then
        operator(pos:pos)='Y'
        do 100 k=pos+1,length-4
          operator(k:k)=operator(k+4:k+4)
 100    continue
        operator(length-3:length)='    '
        length=length-4
        goto 1000
      endif
1100  continue
      pos=0
      pos=index(operator(1:length),'[PXP]')
      if (pos.gt.0) then
        operator(pos:pos)='K'
        do 110 k=pos+1,length-4
          operator(k:k)=operator(k+4:k+4)
 110    continue
        operator(length-3:length)='    '
        length=length-4
        goto 1100
      endif
1500  continue
      pos=0
      pos=index(operator(1:length),'PVP')
      if (pos.gt.0) then
        operator(pos:pos)='D'
        do 150 k=pos+1,length-2
          operator(k:k)=operator(k+2:k+2)
 150    continue
        operator(length-1:length)='  '
        length=length-2
        goto 1500
      endif
1510  continue
      pos=0
      pos=index(operator(1:length),'PXP')
      if (pos.gt.0) then
        operator(pos:pos)='J'
        do 151 k=pos+1,length-2
          operator(k:k)=operator(k+2:k+2)
 151    continue
        operator(length-1:length)='  '
        length=length-2
        goto 1510
      endif
 500  continue
      pos=0
      pos=index(operator(1:length),'PS')
      if (pos.gt.0) then
        idum=dkh_char2int(3,operator(pos+2:pos+4))
        if (operator(pos+5:pos+5).ne.'P') then
          write (stdout,2003)
2003      format (2X,'ERROR in SR "finalize" while substituting ',
     *            'PSxxxP structure.',//2X,'STOP.',/2X)
          CALL Abend
        endif
        ttimes(idum)=ttimes(idum)+1
        stimes(idum)=stimes(idum)-1
        operator(pos:pos+3)=t(idum)
        do 10 k=pos+4,length-2
          operator(k:k)=operator(k+2:k+2)
  10    continue
        operator(length-1:length)='  '
        length=length-2
        goto 500
      endif
2000  continue
      pos=0
      pos=index(operator(1:length),'[V]')
      if (pos.gt.0) then
        operator(pos:pos)='N'
        do 200 k=pos+1,length-2
          operator(k:k)=operator(k+2:k+2)
 200    continue
        operator(length-1:length)='  '
        length=length-2
        goto 2000
      endif
2010  continue
      pos=0
      pos=index(operator(1:length),'[X]')
      if (pos.gt.0) then
        operator(pos:pos)='I'
        do 210 k=pos+1,length-2
          operator(k:k)=operator(k+2:k+2)
 210    continue
        operator(length-1:length)='  '
        length=length-2
        goto 2010
      endif
2500  continue
      pos=0
      pos=index(operator(1:length),'PE01P')
      if (pos.gt.0) then
        operator(pos:pos)='G'
        do 250 k=pos+1,length-4
          operator(k:k)=operator(k+4:k+4)
 250    continue
        operator(length-3:length)='    '
        length=length-4
        goto 2500
      endif
2520  continue
      pos=0
      pos=index(operator(1:length),'PCE0P')
      if (pos.gt.0) then
        operator(pos:pos)='M'
        do 254 k=pos+1,length-4
          operator(k:k)=operator(k+4:k+4)
 254    continue
        operator(length-3:length)='    '
        length=length-4
        goto 2520
      endif
3000  continue
      pos=0
      pos=index(operator(1:length),'E01')
      if (pos.gt.0) then
        operator(pos:pos)='F'
        do 300 k=pos+1,length-2
          operator(k:k)=operator(k+2:k+2)
 300    continue
        operator(length-1:length)='  '
        length=length-2
        goto 3000
      endif
3040  continue
      pos=0
      pos=index(operator(1:length),'CE0')
      if (pos.gt.0) then
        operator(pos:pos)='L'
        do 330 k=pos+1,length-2
          operator(k:k)=operator(k+2:k+2)
 330    continue
        operator(length-1:length)='  '
        length=length-2
        goto 3040
      endif
3500  continue
      pos=0
      pos=index(operator(1:length),'PP')
      if (pos.gt.0) then
        operator(pos:pos)='Z'
        do 350 k=pos+1,length-1
          operator(k:k)=operator(k+1:k+1)
 350    continue
        operator(length:length)=' '
        length=length-1
        goto 3500
      endif
c
      return
      end
