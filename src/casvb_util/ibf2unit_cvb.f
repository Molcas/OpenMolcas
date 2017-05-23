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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine ibf2unit_cvb(ibf,lu,newfile)
      implicit real*8(a-h,o-z)
#include "io_cvb.fh"
      logical newfile,debug
      data debug/.false./

      if(ifilio(ibf).ne.0)then
        newfile=.false.
        ifil=ifilio(ibf)
      else
        newfile=.true.
        ifil=0
        do 100 i=1,mxunits
        if(iorder(i).eq.0)then
          iorder(i)=i
          ifil=i
          goto 200
        endif
100     continue
        do 300 i=1,mxunits
        if(iorder(i).eq.mxunits)then
          ifil=i
          goto 200
        endif
300     continue
        if(ifil.eq.0)then
          write(6,*)' ifil error - iorder :',iorder
          call abend_cvb()
        endif
200     ifilio(ibf)=ifil
      endif
      call touchord_cvb(ifil,iorder,mxunits)
      lu=90+ifil
      if(debug)write(6,*)' buffer:',ibf,' now associated with unit:',lu
      return
      end
