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
* Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
*               1996-2006, David L. Cooper                             *
************************************************************************
      subroutine zz_cvb(act1,zz1,fx,fxbest1,exp1,ip1)
      implicit real*8 (a-h,o-z)
#include "opt_cvb.fh"
#include "locopt1_cvb.fh"
#include "locopt2_cvb.fh"
#include "trst_cvb.fh"
#include "tune_cvb.fh"

#include "formats_cvb.fh"
      save zero,one,thous
      data zero/0.d0/,one/1.d0/,thous/1d3/

      act1=fx-fxbest1
      if(fxbest1.eq.-thous)act1=one
      if((abs(act1).lt.dfxtol.and.abs(exp1).lt.dfxtol)
     >  .or.act1.eq.one)then
        zz1=one
      elseif(exp1.eq.zero)then
        zz1=one
      elseif(abs(exp1).lt.dfxtol)then
        zz1=sign(one,act1)*sign(one,exp1)
      else
        zz1=act1/exp1
      endif
      if(ip1.ge.2)then
        if(act1.ne.one)write(6,formAD)
     >    ' Actual and expected changes :',act1,exp1
        write(6,formAD)' Ratio act/exp    : ',zz1
      endif
      return
      end
