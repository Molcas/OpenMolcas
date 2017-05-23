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
      subroutine main_bdata_cvb()
      implicit real*8 (a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"

#include "seth_cvb.fh"
#include "spinb_cvb.fh"

      md2h(1,1)=1
      md2h(2,1)=2
      md2h(3,1)=3
      md2h(4,1)=4
      md2h(5,1)=5
      md2h(6,1)=6
      md2h(7,1)=7
      md2h(8,1)=8

      md2h(1,2)=2
      md2h(2,2)=1
      md2h(3,2)=4
      md2h(4,2)=3
      md2h(5,2)=6
      md2h(6,2)=5
      md2h(7,2)=8
      md2h(8,2)=7

      md2h(1,3)=3
      md2h(2,3)=4
      md2h(3,3)=1
      md2h(4,3)=2
      md2h(5,3)=7
      md2h(6,3)=8
      md2h(7,3)=5
      md2h(8,3)=6

      md2h(1,4)=4
      md2h(2,4)=3
      md2h(3,4)=2
      md2h(4,4)=1
      md2h(5,4)=8
      md2h(6,4)=7
      md2h(7,4)=6
      md2h(8,4)=5

      md2h(1,5)=5
      md2h(2,5)=6
      md2h(3,5)=7
      md2h(4,5)=8
      md2h(5,5)=1
      md2h(6,5)=2
      md2h(7,5)=3
      md2h(8,5)=4

      md2h(1,6)=6
      md2h(2,6)=5
      md2h(3,6)=8
      md2h(4,6)=7
      md2h(5,6)=2
      md2h(6,6)=1
      md2h(7,6)=4
      md2h(8,6)=3

      md2h(1,7)=7
      md2h(2,7)=8
      md2h(3,7)=5
      md2h(4,7)=6
      md2h(5,7)=3
      md2h(6,7)=4
      md2h(7,7)=1
      md2h(8,7)=2

      md2h(1,8)=8
      md2h(2,8)=7
      md2h(3,8)=6
      md2h(4,8)=5
      md2h(5,8)=4
      md2h(6,8)=3
      md2h(7,8)=2
      md2h(8,8)=1

      zero=0d0
      one=1.d0
      two=2d0
      onem=-1d0
      r3by4=.75d0
      p8=.8d0
      sqp5=.70710678118654752440d0
      sq2=1.41421356237309504880d0
      mxopt=500
      mxcnt=5000
      opposite=.false.
      corenrg=0d0
      spinb(1)='Kotani      '
      spinb(2)='Serber      '
      spinb(3)='Rumer       '
      spinb(4)='Rumer (LT)  '
      spinb(5)='Projected   '
      spinb(6)='Determinants'
      spinb(7)='Determinants'
      spinbkw(1)='KOTANI  '
      spinbkw(2)='SERBER  '
      spinbkw(3)='RUMER   '
      spinbkw(4)='LTRUMER '
      spinbkw(5)='PROJECT '
      spinbkw(6)='DET     '
      spinbkw(7)='DETERM  '
      return
      end
