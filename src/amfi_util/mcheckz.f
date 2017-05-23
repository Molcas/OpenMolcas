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
      logical function mcheckz(m1,m2,m3,m4)
cbs   makes a check, if there is an interaction inbetween cartesian functions
cbs   with m-values m1-m4
      integer m1,m2,m3,m4,int12a,int12b,
     *int34a,int34b
      mcheckz=.true.
      int12a=m1+m2
      int12b=-m1+m2
      int34a=m3+m4
      int34b=-m3+m4
cbs   lots of checks
      if (iabs(int12a+int34a).eq.0) return
      if (iabs(int12a-int34a).eq.0) return
      if (iabs(int12b+int34b).eq.0) return
      if (iabs(int12b-int34b).eq.0) return
      if (iabs(int12a+int34b).eq.0) return
      if (iabs(int12a-int34b).eq.0) return
      if (iabs(int12b+int34a).eq.0) return
      if (iabs(int12b-int34a).eq.0) return
      mcheckz=.false.
      return
      end
