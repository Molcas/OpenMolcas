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
      real*8 function  gammat(x)
************************************************************************
*                                                                      *
* Object: to compute the angular contribution to the multipole integral*
*         between continuum basis functions within an R-matrix sphere  *
*         (theta integration)                                          *
*                                                                      *
************************************************************************
      Implicit real*8 (a-h,o-z)
#include "nrmf.fh"
#include "gam.fh"
*
      lsint=m_gam
      lcost=n_gam
      k=(-1)**lcost
      if(k.eq.(-1)) then
       gammat=0.0d0
      else
       arg1=(DBLE(lcost)+1.0d0)/2.0d0
       arg2=(DBLE(lsint)+2.0d0)/2.0d0
       arg3=(DBLE(lsint)+DBLE(lcost)+3.0d0)/2.0d0
       gammat=dgamma_molcas(arg1)*dgamma_molcas(arg2)/
     >                            dgamma_molcas(arg3)
      Endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(x)
      End
