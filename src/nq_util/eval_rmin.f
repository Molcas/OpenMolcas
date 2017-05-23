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
      Function Eval_RMin(Alpha,m,R_H)
      Implicit Real*8 (a-h,o-z)
#include "itmax.fh"
#include "real.fh"
      Real*8 Eval_RMin, ln_x
*     Write (6,*) 'Alpha,m,R_H=',
*    &            Alpha,m,R_H
*                                                                      *
************************************************************************
*                                                                      *
*---- Compute r_1 as a function of m and R_H
*
*     Eq(25) R. Lindh, P.-A. Malmqvist, L. Gagliardi,
*     TCA, 106:178-187 (2001)
*
      D_m=-4.0D0
      If (m.eq. 4) D_m=-2.3D0
      If (m.eq. 2) D_m=-1.0D0
      If (m.eq. 0) D_m= 1.9D0
      If (m.eq.-2) D_m= 9.1D0
*
*     x = alpha * (r_1)**2
*
      ln_x=(Two/(DBLE(m)+Three))*( D_m- LOG(One/R_H) )
*
      R_Min=Sqrt(Exp(ln_x)/Alpha)
*
      Eval_RMin = R_Min
*     Write (*,*) 'Eval_RMin=',Eval_RMin
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
