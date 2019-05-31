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
* Copyright (C) 1994,2004,2014,2017,2019, Roland Lindh                 *
*               2014,2018, Ignacio Fdez. Galvan                        *
************************************************************************
      Real*8 Function Restriction_Step(q,dq,nInter)
************************************************************************
*                                                                      *
*     Object: External routine to evaluate a general constraint,       *
*             to be used in a constraint optimization. In this case    *
*             the constraint is a step size constraint.                *
*                                                                      *
*     Author: Roland Lindh, Dept. of Chemistry - BMC                   *
*             University of Uppsala, Sweden                            *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer nInter
      Real*8 q(nInter), dq(nInter)
*
      Restriction_Step=DDot_(nInter,dq,1,dq,1)
*
      If (.False.) Call Unused_real_array(q)
      Return
      End
