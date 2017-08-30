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
* Copyright (C) Yannick Carissan                                       *
************************************************************************
*  Cross
*
*> @brief
*>   Performs the cross product \p A &times; \p B and puts the result in \p C
*> @author Y. Carissan
*>
*> @details
*> Performs the cross product \p A &times; \p B and puts the result in \p C.
*>
*> @param[in]  A First (left) input vector
*> @param[in]  B Second (right) input vector
*> @param[out] C Output vector
************************************************************************
      Subroutine Cross(A,B,C)
      Implicit Real*8(a-h,o-z)
#include "WrkSpc.fh"
#include "debug.fh"

      Real*8 A(3),B(3),C(3)

      C(1) = A(2)*B(3) - B(2)*A(3)
      C(2) = A(3)*B(1) - B(3)*A(1)
      C(3) = A(1)*B(2) - B(1)*A(2)

      Return
      End
