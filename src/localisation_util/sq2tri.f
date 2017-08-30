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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  Sq2Tri
*
*> @brief
*>   Convert from square to lower triangular storage
*> @author Thomas Bondo Pedersen
*>
*> @param[in]  Sq  Square storage array
*> @param[out] Tri Lower triangular storage array
*> @param[in]  n   Dimension
*>
*> @details
*> Perform the extraction
*>
*> \code
*> Tri(i*(i-1)/2+j) = Sq(i,j)
*> \endcode
*>
*> where \c i &ge; \c j.
************************************************************************
      SubRoutine Sq2Tri(Sq,Tri,n)
      Implicit None
      Integer n
      Real*8  Sq(n,n), Tri(n*(n+1)/2)

      Integer i,j, iTri

      iTri(i,j)=i*(i-1)/2 + j

      Do j = 1,n
         Do i = j,n
            Tri(iTri(i,j)) = Sq(i,j)
         End Do
      End Do

      End
