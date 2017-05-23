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
*
*     LHR function: convert full Hessian indices in lower triangular index
*
      Integer Function LHR(i,iAtom,j,jAtom)
      ivv1 = (iAtom-1)*3 + i
      ivv2 = (jAtom-1)*3 + j
      ivv  = Max(ivv1,ivv2)
      LHR  = ivv*(ivv-3)/2 + ivv1 + ivv2
      Return
      End
