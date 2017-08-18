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
* Copyright (C) Anders Ohrn                                            *
************************************************************************
*  Hel
*
*> @brief
*>   Couple the electrostatic part of the solvent with the QM-region.
*>   Only include the static part, no polarization at this moment
*> @author A. Ohrn
*>
*> @details
*> (2) The electrostatics.
*>
*> @param[in]  Eint   The static field from the solvent on the QM molecule centers
*> @param[in]  itri   Number of elements in triangular \f$ H \f$-matrix
*> @param[in]  ici    Number of MME-centers
*> @param[in]  ql     MME-charges, obtained from the MME
*> @param[in]  dil    MME-dipoles
*> @param[in]  qqxxyy MME-quadrupoles.
*> @param[out] vmat   The electrostatic part of the solute-solvent interaction matrix
*> @param[in]  iprint Print parameter
************************************************************************
      Subroutine Hel(Eint,itri,ici,ql,dil,qqxxyy,vmat,iprint)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"

      Dimension Ql(MxOT,MxQCen),Dil(MxOT,3,MxQCen)
     &,QQxxyy(MxOT,6,MxQCen),Eint(MxQCen,10),Vmat(MxOT)


*Zeros
      Do 9, i=1,itri
        Vmat(i)=0.0d0
9     Continue

*The electrostatic perturbation: <psi_i|V_el|psi_j>
      Do 10, i=1,itri
        Do 11, k=1,ici
          Vmat(i)=Vmat(i)+Eint(k,1)*Ql(i,k)
          Do 12, j=1,3
            Vmat(i)=Vmat(i)+Eint(k,j+1)*Dil(i,j,k)
12        Continue
          Vmat(i)=Vmat(i)+Eint(k,5)*QQxxyy(i,1,k)
          Vmat(i)=Vmat(i)+Eint(k,7)*QQxxyy(i,3,k)
          Vmat(i)=Vmat(i)+Eint(k,10)*QQxxyy(i,6,k)
          Vmat(i)=Vmat(i)+Eint(k,6)*QQxxyy(i,2,k)*2.0d0
          Vmat(i)=Vmat(i)+Eint(k,8)*QQxxyy(i,4,k)*2.0d0
          Vmat(i)=Vmat(i)+Eint(k,9)*QQxxyy(i,5,k)*2.0d0
11      Continue
10    Continue

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(iprint)
      End
