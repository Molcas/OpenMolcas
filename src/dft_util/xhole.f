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
*  XHole
*
*> @brief
*>   Prepare the kernel for the exchange-hole dipole
*> @author A. Ohrn
*>
*> @details
*> At this moment only for RHF. Collect number of occupied orbitals. Then
*> construct the kernel. Since orbital pair dipoles are triangular, the
*> matrix is unfolded by a factor two for non-diagonal elements. For
*> small values of the density, the exchange hole dipole grid value is
*> set to zero; this is okey since it will only be roughly of size 100,
*> which becomes as good as zero when later multiplied by the density.
*> See \cite Bec2005-JCP-122-154104 for details.
*>
*> @param[in]  nRho      Should equal one, since no derivatives
*> @param[in]  mGrid     Number of grid-points
*> @param[in]  Rho       The density in the grid-points
*> @param[in]  Grid      The coordinates of the grid-points
*> @param[in]  mAO       Should equal one, since no derivatives
*> @param[in]  nMOs      Total number of molecular orbitals
*> @param[in]  TabMO     The orbitals in the grid-points
*> @param[in]  ndF_dRho  Should equal one
*> @param[in]  nD        Should equal one
*> @param[out] dF_dRho   The value of the kernel. This is used as input later to compute the actual matrix elements,
*>                       but we do so by fooling the program to believe that it is the derivate of the functional
*>                       with respect to density
*> @param[in]  Weights   Quadrature weights
*> @param[in]  ip_OrbDip Dipole moments computed by LoProp routine
*> @param[out] Func      The expectation value of the kernel
************************************************************************
#ifdef _NOT_USED_OR_TESTED_
      Subroutine XHole(nRho,mGrid,Rho,Grid,mAO,nMOs,TabMO,ndF_dRho,nD,
     &                 dF_dRho,Weights,ip_OrbDip,Func)
      Implicit Real*8 (a-h,o-z)

#include "real.fh"
#include "nq_index.fh"
#include "WrkSpc.fh"

      Real*8 Rho(nRho,mGrid),Grid(3,mGrid),Weights(mGrid)
      Real*8 TabMO(mAO,mGrid,nMOs)
      Real*8 dF_dRho(ndF_dRho,mGrid)
      Real*8 dx(3)
      Integer nOcc(1),ip_OrbDip(3)

*
*-- Some crazy consistency tests.
*
      If(nRho.ne.1.or.nD.ne.1.or.mAO.ne.1) then
        Call WarningMessage(2,
     &            'How did you manage this! nRho, nD or mAO is not '
     &          //'equal to one in Do_XHole!')
        Call Abend()
      Endif

*
*-- Some numbers.
*
      Call Get_iArray('nIsh',nOcc,1)
      nOccMO=nOcc(1)

*
*-- Now evaluate dx**2. Put it in dF_dRho for later processing.
*
      Do iGrid=1,mGrid
        If(Rho(1,iGrid).gt.1d-14) then
          Do k=1,3
            kaunter=0
            SumOrb=Zero
            Do iMO=1,nOccMO
              Do jMO=1,iMO
                kaunter=kaunter+1
                dUnfoldToSq=One
                If(iMO.ne.jMO) dUnfoldToSq=Two
                dMOPr=dUnfoldToSq*TabMO(1,iGrid,iMO)*TabMO(1,iGrid,jMO)
                SumOrb=SumOrb+Two*Work(ip_OrbDip(k)+kaunter-1)*dMOPr
              Enddo
            Enddo
            dx(k)=(One/(Rho(1,iGrid)))*SumOrb-Grid(k,iGrid)
          Enddo
        Else
          Do k=1,3
            dx(k)=Zero
          Enddo
        Endif
        dSquare=dx(1)**2+dx(2)**2+dx(3)**2
        dF_dRho(ipR,iGrid)=dSquare
        Func=Func+Weights(iGrid)*dSquare*Rho(1,iGrid)
      Enddo

      Return
      End
#else
      Subroutine XHole()
      End
#endif
