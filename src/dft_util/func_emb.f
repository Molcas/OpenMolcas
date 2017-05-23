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
* Copyright (C) 2010, Francesco Aquilante                              *
************************************************************************
      Real*8 Function Fexp(rho,drho)
************************************************************************
************************************************************************
***                                                                  ***
*** Switching function of the NDSD potential:                        ***
***     J.-M. Garcia Lastra, J. W. Kaminski, T. A. Wesolowski,       ***
***                               J. Chem. Phys.  129 (2008) 074107. ***
***                                                                  ***
*** Author: F. Aquilante, Geneva July 2010                           ***
***                                                                  ***
************************************************************************
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Real*8 lambda, rho, drho(3)
#include "real.fh"
      Parameter ( One3 = One/Three )
      Parameter ( sBmin = 0.3d0 )
      Parameter ( sBmax = 0.9d0 )
      Parameter ( rBmin = 0.7d0 )
      Parameter ( lambda= 5.0d2 )


      rhoinv= One/rho
      rhoinv13= rhoinv**One3
      fact= Two*(Three*Pi**2)**One3
      factinv= One/fact
      xnorm = drho(1)**2 + drho(2)**2 + drho(3)**2
      xnorm_= sqrt(xnorm)
      sB = factinv*rhoinv*rhoinv13*xnorm_

      es_sBmin = exp(lambda*(sBmin-sB))
      es_sBmax = exp(lambda*(sBmax-sB))
      er_rBmin = exp(lambda*(rBmin-rho))

      eis_sBmin= One/(es_sBmin+One)
      eis_sBmax= One/(es_sBmax+One)
      eir_rBmin= One/(er_rBmin+One)

      Fexp = eis_sBmin*(One-eis_sBmax)*eir_rBmin

      End
*                                                                      *
************************************************************************
*                                                                      *
      Real*8 Function Vt_lim(rho,drho,ddrho)

      Implicit Real*8 (a-h,o-z)
      Real*8 rho, drho(3), ddrho
#include "real.fh"
      Parameter ( One8 = One/Eight )
      Parameter ( One4 = One/Four  )


      rhoinv= One/rho
      rhoinv2= rhoinv**Two
      xnorm = drho(1)**2 + drho(2)**2 + drho(3)**2

      Vt_lim = One8*xnorm*rhoinv2 - One4*ddrho*rhoinv

      End
*                                                                      *
************************************************************************
*                                                                      *
