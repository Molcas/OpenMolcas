!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************
!  Cho_X_GetTol
!
!> @brief
!>   Set tolerance integer for use with ::Add_Info (for verification).
!> @author Thomas Bondo Pedersen
!> @modified_by Thomas Bondo Pedersen, October 2010: LDF support.
!>
!> @details
!> The tolerance in verification might depend on the
!> threshold of the Cholesky decomposition. The integer
!> used to specify tolerance in ::Add_Info is computed
!> here according to the formula:
!>
!> \f[ \text{iTol} = -\log(\text{Thr}) \f]
!>
!> where \p Thr is the threshold and \f$ \log \f$ is the logarithm
!> (base 10).
!>
!> If LDF (local DF) is used, \p Thr is the LDF target
!> accuracy.
!> If the integrals have not been Cholesky decomposed
!> (or represented with DF or LDF), this function simply
!> returns \p iTolDef.
!>
!> @param[in] iTolDef Default tolerance
!>
!> @return Tolerance integer for use with ::Add_Info
!***********************************************************************

integer function Cho_X_GetTol(iTolDef)

use ChoIni

implicit none
integer iTolDef
#include "cholesky.fh"
real*8, external :: Get_LDFAccuracy
logical DidCholesky, DidLDF
real*8 ThrAbs, d
integer ChoIsIni

call DecideOnCholesky(DidCholesky)
if (DidCholesky) then
  call DecideOnLocalDF(DidLDF)
  if (DidLDF) then
    ThrAbs = abs(Get_LDFAccuracy())
  else
    call Get_iScalar('ChoIni',ChoIsIni)
    if (ChoIsIni /= ChoIniCheck) call Get_dScalar('Cholesky Threshold',ThrCom) ! not initialized
    ThrAbs = abs(ThrCom)
  end if
  d = -log(ThrAbs)/log(1.0d1)
  Cho_X_GetTol = nint(d)
else
  Cho_X_GetTol = iTolDef
end if

end function Cho_X_GetTol
