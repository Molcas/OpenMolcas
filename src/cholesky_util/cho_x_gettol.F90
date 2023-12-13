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
!> If the integrals have not been Cholesky decomposed
!> (or represented with DF), this function simply
!> returns \p iTolDef.
!>
!> @param[in] iTolDef Default tolerance
!>
!> @return Tolerance integer for use with ::Add_Info
!***********************************************************************

function Cho_X_GetTol(iTolDef)

use Cholesky, only: ChoIniCheck, ThrCom
use Constants, only: Ten
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: Cho_X_GetTol
integer(kind=iwp), intent(in) :: iTolDef
integer(kind=iwp) :: ChoIsIni
real(kind=wp) :: d, ThrAbs
logical(kind=iwp) :: DidCholesky

call DecideOnCholesky(DidCholesky)
if (DidCholesky) then
  call Get_iScalar('ChoIni',ChoIsIni)
  if (ChoIsIni /= ChoIniCheck) call Get_dScalar('Cholesky Threshold',ThrCom) ! not initialized
  ThrAbs = abs(ThrCom)
  d = -log(ThrAbs)/log(Ten)
  Cho_X_GetTol = nint(d)
else
  Cho_X_GetTol = iTolDef
end if

end function Cho_X_GetTol
