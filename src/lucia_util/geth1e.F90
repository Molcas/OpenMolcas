!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

function GETH1E(IORB,ITP,ISM,JORB,JTP,JSM)
! One-electron integral for active
! orbitals (IORB,ITP,ISM),(JORB,JTP,JSM)
!
! The orbital symmetries are used to obtain the
! total symmetry of the operator

use Symmetry_Info, only: Mul
use lucia_data, only: IBSO, IH1FORM, INT1, IOBPTS, IREOTS, MXPNGAS, NACOBS, PGINT1, PGINT1A, PINT1
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: GETH1E
integer(kind=iwp), intent(in) :: IORB, ITP, ISM, JORB, JTP, JSM
integer(kind=iwp) :: IJSM
real(kind=wp), external :: GTH1ES

IJSM = Mul(ISM,JSM)
GETH1E = Zero
if (IH1FORM == 1) then
  ! Normal integrals, lower triangular packed
  if (IJSM == 1) then
    !write(u6,*) ' GETH1E, old route'
    GETH1E = GTH1ES(IREOTS,PINT1,INT1,IBSO,MXPNGAS,IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
  else
    GETH1E = GTH1ES(IREOTS,PGINT1(IJSM)%A,INT1,IBSO,MXPNGAS,IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,1)
  end if
else
  ! Integrals are in full blocked form
  GETH1E = GTH1ES(IREOTS,PGINT1A(IJSM)%A,INT1,IBSO,MXPNGAS,IOBPTS,NACOBS,IORB,ITP,ISM,JORB,JTP,JSM,0)
end if

end function GETH1E
