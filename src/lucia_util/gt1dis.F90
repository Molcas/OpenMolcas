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

!#define _DEBUGPRINT_
subroutine GT1DIS(H1DIA,IREOTS,IPNT,H,ISMFTO,IBSO,NACOB)
! diagonal of one electron integrals over active orbitals

use Index_Functions, only: nTri_Elem
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: H1DIA(*)
integer(kind=iwp), intent(in) :: NACOB, IREOTS(NACOB), IPNT(*), ISMFTO(NACOB), IBSO(*)
real(kind=wp), intent(in) :: H(*)
integer(kind=iwp) :: IIOB, IOB, IOBREL, ISM

do IIOB=1,NACOB
  IOB = IREOTS(IIOB)
  ISM = ISMFTO(IIOB)
  IOBREL = IOB-IBSO(ISM)+1
  !write(u6,*) ' IIOB IOB ISM IOBREL'
  !write(u6,*) IIOB,IOB,ISM,IOBREL
  H1DIA(IIOB) = H(IPNT(ISM)-1+nTri_Elem(IOBREL))
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' Diagonal one electron integrals'
call WRTMAT(H1DIA,1,NACOB,1,NACOB)
#endif

end subroutine GT1DIS
