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

subroutine GetTau(Tau,T1,aGrp,bGrp,dima,dimb,adda,addb,lunTau)
! this routine does:
! 1) read the block of T2((a,b)',ij) from T2Name(aGrp,bGrp)
! 2) Make Tau ((a,b)',ij) in T2((a,b)',ij) array
!
! I/O parameter description:
! Tau    - array for Tau((a,b),ij) (O)
! aGrp   - group of a index
! bGrp   - group of b index
! dima   - dimension group of a' index
! dimb   - dimension group of b' index
! adda   - shift of a' in full virtual space
! addb   - shift of b' in full virtual space
! lunTau - Lun of opened file, where Tau is stored

use Index_Functions, only: nTri_Elem
use chcc_global, only: no, nv, T2Name
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(_OUT_) :: Tau(*), T1(*)
integer(kind=iwp), intent(in) :: aGrp, bGrp, dima, dimb, adda, addb, lunTau
integer(kind=iwp) :: length

!1 def length

if (aGrp == bGrp) then
  ! groups of a and b are equal, reading for a'>=b'
  length = no*no*nTri_Elem(dima)
else
  ! aGrp>bGrp, reading for a',b' in given groups
  length = no*no*dima*dimb
end if

!2 read block of T2 amplitudes

call GetX(Tau,length,LunTau,T2Name(aGrp,bGrp),1,1)

!3 make Tau

if (aGrp /= bGrp) then
  call GetTauHlp1(Tau,T1,dima,dimb,adda,addb,no,nv)
else
  call GetTauHlp2(Tau,T1,dima,adda,no,nv)
end if

return

end subroutine GetTau
