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
subroutine TRACID(T,LUCIN,LUCOUT,LUSC1,LUSC2,LUSC3,VEC1,VEC2)
! Transform CI vector on LUCIN with T matrix after
! Docent Malmquist's recipe. Place result as next vector on LUOUT
!
! The transformation is done as a sequence of one-electron transformations
!
! with each orbital transformation being
!
! Sum(k=0,2) ( 1/k! sum(n' /= n) S(n'n) E_{n'n} ) Tnn^N_n
!
! with Sn'n = T(n'n)/Tnn
!
! each transformation is

use CandS, only: ISSM, ISSPC
use lucia_data, only: I12, I_RES_AB, IDISK, IH1FORM, INT1, NTOOB
use Constants, only: One, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: T(*)
real(kind=wp), intent(_OUT_) :: VEC1(*), VEC2(*)
integer(kind=iwp), intent(in) :: LUCIN, LUCOUT, LUSC1, LUSC2, LUSC3
integer(kind=iwp) :: K, LBLK
real(kind=wp) :: TKK
#ifdef _DEBUGPRINT_
real(kind=wp) :: CNORM
real(kind=wp), external :: INPRDD
#endif

LBLK = -1
! Transfer vector on LUCIN to LUSC1
!    COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
call COPVCD(LUCIN,LUSC1,VEC1,1,LBLK)
! A bit of info for the sigma routine
I_RES_AB = 0
! Do the one-electron update
I12 = 1
! With 1-electron integrals in complete block form
IH1FORM = 2
! Transform each orbital separately
do K=1,NTOOB
  ! Place (T(P,K)/S(K,K)   in one-electron integral list
  !    T_ROW_TO_H(T,H,K)
  call T_ROW_TO_H(T,INT1,K,TKK)
  ! T_{kk}^Nk
  !    T_TO_NK_VEC(T,KORB,ISM,ISPC,LUCIN,LUCOUT,C)
  call T_TO_NK_VEC(TKK,K,ISSM,ISSPC,LUSC1,LUSC2,VEC1)
  call COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
# ifdef _DEBUGPRINT_
  write(u6,*) ' output from T_TO_NK'
  call WRTVCD(VEC1,LUSC1,1,LBLK)
# endif
  ! For each orbital calculate (1+T+1/2 T^2)|0>
  ! + T
  call MV7(VEC1,VEC2,LUSC1,LUSC2)
# ifdef _DEBUGPRINT_
  write(u6,*) ' Correction vector'
  call WRTVCD(VEC1,LUSC2,1,LBLK)
# endif
  call VECSMDP(VEC1,VEC2,One,One,LUSC1,LUSC2,LUSC3,1,LBLK)
  call COPVCD(LUSC3,LUSC1,VEC1,1,LBLK)
# ifdef _DEBUGPRINT_
  write(u6,*) ' Updated vector'
  call WRTVCD(VEC1,LUSC1,1,LBLK)
# endif
  ! + 1/2 T^2
  call MV7(VEC1,VEC2,LUSC2,LUSC3)
# ifdef _DEBUGPRINT_
  write(u6,*) ' Correction vector'
  call WRTVCD(VEC1,LUSC3,1,LBLK)
# endif
  call VECSMDP(VEC1,VEC2,One,Half,LUSC1,LUSC3,LUSC2,1,LBLK)
  ! and transfer back to LUSC1
  call COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
# ifdef _DEBUGPRINT_
  write(u6,*) ' Updated vector'
  call WRTVCD(VEC1,LUSC1,1,LBLK)
# endif
end do
! And transfer to LUCOUT
#ifdef _DEBUGPRINT_
CNORM = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
write(u6,*) ' Norm of transformed vector',CNORM
#endif
!write(u6,*) ' Transformed vector'
!call WRTVCD(VEC1,LUSC1,1,LBLK)
IDISK(LUSC1) = 0
!write(u6,*) ' LUCOUT LUSC1 = ',LUCOUT,LUSC1
call COPVCD(LUSC1,LUCOUT,VEC1,0,LBLK)

end subroutine TRACID
