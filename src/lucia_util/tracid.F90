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

use GLBBAS, only: INT1
use CandS, only: ISSM, ISSPC
use lucia_data, only: I12, I_RES_AB, IDISK, IH1FORM, NTOOB
use Constants, only: One, Half
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: T(*)
real(kind=wp), intent(_OUT_) :: VEC1(*), VEC2(*)
integer(kind=iwp), intent(in) :: LUCIN, LUCOUT, LUSC1, LUSC2, LUSC3
integer(kind=iwp) :: K, LBLK, NTEST
real(kind=wp) :: CNORM, INPRDD, TKK

NTEST = 0
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
  if (NTEST >= 1000) then
    write(u6,*) ' output from T_TO_NK'
    call WRTVCD(VEC1,LUSC1,1,LBLK)
  end if
  ! For each orbital calculate (1+T+1/2 T^2)|0>
  ! + T
  call MV7(VEC1,VEC2,LUSC1,LUSC2)
  if (NTEST >= 1000) then
    write(u6,*) ' Correction vector'
    call WRTVCD(VEC1,LUSC2,1,LBLK)
  end if
  call VECSMDP(VEC1,VEC2,One,One,LUSC1,LUSC2,LUSC3,1,LBLK)
  call COPVCD(LUSC3,LUSC1,VEC1,1,LBLK)
  if (NTEST >= 1000) then
    write(u6,*) ' Updated vector'
    call WRTVCD(VEC1,LUSC1,1,LBLK)
  end if
  ! + 1/2 T^2
  call MV7(VEC1,VEC2,LUSC2,LUSC3)
  if (NTEST >= 1000) then
    write(u6,*) ' Correction vector'
    call WRTVCD(VEC1,LUSC3,1,LBLK)
  end if
  call VECSMDP(VEC1,VEC2,One,Half,LUSC1,LUSC3,LUSC2,1,LBLK)
  ! and transfer back to LUSC1
  call COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
  if (NTEST >= 1000) then
    write(u6,*) ' Updated vector'
    call WRTVCD(VEC1,LUSC1,1,LBLK)
  end if
end do
! And transfer to LUCOUT
CNORM = INPRDD(VEC1,VEC2,LUSC1,LUSC1,1,LBLK)
if (NTEST > 0) write(u6,*) ' Norm of transformed vector',CNORM
!write(u6,*) ' Transformed vector'
!call WRTVCD(VEC1,LUSC1,1,LBLK)
IDISK(LUSC1) = 0
!write(u6,*) ' LUCOUT LUSC1 = ',LUCOUT,LUSC1
call COPVCD(LUSC1,LUCOUT,VEC1,0,LBLK)

end subroutine TRACID
