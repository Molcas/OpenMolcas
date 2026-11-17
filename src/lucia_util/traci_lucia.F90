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

subroutine TRACI_Lucia(T,nT,CI_Vector,nCSF)
! Transform CI vector on LUHC with T matrix after
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

use stdalloc, only: mma_allocate, mma_deallocate
use CandS, only: ICSPC, ICSM, ISSM, ISSPC
use lucia_data, only: INT1, IH1FORM, I12, I_RES_AB, IDISK, IH1FORM, INT1, MXNTTS
use lucia_data, only: LUSC1, LUSC2, LUSC3, LUDIA, LUHC, LUC, IREFSM, NTOOB
use lucia_data, only: kvec3_length, VEC3, XISPSM, MXSOOB, PSSIGN
use Constants, only: Zero, One, Half, Two
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in):: nT, nCSF
real(kind=wp), intent(in) :: T(nT)
real(kind=wp), intent(inout) :: CI_Vector(nCSF)
integer(kind=iwp) :: K, I12_SAVE, I_RES_AB_SAVE, nREC, I_Dummy(1)
integer(kind=iwp), parameter :: LBLK=-1
integer(kind=iwp), allocatable :: lRec(:)
real(kind=wp), allocatable :: VEC1(:), VEC2(:), LH1SAVE(:)
real(kind=wp) :: TKK
integer(kind=iwp) LBLOCK, IH1FORM_SAVE
#ifdef _DEBUGPRINT_
real(kind=wp), external :: INPRDD
#endif

! Copy of one-electron integrals
call mma_allocate(LH1SAVE,nT,Label='LH1SAVE')
! We are going to mess with the one-electron integrals, take a copy
LH1SAVE(:) = INT1(:)
IH1FORM_SAVE = IH1FORM

LBLOCK = max(int(XISPSM(IREFSM,1)),MXSOOB)
if (PSSIGN /= Zero) LBLOCK = int(Two*XISPSM(IREFSM,1))

! The two scratch  blocks
call mma_allocate(VEC1,LBLOCK,Label='VEC1')
call mma_allocate(VEC2,LBLOCK,Label='VEC2')

Call mma_allocate(lRec,MXNTTS,Label='lRec')
call BLKFO_MIN(IREFSM,NREC,LREC)

call mma_allocate(VEC3,KVEC3_LENGTH,Label='VEC3')

IDISK(LUC) = 0
! Convert to internal Lucia SD format
call CSDTVC(CI_Vector,VEC1,1,IREFSM,0)
! Write CI vector to the Lucia file LUC
call TODSCN(VEC1,NREC,LREC,LBLK,LUC)
I_DUMMY(1) = -1
call ITODS(I_DUMMY,1,LBLK,LUC)
IDISK(LUC)=0
IDISK(LUSC1) = 0
! Copy the CI vector from file LUC to LUSC1
call COPVCD(LUC,LUSC1,VEC1,0,LBLK) ! From the current position
IDISK(LUC)=0
IDISK(LUHC)=0
call COPVCD(LUSC1,LUHC,VEC1,1,LBLK)! Rewind first
IDISK(LUDIA)=0

ICSPC = ISSPC
ICSM = ISSM
! Transfer vector on LUHC to LUSC1
call COPVCD(LUHC,LUSC1,VEC1,1,LBLK)

! A bit of info for the sigma routine
I_RES_AB_SAVE = I_RES_AB
I_RES_AB = 0
! Do one-electron update, i.e. don't use any two-electron integrals
I12_SAVE=I12
I12 = 1
! With 1-electron integrals in complete block form
IH1FORM = 2
! Transform each orbital separately

do K=1,NTOOB
  ! Place (T(P,K)/S(K,K)   in one-electron integral list
  call T_ROW_TO_H(T,INT1,K,TKK)
  ! T_{kk}^Nk
  call T_TO_NK_VEC(TKK,K,ISSM,ISSPC,LUSC1,LUSC2,VEC1)

  call COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)

  ! VEC1 = |Psi>
  ! For each orbital calculate (1+T+1/2 T^2)|0>
  ! + T
  ! VEC2 = T|Psi>
  call MV7(VEC1,VEC2,LUSC1,LUSC2,1,[Zero]) ! Dummy ERIs

  ! Form |Psi> + T|Psi>
  call VECSMDP(VEC1,VEC2,One,One,LUSC1,LUSC2,LUSC3,1,LBLK)

  call COPVCD(LUSC3,LUSC1,VEC1,1,LBLK)
  ! + 1/2 T^2
  call MV7(VEC1,VEC2,LUSC2,LUSC3,1,[Zero]) ! Dummy ERIs

  call VECSMDP(VEC1,VEC2,One,Half,LUSC1,LUSC3,LUSC2,1,LBLK)

  ! and transfer back to LUSC1
  call COPVCD(LUSC2,LUSC1,VEC1,1,LBLK)
end do

! And transfer to LUDIA
IDISK(LUSC1) = 0
call COPVCD(LUSC1,LUDIA,VEC1,0,LBLK)
IDISK(LUDIA)=0
! Copy CI-vector back to MOLCAS JOBIPH file
call FRMDSCN(VEC1,NREC,LBLK,LUDIA)
IDISK(LUDIA)=0
call CSDTVC(VEC2,VEC1,2,IREFSM,0)
call IFRMDS(I_DUMMY,1,LBLK,LUDIA)
IDISK(LUDIA)=0

I12=I12_SAVE
I_RES_AB = I_RES_AB_SAVE
call mma_deallocate(VEC3)

Call mma_deallocate(lRec)
CI_Vector(1:nCSF)=Vec2(1:nCSF)

call mma_deallocate(VEC2)
call mma_deallocate(VEC1)

! clean up time : copy 1-e integrals back in place
INT1(:) = LH1SAVE(:)
IH1FORM = IH1FORM_SAVE
call mma_deallocate(LH1SAVE)

end subroutine TRACI_Lucia
