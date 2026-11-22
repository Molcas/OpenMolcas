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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!               2026, Roland Lindh                                     *
!***********************************************************************

!#define _DEBUGPRINT_
! Matrix elements <I|E_pq|J>, <I|E_pq E_rs|J> and <I|H|J> between CSFs
! of the split graph (walk-pair evaluation, coupling-block lookups).
! (Default public, as SGUGA_Core; the interface is the module SGUGA.)
module SGUGA_MatEl

use SGUGA_Core

implicit none

contains


!=======================================================================
! Matrix elements <I|E_pq|J>, <I|E_pq E_rs|J> and <I|H|J> between single
! CSFs (module level, used by SG_Test_EpqErs_Ref and by SG_MatEl_H). The
! state and the CSF symmetry are ME_iState and ME_LSM.
!=======================================================================

!
! Quick zero test for <I|E_pq E_rs|J>, from the orbital occupations OccI and
! OccJ of I and J (level by level). Returns .true. only if the matrix element is exactly zero:
!
!   1. |K> = E_rs|J> must exist: n_J(s) >= 1 and n_J(r) <= 1 for r /= s,
!      n_J(r) >= 1 for r == s (E_rr|J> = n_r|J>). Likewise E_pq|K> needs
!      n_K(p) >= 1 for p == q (for p /= q this follows from point 2).
!   2. n_I(t) = n_J(t) + d(t,p) - d(t,q) + d(t,r) - d(t,s) for every level t.
!
! A .false. result does not mean the element is nonzero (spin coupling can
! still make it vanish). Spatial symmetry is implied by point 2.
!
function Quick_Zero(OccI,OccJ,p,q,r,s) result(IsZero)

integer(kind=iwp), intent(in) :: OccI(:), OccJ(:), p, q, r, s
logical :: IsZero
integer(kind=iwp) :: t, nK, nI

IsZero = .true.

! Intermediate state K = E_rs|J>.
if (r == s) then
   if (OccJ(r) == 0) return
else
   if (OccJ(s) == 0) return
   if (OccJ(r) == 2) return
end if

! Occupations of I must equal those of E_pq E_rs |J>.
do t=1,SGS(ME_iState)%nLev
   nK = OccJ(t)
   if (t == r) nK = nK+1
   if (t == s) nK = nK-1
   if (p == q .and. t == p .and. nK == 0) return
   nI = nK
   if (t == p) nI = nI+1
   if (t == q) nI = nI-1
   if (nI /= OccI(t)) return
end do

IsZero = .false.

end function Quick_Zero

function MatEl_Upper_Half( &
           ISYUP,MV, &
           IUPI,IUPJ, &
           p,q ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYUP
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: IUPI
integer(kind=iwp), intent(in) :: IUPJ
integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q

real(kind=wp) :: Value

integer(kind=iwp) :: INDEO
integer(kind=iwp) :: NCP
integer(kind=iwp) :: IOFF
integer(kind=iwp) :: ICP
integer(kind=iwp) :: IFrom
integer(kind=iwp) :: ITo
integer(kind=iwp) :: iValue

logical :: Reverse

Value = Zero

if (p <= SGS(ME_iState)%MidLev) return
if (q <= SGS(ME_iState)%MidLev) return
if (p >  SGS(ME_iState)%nLev) return
if (q >  SGS(ME_iState)%nLev) return

!
! Diagonal upper-half operator.
!
if (p == q) then
   if (IUPI /= IUPJ) return
   Value = Upper_Diagonal_Value( &
              ISYUP,MV,IUPI,p )
   return
end if

Reverse = p < q

call Decode_Upper_Operator( &
     SGS(ME_iState),p,q,INDEO )

NCP = EXS(ME_iState)%NOCP(INDEO,ISYUP,MV)
if (NCP == 0) return

IOFF = EXS(ME_iState)%IOCP(INDEO,ISYUP,MV)

do ICP=1,NCP
   if (.not. Reverse) then
      IFrom = EXS(ME_iState)%ICoup(1,IOFF+ICP)
      ITo   = EXS(ME_iState)%ICoup(2,IOFF+ICP)
   else
      IFrom = EXS(ME_iState)%ICoup(2,IOFF+ICP)
      ITo   = EXS(ME_iState)%ICoup(1,IOFF+ICP)
   end if

   if (IFrom /= IUPJ) cycle
   if (ITo   /= IUPI) cycle

   iValue = EXS(ME_iState)%ICoup(3,IOFF+ICP)
   Value = Value + EXS(ME_iState)%VTab(iValue)
end do

end function MatEl_Upper_Half

function MatEl_Lower_Half( &
           ISYCI,ISYUP,MV, &
           IDWNI,IDWNJ, &
           p,q ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: ISYUP
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: IDWNI
integer(kind=iwp), intent(in) :: IDWNJ
integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q

real(kind=wp) :: Value

integer(kind=iwp) :: ISYDWN
integer(kind=iwp) :: INDEO
integer(kind=iwp) :: NCP
integer(kind=iwp) :: IOFF
integer(kind=iwp) :: ICP
integer(kind=iwp) :: IFrom
integer(kind=iwp) :: ITo
integer(kind=iwp) :: iValue

logical :: Reverse


Value = Zero

if (p < 1) return
if (q < 1) return
if (p > SGS(ME_iState)%MidLev) return
if (q > SGS(ME_iState)%MidLev) return

if (p == q) then
   if (IDWNI /= IDWNJ) return

   Value = Lower_Diagonal_Value( &
              ISYUP,MV,IDWNI,p )

   return
end if

!
! The lower coupling tables are indexed by the symmetry of
! the lower partial walk.
!
ISYDWN = MUL(ISYCI,ISYUP)

Reverse = p < q
call Decode_Lower_Operator( SGS(ME_iState),p,q,INDEO )

!
! Off-diagonal lower coupling paths are indexed by the
! symmetry carried by the lower partial walk.
!
NCP = EXS(ME_iState)%NOCP(INDEO,ISYDWN,MV)
if (NCP == 0) return

IOFF = EXS(ME_iState)%IOCP(INDEO,ISYDWN,MV)

do ICP=1,NCP


   if (Reverse) then
      IFrom = EXS(ME_iState)%ICoup(2,IOFF+ICP)
      ITo   = EXS(ME_iState)%ICoup(1,IOFF+ICP)
   else
      IFrom = EXS(ME_iState)%ICoup(1,IOFF+ICP)
      ITo   = EXS(ME_iState)%ICoup(2,IOFF+ICP)
   end if

   if (IFrom /= IDWNJ) cycle
   if (ITo   /= IDWNI) cycle

   iValue = EXS(ME_iState)%ICoup(3,IOFF+ICP)
   Value = Value + EXS(ME_iState)%VTab(iValue)

end do

end function MatEl_Lower_Half

function Lower_Diagonal_Value( &
           ISYUP,MV,IDWN,p ) result(Value)

use Constants,     only: Zero
use Definitions,   only: wp, iwp
use Symmetry_Info, only: Mul

implicit none

integer(kind=iwp), intent(in) :: ISYUP
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: IDWN
integer(kind=iwp), intent(in) :: p

real(kind=wp) :: Value

integer(kind=iwp) :: ISYDOWN
integer(kind=iwp) :: IOLW
integer(kind=iwp) :: IPSHFT
integer(kind=iwp) :: IPPOW
integer(kind=iwp) :: LLW
integer(kind=iwp) :: JC
integer(kind=iwp) :: ICS

Value = Zero

!
! Evaluate
!
!     <I_L|E_pp|I_L>
!
! for a lower-half walk.
!
if (p < 1) return
if (p > SGS(ME_iState)%MidLev) return

!
! The symmetry of the lower walk follows from the upper-walk
! symmetry and the symmetry of the complete CI state.
!
ISYDOWN = Mul(ISYUP,ME_LSM)

!
! Starting position of the lower walks in this split-graph block.
!
IOLW = CIS(ME_iState)%IOW(2,ISYDOWN,MV)

!
! Each step vector occupies two packed bits.
!
IPSHFT = 2*(p-1)

!
! Locate the packed integer containing the step vector at p.
!
LLW = 1 + IOLW - CIS(ME_iState)%nIpWlk + &
      IPSHFT/(2*nPack)

IPSHFT = mod(IPSHFT,2*nPack)
IPPOW  = 2**IPSHFT

!
! Extract the step vector belonging to lower walk IDWN.
!
JC = CIS(ME_iState)%ICase( &
        LLW + IDWN*CIS(ME_iState)%nIpWlk )

ICS = mod(JC/IPPOW,4)

!
! Diagonal segment values:
!
!     step 0 -> 0
!     step 1 -> 1
!     step 2 -> 1
!     step 3 -> 2
!
if (ICS /= 0) then
   Value = real((1+ICS)/2,kind=wp)
end if

end function Lower_Diagonal_Value

function MatEl_Lower_One( &
           ISYCI, &
           I,J, &
           p,q ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q

real(kind=wp) :: Value

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ
integer(kind=iwp) :: ISYDWN

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

integer(kind=iwp) :: INDEO
integer(kind=iwp) :: NCP
integer(kind=iwp) :: IOFF


logical :: Reverse

Value = Zero

!
! This helper handles one non-diagonal lower operator.
!
if (p == q) return

if (p < 1) return
if (q < 1) return

if (p > SGS(ME_iState)%MidLev) return
if (q > SGS(ME_iState)%MidLev) return

call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! A lower operator leaves the upper walk unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IUPI   /= IUPJ)   return

!
! Lower coupling tables are indexed by the symmetry carried
! by the lower partial walk.
!
ISYDWN = Mul(ISYCI,ISYUPI)

Reverse = p < q

call Decode_Lower_Operator( &
     SGS(ME_iState),p,q,INDEO )

NCP = EXS(ME_iState)%NOCP(INDEO,ISYDWN,MVI)
if (NCP == 0) return

IOFF = EXS(ME_iState)%IOCP(INDEO,ISYDWN,MVI)

! the transition IDWNJ -> IDWNI (walk index)
Value = Coup_Value(EXS(ME_iState),IOFF,NCP,merge(2,1,Reverse),IDWNJ,IDWNI)

end function MatEl_Lower_One


function Upper_Diagonal_Value( &
           ISYUP,MV,IUP,p ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYUP
integer(kind=iwp), intent(in) :: MV
integer(kind=iwp), intent(in) :: IUP
integer(kind=iwp), intent(in) :: p

real(kind=wp) :: Value

integer(kind=iwp) :: IOUW
integer(kind=iwp) :: IPSHFT
integer(kind=iwp) :: IPPOW
integer(kind=iwp) :: LUW
integer(kind=iwp) :: IC
integer(kind=iwp) :: ICS

Value = Zero

!
! This routine evaluates the upper-walk diagonal matrix element
!
!     <I_U|E_pp|I_U>.
!
if (p <= SGS(ME_iState)%MidLev) return
if (p >  SGS(ME_iState)%nLev)   return

!
! Starting position of the upper walks belonging to this
! split-graph block.
!
IOUW = CIS(ME_iState)%IOW(1,ISYUP,MV)

!
! Each step vector occupies two packed bits.
!
IPSHFT = 2*(p-1-SGS(ME_iState)%MidLev)

!
! Select the packed integer containing the step vector at level p.
!
LUW = 1 + IOUW - CIS(ME_iState)%nIpWlk + &
      IPSHFT/(2*nPack)

IPSHFT = mod(IPSHFT,2*nPack)
IPPOW  = 2**IPSHFT

!
! Extract the step vector of upper walk IUP at level p.
!
IC = CIS(ME_iState)%ICase( &
        LUW + IUP*CIS(ME_iState)%nIpWlk )

ICS = mod(IC/IPPOW,4)

!
! The diagonal segment values are:
!
!     step 0 -> 0
!     step 1 -> 1
!     step 2 -> 1
!     step 3 -> 2
!
if (ICS /= 0) then
   Value = real((1+ICS)/2,kind=wp)
end if

end function Upper_Diagonal_Value

function MatEl_Upper_One( &
           ISYCI, &
           I,J, &
           p,q ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q

real(kind=wp) :: Value

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

integer(kind=iwp) :: INDEO
integer(kind=iwp) :: NCP
integer(kind=iwp) :: IOFF


logical :: Reverse

Value = Zero

!
! This helper handles only a non-diagonal upper operator.
!
if (p == q) return

if (p <= SGS(ME_iState)%MidLev) return
if (q <= SGS(ME_iState)%MidLev) return

call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! An upper operator leaves the lower walk unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IDWNI  /= IDWNJ)  return

Reverse = p < q

call Decode_Upper_Operator( &
     SGS(ME_iState),p,q,INDEO )

NCP = EXS(ME_iState)%NOCP(INDEO,ISYUPI,MVI)

if (NCP == 0) return

IOFF = EXS(ME_iState)%IOCP(INDEO,ISYUPI,MVI)

! the transition IUPJ -> IUPI (walk index)
Value = Coup_Value(EXS(ME_iState),IOFF,NCP,merge(2,1,Reverse),IUPJ,IUPI)

end function MatEl_Upper_One

function MatEl_EpqErs_Ref(ISYCI,I,J,p,q,r,s) result(Value)

   integer(kind=iwp), intent(in) :: ISYCI,I,J
   integer(kind=iwp), intent(in) :: p,q,r,s

   real(kind=wp) :: Value

   integer(kind=iwp) :: nCSFI,nCSFK
#ifdef _DEBUGPRINT_
   integer(kind=iwp) :: K
   integer(kind=iwp) :: MVK,ISYUPK,IUPK,IDWNK
#endif
   integer(kind=iwp) :: ISYK

   real(kind=wp), allocatable :: Ket(:)
   real(kind=wp), allocatable :: Tmp(:)
   real(kind=wp), allocatable :: Res(:)
   logical(kind=iwp) :: Reuse_Save

   Reuse_Save=EXS(ME_iState)%ReUse_SGTMP
   EXS(ME_iState)%ReUse_SGTMP=.False.
   !----------------------------------------------------------
   ! Symmetry after application of E_rs.
   !
   ! |K> = E_rs |J>
   !
   ! Gamma(K)=Gamma(r)Gamma(s)Gamma(J)
   !----------------------------------------------------------

   ISYK = Mul(SGS(ME_iState)%ISm(r), &
              Mul(SGS(ME_iState)%ISm(s),ISYCI))

   nCSFI = CIS(ME_iState)%NCSF(ISYCI)
   nCSFK = CIS(ME_iState)%NCSF(ISYK)

   call mma_allocate(Ket,nCSFI,Label='Ket')
   call mma_allocate(Tmp,nCSFK,Label='Tmp')
   call mma_allocate(Res,nCSFI,Label='Res')

   Ket(:) = Zero
   Ket(J) = One

   Tmp(:) = Zero

   call SG_Epq_Psi( SGS(ME_iState), &
                    CIS(ME_iState), &
                    EXS(ME_iState), &
                    r,s,One,ISYCI, &
                    Ket,Tmp )

#ifdef _DEBUGPRINT_
   write(6,*)
   write(6,*) 'INTERMEDIATE STATES AFTER E_rs'
   write(6,*) 'I,J,p,q,r,s=',I,J,p,q,r,s
   write(6,*) 'ISYK=',ISYK

   do K=1,nCSFK

      if (abs(Tmp(K)) > 1.0d-12) then

         call Decode_CSF( CIS(ME_iState),ISYK,K, &
                          MVK,ISYUPK,IUPK,IDWNK )

         write(6,'(I6,F16.8,4I6)') &
              K,Tmp(K),MVK,ISYUPK,IUPK,IDWNK

      end if

   end do
#endif

   Res(:) = Zero

   call SG_Epq_Psi( SGS(ME_iState), &
                    CIS(ME_iState), &
                    EXS(ME_iState), &
                    p,q,One,ISYK, &
                    Tmp,Res )

   Value = Res(I)

   call mma_deallocate(Ket)
   call mma_deallocate(Tmp)
   call mma_deallocate(Res)

   EXS(ME_iState)%ReUse_SGTMP=Reuse_Save

end function MatEl_EpqErs_Ref


function MatEl_LL( &
           I,J, &
           INDEO_Left,Reverse_Left, &
           INDEO_Right,Reverse_Right, &
           ISYRS ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: INDEO_Left
integer(kind=iwp), intent(in) :: INDEO_Right

logical, intent(in) :: Reverse_Left
logical, intent(in) :: Reverse_Right

! Symmetry of the right operator, Gamma(r)Gamma(s)
integer(kind=iwp), intent(in) :: ISYRS

real(kind=wp) :: Value

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ
integer(kind=iwp) :: ISYDWN

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ




Value = Zero

call Decode_ME( &
     ME_LSM, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ME_LSM, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! Both operators act entirely on the lower half.
! The upper half-walk is a spectator.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IUPI   /= IUPJ)   return

!
! For
!
!     <I|E_Left E_Right|J>,
!
! the rightmost operator acts first:
!
!     J_L --E_Right--> K_L --E_Left--> I_L.
!
!
! Lower coupling blocks are indexed by the symmetry carried
! by the lower partial walk, not by the upper-walk symmetry.
!
ISYDWN = Mul(ME_LSM,ISYUPI)

Value = Walk_Path_Value( &
           INDEO_Right,Reverse_Right,merge(Mul(ISYRS,ISYDWN),ISYDWN,Reverse_Right), &
           INDEO_Left, Reverse_Left, merge(ISYDWN,Mul(ISYRS,ISYDWN),Reverse_Left), &
           MVI,IDWNJ,IDWNI )

end function MatEl_LL

function MatEl_UU( &
           ISYCI, &
           I,J, &
           INDEO_Left,Reverse_Left, &
           INDEO_Right,Reverse_Right, &
           ISYRS ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: INDEO_Left
integer(kind=iwp), intent(in) :: INDEO_Right

logical, intent(in) :: Reverse_Left
logical, intent(in) :: Reverse_Right

! Symmetry of the right operator, Gamma(r)Gamma(s)
integer(kind=iwp), intent(in) :: ISYRS

real(kind=wp) :: Value

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ




Value = Zero

call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! Both operators act entirely on the upper half.
! The lower half-walk must therefore be unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IDWNI  /= IDWNJ)  return

!
! For
!
!     <I|E_Left E_Right|J>,
!
! the rightmost operator acts first:
!
!     J_U --E_Right--> K_U --E_Left--> I_U.
!
Value = Walk_Path_Value( &
           INDEO_Right,Reverse_Right,merge(Mul(ISYRS,ISYUPI),ISYUPI,Reverse_Right), &
           INDEO_Left, Reverse_Left, merge(ISYUPI,Mul(ISYRS,ISYUPI),Reverse_Left), &
           MVI,IUPJ,IUPI )

end function MatEl_UU

!
! Sum over the intermediate walk K of the coupling-coefficient products
! along J --E_Right--> K --E_Left--> I, for the walks IFromJ and IToI of a
! half graph (upper or lower) at midvertex MV. The blocks of both operators
! are given by INDEO, orientation (Reverse: transposed stored transitions)
! and the symmetry under which the block is stored (that of its ICoup(1)
! walk). Only the transitions starting at IFromJ are followed, located
! by the walk index (Path_Value).
!
function Walk_Path_Value(INDEO_Right,Reverse_Right,ISYRight,INDEO_Left,Reverse_Left,ISYLeft,MV,IFromJ,IToI) result(Value)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: INDEO_Right, ISYRight, INDEO_Left, ISYLeft, MV, IFromJ, IToI
logical, intent(in) :: Reverse_Right, Reverse_Left
real(kind=wp) :: Value
integer(kind=iwp) :: NRight, NLeft

NRight = EXS(ME_iState)%NOCP(INDEO_Right,ISYRight,MV)
NLeft  = EXS(ME_iState)%NOCP(INDEO_Left, ISYLeft, MV)
Value = Path_Value(EXS(ME_iState), &
                   EXS(ME_iState)%IOCP(INDEO_Right,ISYRight,MV),NRight,merge(2,1,Reverse_Right),IFromJ, &
                   EXS(ME_iState)%IOCP(INDEO_Left, ISYLeft, MV),NLeft, merge(2,1,Reverse_Left), IToI)

end function Walk_Path_Value

subroutine Decode_Upper_Operator(SGS,p,q,INDEO)

type(SGStruct), intent(in) :: SGS
integer(kind=iwp), intent(in)  :: p,q
integer(kind=iwp), intent(out) :: INDEO

logical(kind=iwp) :: Swap

Swap = q > p

INDEO = merge( 2*SGS%nLev + (q*(q-1))/2 + p, &
               2*SGS%nLev + (p*(p-1))/2 + q, &
               Swap )

end subroutine Decode_Upper_Operator

subroutine Decode_Lower_Operator(SGS,p,q,INDEO)

implicit none

type(SGStruct), intent(in) :: SGS

integer(kind=iwp), intent(in)  :: p
integer(kind=iwp), intent(in)  :: q
integer(kind=iwp), intent(out) :: INDEO

integer(kind=iwp) :: iHigh
integer(kind=iwp) :: iLow

!
! The closed-loop table is indexed by the unordered orbital pair.
! Store the larger orbital index in iHigh and the smaller one
! in iLow.
!
iHigh = max(p,q)
iLow  = min(p,q)

INDEO = 2*SGS%nLev + &
        (iHigh*(iHigh-1))/2 + iLow

end subroutine Decode_Lower_Operator

function MatEl_EpqErs_UU( &
           ISYCI, &
           I,J, &
           p,q,r,s ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI

integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q
integer(kind=iwp), intent(in) :: r
integer(kind=iwp), intent(in) :: s

integer(kind=iwp) :: INDEO_Left
integer(kind=iwp) :: INDEO_Right

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

logical :: Diagonal_Left
logical :: Diagonal_Right

logical :: Reverse_Left
logical :: Reverse_Right

real(kind=wp) :: XClosed
real(kind=wp) :: WLeft
real(kind=wp) :: WRight

real(kind=wp) :: Value

Value = Zero

!
! Both operators must lie completely in the upper graph.
!
if (p <= SGS(ME_iState)%MidLev) return
if (q <= SGS(ME_iState)%MidLev) return
if (r <= SGS(ME_iState)%MidLev) return
if (s <= SGS(ME_iState)%MidLev) return

Diagonal_Left  = p == q
Diagonal_Right = r == s

!
! Decode the bra and ket CSFs. The decoded upper walks are
! required for diagonal weights.
!
call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! Every UU operator product leaves the lower walk unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IDWNI  /= IDWNJ)  return

!
! ============================================================
! Case 1: diagonal x diagonal
!
!     <I|E_pp E_rr|J>
!
! Both operators are diagonal, so the complete upper walk must
! also be unchanged.
! ============================================================
!
if (Diagonal_Left .and. Diagonal_Right) then

   if (IUPI /= IUPJ) return

   WLeft = Upper_Diagonal_Value( &
              ISYUPI,MVI,IUPI,p )

   if (WLeft == Zero) return

   WRight = Upper_Diagonal_Value( &
               ISYUPI,MVI,IUPI,r )

   Value = WLeft*WRight

   return

end if

!
! ============================================================
! Case 2: diagonal x closed
!
!     <I|E_pp E_rs|J>
!
! Since E_pp is the left operator, its weight is evaluated on
! the final bra-side upper walk I_U:
!
!     w_p(I_U) <I|E_rs|J>.
! ============================================================
!
if (Diagonal_Left) then

   WLeft = Upper_Diagonal_Value( &
              ISYUPI,MVI,IUPI,p )

   if (WLeft == Zero) return

   XClosed = MatEl_Upper_One( &
                ISYCI, &
                I,J, &
                r,s )

   Value = WLeft*XClosed

   return

end if

!
! ============================================================
! Case 3: closed x diagonal
!
!     <I|E_pq E_rr|J>
!
! Since E_rr is the rightmost operator, its weight is evaluated
! on the initial ket-side upper walk J_U:
!
!     <I|E_pq|J> w_r(J_U).
! ============================================================
!
if (Diagonal_Right) then

   XClosed = MatEl_Upper_One( &
                ISYCI, &
                I,J, &
                p,q )

   if (XClosed == Zero) return

   WRight = Upper_Diagonal_Value( &
               ISYUPJ,MVJ,IUPJ,r )

   Value = XClosed*WRight

   return

end if

!
! ============================================================
! Case 4: closed x closed
!
! Use the already validated general ordered contraction.
! ============================================================
!
Reverse_Left  = p < q
Reverse_Right = r < s

call Decode_Upper_Operator( &
     SGS(ME_iState),p,q,INDEO_Left )

call Decode_Upper_Operator( &
     SGS(ME_iState),r,s,INDEO_Right )

Value = MatEl_UU( &
           ISYCI, &
           I,J, &
           INDEO_Left, Reverse_Left, &
           INDEO_Right,Reverse_Right, &
           Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)) )

end function MatEl_EpqErs_UU

function MatEl_EpqErs_LL( &
           ISYCI, &
           I,J, &
           p,q,r,s ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J

integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q
integer(kind=iwp), intent(in) :: r
integer(kind=iwp), intent(in) :: s

integer(kind=iwp) :: INDEO_Left
integer(kind=iwp) :: INDEO_Right

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ

integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ

integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ

integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

logical :: Diagonal_Left
logical :: Diagonal_Right

logical :: Reverse_Left
logical :: Reverse_Right

real(kind=wp) :: XClosed
real(kind=wp) :: WLeft
real(kind=wp) :: WRight

real(kind=wp) :: Value

Value = Zero

!
! Both operators must lie completely inside the lower graph.
!
if (p < 1 .or. q < 1) return
if (r < 1 .or. s < 1) return

if (p > SGS(ME_iState)%MidLev) return
if (q > SGS(ME_iState)%MidLev) return

if (r > SGS(ME_iState)%MidLev) return
if (s > SGS(ME_iState)%MidLev) return

Diagonal_Left  = p == q
Diagonal_Right = r == s

!
! Decode the bra and ket CSFs.
!
call Decode_ME( &
     ISYCI, &
     I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI, &
     J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! Every LL product leaves the upper walk unchanged.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= ISYUPJ) return
if (IUPI   /= IUPJ)   return

!
! ============================================================
! Case 1: diagonal x diagonal
!
!     <I|E_pp E_rr|J>
!
! Both operators are diagonal. Consequently, the lower walk
! must also remain unchanged.
! ============================================================
!
if (Diagonal_Left .and. Diagonal_Right) then

   if (IDWNI /= IDWNJ) return

   WLeft = Lower_Diagonal_Value( &
              ISYUPI,MVI,IDWNI,p )

   if (WLeft == Zero) return

   WRight = Lower_Diagonal_Value( &
               ISYUPI,MVI,IDWNI,r )

   Value = WLeft*WRight

   return

end if

!
! ============================================================
! Case 2: diagonal x closed
!
!     <I|E_pp E_rs|J>
!
! E_pp is the left operator. Its weight is therefore evaluated
! on the final bra-side lower walk:
!
!     w_p(I_L) <I|E_rs|J>.
! ============================================================
!
if (Diagonal_Left) then

   WLeft = Lower_Diagonal_Value( &
              ISYUPI,MVI,IDWNI,p )

   if (WLeft == Zero) return

   XClosed = MatEl_Lower_One( &
                ISYCI, &
                I,J, &
                r,s )

   Value = WLeft*XClosed

   return

end if

!
! ============================================================
! Case 3: closed x diagonal
!
!     <I|E_pq E_rr|J>
!
! E_rr is the rightmost operator and acts first. Its weight is
! evaluated on the initial ket-side lower walk:
!
!     <I|E_pq|J> w_r(J_L).
! ============================================================
!
if (Diagonal_Right) then

   XClosed = MatEl_Lower_One( &
                ISYCI, &
                I,J, &
                p,q )

   if (XClosed == Zero) return

   WRight = Lower_Diagonal_Value( &
               ISYUPJ,MVJ,IDWNJ,r )

   Value = XClosed*WRight

   return

end if

!
! ============================================================
! Case 4: closed x closed
!
! Use the already validated oriented LL contraction.
! ============================================================
!
Reverse_Left  = p < q
Reverse_Right = r < s

call Decode_Lower_Operator( &
     SGS(ME_iState),p,q,INDEO_Left )

call Decode_Lower_Operator( &
     SGS(ME_iState),r,s,INDEO_Right )

Value = MatEl_LL( &
           I,J, &
           INDEO_Left, Reverse_Left, &
           INDEO_Right,Reverse_Right, &
           Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)) )

end function MatEl_EpqErs_LL

!
! Compute <I|EpqErs|J> as Sum_K <I|Epq|K><K|Ers|J>
!
! Where <I|Epq|K> is an upper closed loop: Upq(Iu,Ku) x Delta(Il,Kl)
! and   <K|Ers|J> is an open loop: A_r(Ku,Ju) x A_s(Kl,Jl)
!
! Which reduce to
!
!  A_s(Il,Jl) x Sum_Ku Upq(Iu,Ku) A_r(Ku,Ju)
!
function MatEl_EpqErs_UO(ISYCJ,I,J,p,q,r_in,s_in) result(Value)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6
use Symmetry_Info, only: Mul

implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r_in,s_in

real(kind=wp) :: Value
real(kind=wp) :: XUpper,XLower,XRef

integer(kind=iwp) :: r,s
integer(kind=iwp) :: MVI,MVJ,MVK,MVX,MV
integer(kind=iwp) :: ISYUI,ISYUJ,ISYUK
integer(kind=iwp) :: IUI,IUJ
integer(kind=iwp) :: IDI,IDJ
integer(kind=iwp) :: ISYDI,ISYDJ
integer(kind=iwp) :: iBand
integer(kind=iwp) :: INDEO_R,INDEO_S,INDEO_PQ
integer(kind=iwp) :: NOU,NOL,NU
integer(kind=iwp) :: OffOU,OffOL,OffU
integer(kind=iwp) :: ISYCI, ISYP, ISYQ, ISYPQ, ISYR, ISYS, ISYRS, ISYPQRS

logical(kind=iwp) :: Reverse_U
logical(kind=iwp) :: Swap_rs, Verbose

Verbose=.False.
Value=Zero

Swap_rs = s_in > r_in
r = merge(s_in,r_in,Swap_rs)
s = merge(r_in,s_in,Swap_rs)

ISYP =SGS(ME_iState)%ISm(p)
ISYQ =SGS(ME_iState)%ISm(q)
ISYPQ=Mul(ISYP,ISYQ)
ISYR =SGS(ME_iState)%ISm(r)
ISYS =SGS(ME_iState)%ISm(s)
ISYRS=Mul(ISYR,ISYS)
ISYPQRS=Mul(ISYPQ,ISYRS)
ISYCI = Mul(ISYPQRS,ISYCJ)

if (I < 1 .or. I > CIS(ME_iState)%NCSF(ISYCI)) return
if (J < 1 .or. J > CIS(ME_iState)%NCSF(ISYCJ)) return

if (p<=SGS(ME_iState)%MidLev .or. q<=SGS(ME_iState)%MidLev) return
if (r<=SGS(ME_iState)%MidLev .or. s<1 .or. s>SGS(ME_iState)%MidLev) return
if (max(p,q,r)>SGS(ME_iState)%nLev) return

call Decode_ME(ISYCI,I,MVI,ISYUI,IUI,IDI)
call Decode_ME(ISYCJ,J,MVJ,ISYUJ,IUJ,IDJ)
ISYUK=Mul(ISYPQ,ISYUI)
ISYDI=Mul(ISYCI,ISYUI)
ISYDJ=Mul(ISYCJ,ISYUJ)

If (Verbose) Then
   Write (u6,*) '        ISYP, ISYQ  =', ISYP, ISYQ
   Write (u6,*) '        MVI, MVJ  =', MVI, MVJ
   Write (u6,*) '        IUP, JUP  =', IUI, IUJ
   Write (u6,*) '        IDN, JDN  =', IDI, IDJ
   Write (u6,*) '        ISYCI, ISYCJ=', ISYCI, ISYCJ
End If


if (p==q) then
   ! <I|EppErs|J> = Sum_K <I|Epp|K><K|Ers|J>
   !
   ! Epp|K> = np|K> gives
   !
   ! <I|EppErs> = <I|Epp|I> <I|Ers|J>
   Value=MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),p,p,ISYCI,I,I,.False.)* &
         MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),r_in,s_in,ISYCJ,I,J,.False.)
   if (SG_Checks) then  ! sanity check, only when verifying CI developments
   XRef=MatEl_EpqErs_Ref(ISYCJ,I,J,p,q,r_in,s_in)
   if (abs(Value-XRef) > 1.0e-12_wp) then
      write(u6,*) '<I|EppErs|J> discrepancy:',I,J,p,r_in,s_in,Value,XRef
      call Abend()
   end if
   end if
   return
end if

! For the open-loop, A_r(Ku(m),Ju(m')) and m=MVK find the connecting m'=MVJ
! Since we have delta(Il,Kl) from the closed loop this implies that MVK=MVI.

MVK=MVI
iBand=0
do MV = 1, 2
   MVX=MERGE(EXS(ME_iState)%MVR(MVK,MV),EXS(ME_iState)%MVL(MVK,MV),Swap_rs)
   If (MVX == 0 .or. MVX /= MVJ) cycle
   iBand = MV
end do

If (iBand==0) Return
MV=iBand

INDEO_R=(MV-1)*SGS(ME_iState)%nLev+r  ! Open-Loop Upper

INDEO_S=(MV-1)*SGS(ME_iState)%nLev+s  ! Open-Loop Lower
If (Verbose) Write (u6,*) 'INDEO_R,INDEO_S=',INDEO_R,INDEO_S

! Symmetry relations for the two canonical open-loop factors.

! The open loop E_rs connects J to the intermediate K (not to I), whose
! upper walk has the symmetry ISYUK = Gamma(pq) Gamma(Iu).
if (ISYUJ /= Mul(SGS(ME_iState)%ISm(r),ISYUK)) return

ISYDI=Mul(ISYCI,ISYUI)
ISYDJ=Mul(ISYCJ,ISYUJ)

if (ISYDJ /= Mul(SGS(ME_iState)%ISm(s),ISYDI)) return

! ---------------------------------------------------------------
! Upper contraction
!
!   XUpper(Iu(m),Ju(m') =
!
!     sum_Ku U_pq(Iu(m),Ku(m)) x A_r(Ku(m),Ju(m'))
! ---------------------------------------------------------------

! Process: A_r(Ku(m),Ju(m'))
NOU   = MERGE(EXS(ME_iState)%NOCP(INDEO_R,ISYUK,MVK), &
              EXS(ME_iState)%NOCP(INDEO_R,ISYUJ,MVJ),Swap_rs)
OffOU = MERGE(EXS(ME_iState)%IOCP(INDEO_R,ISYUK,MVK), &
              EXS(ME_iState)%IOCP(INDEO_R,ISYUJ,MVJ),Swap_rs)
If (Verbose) Write (u6,*) 'NOU=',NOU
if (NOU==0) return

! Retrieve the upper closed-loop operator E_pq: Upq(Iu,Ku)
Reverse_U=p<q

call Decode_Upper_Operator(SGS(ME_iState),p,q,INDEO_PQ)

! A closed-loop block is indexed by the symmetry of its ICoup(1) (left)
! walk: K in canonical orientation, I if reversed (they differ for
! Gamma(pq) /= 1).
NU  =EXS(ME_iState)%NOCP(INDEO_PQ,merge(ISYUI,ISYUK,Reverse_U),MVI)
OffU=EXS(ME_iState)%IOCP(INDEO_PQ,merge(ISYUI,ISYUK,Reverse_U),MVI)
if (NU==0) return
If (Verbose) Write (u6,*) 'NU=',NU

! Sum_Ku U_pq(Iu,Ku) A_r(Ku,Ju): only the transitions of A_r starting at
! Ju and of U_pq starting at Ku (walk index, Path_Value).
XUpper=Path_Value(EXS(ME_iState),OffOU,NOU,merge(2,1,Swap_rs),IUJ,OffU,NU,merge(2,1,Reverse_U),IUI)
If (Verbose) Write (u6,*) 'XUpper=',XUpper
if (XUpper==Zero) return

! ---------------------------------------------------------------
! Lower canonical open-loop factor
!
!   XLower = A_s(Il(m),Jl(m')).
!
! ---------------------------------------------------------------

NOL  = MERGE(EXS(ME_iState)%NOCP(INDEO_S,ISYDI,MVI), &
             EXS(ME_iState)%NOCP(INDEO_S,ISYDJ,MVJ),Swap_rs)
OffOL= MERGE(EXS(ME_iState)%IOCP(INDEO_S,ISYDI,MVI), &
             EXS(ME_iState)%IOCP(INDEO_S,ISYDJ,MVJ),Swap_rs)
If (NOL==0) return
If (Verbose) Write (u6,*) 'NOL=',NOL

XLower=Coup_Value(EXS(ME_iState),OffOL,NOL,merge(2,1,Swap_rs),IDJ,IDI)
If (Verbose) Write (u6,*) 'XLower=',XLower

Value=XUpper*XLower
If (Verbose) Write (u6,*) 'Value=',Value

end function MatEl_EpqErs_UO

!
! <I|E_pq E_rs|J> with E_pq open (crossing the midlevel) and E_rs in the
! upper half.
!
! For real CI vectors <I|E_pq E_rs|J> = <J|(E_pq E_rs)^+|I> = <J|E_sr E_qp|I>,
! and E_sr E_qp is an upper x open product, so this is MatEl_EpqErs_UO with
! bra and ket exchanged and the indices (p,q,r,s) -> (s,r,q,p). (Unlike the
! UL/LU case the two operators may share an upper level, so they do not
! commute in general.)
!
! MatEl_EpqErs_UO takes the symmetry of its ket, here I:
! Gamma(I) = Gamma(p) Gamma(q) Gamma(r) Gamma(s) Gamma(J).
!
function MatEl_EpqErs_OU(ISYCJ,I,J,p,q,r,s) result(Value)

use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r,s
real(kind=wp) :: Value
integer(kind=iwp) :: ISYCI

ISYCI = Mul(Mul(SGS(ME_iState)%ISm(p),SGS(ME_iState)%ISm(q)), &
            Mul(Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)),ISYCJ))
Value = MatEl_EpqErs_UO(ISYCI,J,I,s,r,q,p)

end function MatEl_EpqErs_OU

!
! Compute <I|EpqErs|J> as Sum_K <I|Epq|K><K|Ers|J>
!
! Where <I|Epq|K> is a lower closed loop: Lpq(Il,Kl) x Delta(Iu,Ku)
! and   <K|Ers|J> is an open loop: A_r(Ku,Ju) x A_s(Kl,Jl)
!
! Which reduce to
!
!  A_r(Iu,Ju) x Sum_Kl Lpq(Il,Kl) A_s(Kl,Jl)
!
! This is the mirror image of MatEl_EpqErs_UO with the roles of the upper
! and lower halves exchanged. ISYCJ is the symmetry of the ket J.
!
function MatEl_EpqErs_LO(ISYCJ,I,J,p,q,r_in,s_in) result(Value)

use Constants, only: Zero
use Definitions, only: wp, iwp, u6
use Symmetry_Info, only: Mul

implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r_in,s_in

real(kind=wp) :: Value
real(kind=wp) :: XUpper,XLower,XRef

integer(kind=iwp) :: r,s
integer(kind=iwp) :: MVI,MVJ,MVK,MVX,MV
integer(kind=iwp) :: ISYUI,ISYUJ,ISYUK
integer(kind=iwp) :: IUI,IUJ
integer(kind=iwp) :: IDI,IDJ
integer(kind=iwp) :: ISYDI,ISYDJ,ISYDK
integer(kind=iwp) :: iBand
integer(kind=iwp) :: INDEO_R,INDEO_S,INDEO_PQ
integer(kind=iwp) :: NOU,NOL,NL
integer(kind=iwp) :: OffOU,OffOL,OffL
integer(kind=iwp) :: ISYCI, ISYP, ISYQ, ISYPQ, ISYR, ISYS, ISYRS, ISYPQRS

logical(kind=iwp) :: Reverse_L
logical(kind=iwp) :: Swap_rs, Verbose

Verbose=.False.
Value=Zero

! Canonical open loop: r in the upper, s in the lower half.
Swap_rs = s_in > r_in
r = merge(s_in,r_in,Swap_rs)
s = merge(r_in,s_in,Swap_rs)

ISYP =SGS(ME_iState)%ISm(p)
ISYQ =SGS(ME_iState)%ISm(q)
ISYPQ=Mul(ISYP,ISYQ)
ISYR =SGS(ME_iState)%ISm(r)
ISYS =SGS(ME_iState)%ISm(s)
ISYRS=Mul(ISYR,ISYS)
ISYPQRS=Mul(ISYPQ,ISYRS)
ISYCI = Mul(ISYPQRS,ISYCJ)

if (I < 1 .or. I > CIS(ME_iState)%NCSF(ISYCI)) return
if (J < 1 .or. J > CIS(ME_iState)%NCSF(ISYCJ)) return

if (p<1 .or. q<1 .or. p>SGS(ME_iState)%MidLev .or. q>SGS(ME_iState)%MidLev) return
if (r<=SGS(ME_iState)%MidLev .or. r>SGS(ME_iState)%nLev .or. s<1 .or. s>SGS(ME_iState)%MidLev) return

call Decode_ME(ISYCI,I,MVI,ISYUI,IUI,IDI)
call Decode_ME(ISYCJ,J,MVJ,ISYUJ,IUJ,IDJ)
ISYDI=Mul(ISYCI,ISYUI)
ISYDJ=Mul(ISYCJ,ISYUJ)

! The intermediate K has the upper walk and the midvertex of I, and the
! lower walk symmetry Gamma(Kl) = Gamma(pq) Gamma(Il).
MVK  =MVI
ISYUK=ISYUI
ISYDK=Mul(ISYPQ,ISYDI)

If (Verbose) Then
   Write (u6,*) '        ISYP, ISYQ  =', ISYP, ISYQ
   Write (u6,*) '        MVI, MVJ  =', MVI, MVJ
   Write (u6,*) '        IUP, JUP  =', IUI, IUJ
   Write (u6,*) '        IDN, JDN  =', IDI, IDJ
   Write (u6,*) '        ISYCI, ISYCJ=', ISYCI, ISYCJ
End If

if (p==q) then
   ! E_pp is diagonal in the CSF basis: <I|EppErs|J> = <I|Epp|I> <I|Ers|J>
   Value=MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),p,p,ISYCI,I,I,.False.)* &
         MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),r_in,s_in,ISYCJ,I,J,.False.)
   if (SG_Checks) then  ! sanity check, only when verifying CI developments
   XRef=MatEl_EpqErs_Ref(ISYCJ,I,J,p,q,r_in,s_in)
   if (abs(Value-XRef) > 1.0e-12_wp) then
      write(u6,*) '<I|EppErs|J> discrepancy (LO):',I,J,p,r_in,s_in,Value,XRef
      call Abend()
   end if
   end if
   return
end if

! For the open loop A_r(Ku(m),Ju(m')) with m=MVK find the connecting m'=MVJ.
iBand=0
do MV = 1, 2
   MVX=MERGE(EXS(ME_iState)%MVR(MVK,MV),EXS(ME_iState)%MVL(MVK,MV),Swap_rs)
   If (MVX == 0 .or. MVX /= MVJ) cycle
   iBand = MV
end do
If (iBand==0) Return
MV=iBand

INDEO_R=(MV-1)*SGS(ME_iState)%nLev+r  ! Open-Loop Upper
INDEO_S=(MV-1)*SGS(ME_iState)%nLev+s  ! Open-Loop Lower
If (Verbose) Write (u6,*) 'INDEO_R,INDEO_S=',INDEO_R,INDEO_S

! Symmetry relations for the two open-loop factors.
if (ISYUJ /= Mul(ISYR,ISYUK)) return
if (ISYDJ /= Mul(ISYS,ISYDK)) return

! ---------------------------------------------------------------
! Upper open-loop factor
!
!   XUpper = A_r(Iu(m),Ju(m'))      (Ku = Iu)
! ---------------------------------------------------------------

NOU  = MERGE(EXS(ME_iState)%NOCP(INDEO_R,ISYUK,MVK), &
             EXS(ME_iState)%NOCP(INDEO_R,ISYUJ,MVJ),Swap_rs)
OffOU= MERGE(EXS(ME_iState)%IOCP(INDEO_R,ISYUK,MVK), &
             EXS(ME_iState)%IOCP(INDEO_R,ISYUJ,MVJ),Swap_rs)
If (Verbose) Write (u6,*) 'NOU=',NOU
If (NOU==0) return

XUpper=Coup_Value(EXS(ME_iState),OffOU,NOU,merge(2,1,Swap_rs),IUJ,IUI)
If (Verbose) Write (u6,*) 'XUpper=',XUpper
if (XUpper==Zero) return

! ---------------------------------------------------------------
! Lower contraction
!
!   XLower = sum_Kl L_pq(Il(m),Kl(m)) x A_s(Kl(m),Jl(m'))
! ---------------------------------------------------------------

! Process: A_s(Kl(m),Jl(m'))
NOL  = MERGE(EXS(ME_iState)%NOCP(INDEO_S,ISYDK,MVK), &
             EXS(ME_iState)%NOCP(INDEO_S,ISYDJ,MVJ),Swap_rs)
OffOL= MERGE(EXS(ME_iState)%IOCP(INDEO_S,ISYDK,MVK), &
             EXS(ME_iState)%IOCP(INDEO_S,ISYDJ,MVJ),Swap_rs)
If (Verbose) Write (u6,*) 'NOL=',NOL
if (NOL==0) return

! The lower closed-loop operator E_pq: Lpq(Il,Kl). As in MatEl_Epq, the
! table is indexed by the bra (Il) symmetry for p < q and by the ket (Kl)
! symmetry otherwise.
Reverse_L=p<q
call Decode_Lower_Operator(SGS(ME_iState),p,q,INDEO_PQ)
NL  =merge(EXS(ME_iState)%NOCP(INDEO_PQ,ISYDI,MVI), &
           EXS(ME_iState)%NOCP(INDEO_PQ,ISYDK,MVI),Reverse_L)
OffL=merge(EXS(ME_iState)%IOCP(INDEO_PQ,ISYDI,MVI), &
           EXS(ME_iState)%IOCP(INDEO_PQ,ISYDK,MVI),Reverse_L)
If (Verbose) Write (u6,*) 'NL=',NL
if (NL==0) return

! Sum_Kl L_pq(Il,Kl) A_s(Kl,Jl): walk index, Path_Value.
XLower=Path_Value(EXS(ME_iState),OffOL,NOL,merge(2,1,Swap_rs),IDJ,OffL,NL,merge(2,1,Reverse_L),IDI)
If (Verbose) Write (u6,*) 'XLower=',XLower

Value=XUpper*XLower
If (Verbose) Write (u6,*) 'Value=',Value

end function MatEl_EpqErs_LO

!
! <I|E_pq E_rs|J> with E_pq open (crossing the midlevel) and E_rs in the
! lower half.
!
! For real CI vectors <I|E_pq E_rs|J> = <J|E_sr E_qp|I>, and E_sr E_qp is a
! lower x open product, so this is MatEl_EpqErs_LO with bra and ket
! exchanged and the indices (p,q,r,s) -> (s,r,q,p). The operators may share
! a lower level, so they do not commute in general.
!
! MatEl_EpqErs_LO takes the symmetry of its ket, here I:
! Gamma(I) = Gamma(p) Gamma(q) Gamma(r) Gamma(s) Gamma(J).
!
function MatEl_EpqErs_OL(ISYCJ,I,J,p,q,r,s) result(Value)

use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r,s
real(kind=wp) :: Value
integer(kind=iwp) :: ISYCI

ISYCI = Mul(Mul(SGS(ME_iState)%ISm(p),SGS(ME_iState)%ISm(q)), &
            Mul(Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)),ISYCJ))
Value = MatEl_EpqErs_LO(ISYCI,J,I,s,r,q,p)

end function MatEl_EpqErs_OL

!
! Compute <I|EpqErs|J> with both E_pq and E_rs open (crossing the midlevel)
! as Sum_K <I|Epq|K><K|Ers|J>, where each open loop factorises into an
! upper and a lower coupling coefficient:
!
!   <I|Epq|K> = A(Iu,Ku) x A(Il,Kl),   <K|Ers|J> = A(Ku,Ju) x A(Kl,Jl)
!
! so that the sum over K separates for each intermediate midvertex m_K:
!
!   Sum_mK [ Sum_Ku A(Iu,Ku) A(Ku,Ju) ] x [ Sum_Kl A(Il,Kl) A(Kl,Jl) ]
!
! m_K runs over the (at most two) midvertices reached from m_I by the band
! of E_pq, from which m_J is reached by the band of E_rs. Both operators may
! have either orientation (upper index first or second); the coupling
! tables are accessed as in the open-loop case of MatEl_Epq. The
! symmetries of the bra, the ket and K are kept separate, so this also holds
! for symmetry-changing operators. ISYCJ is the symmetry of the ket J.
!
function MatEl_EpqErs_OO(ISYCJ,I,J,p,q,r,s) result(Value)

use Constants, only: Zero
use Definitions, only: wp, iwp
use Symmetry_Info, only: Mul

implicit none

integer(kind=iwp), intent(in) :: ISYCJ,I,J,p,q,r,s

real(kind=wp) :: Value
real(kind=wp) :: XU, XL

integer(kind=iwp) :: MidLev, nLev
integer(kind=iwp) :: ISYCI
integer(kind=iwp) :: MVI,MVJ,MVK,MV1,MV2
integer(kind=iwp) :: ISYUI,ISYUJ,ISYUK,ISYDI,ISYDJ,ISYDK
integer(kind=iwp) :: IUI,IUJ,IDI,IDJ
integer(kind=iwp) :: U1,L1,U2,L2
integer(kind=iwp) :: iK
logical(kind=iwp) :: Swap1, Swap2

Value=Zero
MidLev=SGS(ME_iState)%MidLev
nLev  =SGS(ME_iState)%nLev

! Both operators must be open.
if (min(p,q,r,s) < 1 .or. max(p,q,r,s) > nLev) return
if ((p > MidLev) .eqv. (q > MidLev)) return
if ((r > MidLev) .eqv. (s > MidLev)) return

! Upper (U) and lower (L) index of each operator, and its orientation.
Swap1 = q > p
U1 = max(p,q); L1 = min(p,q)
Swap2 = s > r
U2 = max(r,s); L2 = min(r,s)

ISYCI = Mul(Mul(SGS(ME_iState)%ISm(p),SGS(ME_iState)%ISm(q)), &
            Mul(Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)),ISYCJ))
if (I < 1 .or. I > CIS(ME_iState)%NCSF(ISYCI)) return
if (J < 1 .or. J > CIS(ME_iState)%NCSF(ISYCJ)) return

call Decode_ME(ISYCI,I,MVI,ISYUI,IUI,IDI)
call Decode_ME(ISYCJ,J,MVJ,ISYUJ,IUJ,IDJ)
ISYDI=Mul(ISYCI,ISYUI)
ISYDJ=Mul(ISYCJ,ISYUJ)

! Symmetries of the halves of K: E_pq takes K to I.
ISYUK=Mul(SGS(ME_iState)%ISm(U1),ISYUI)
ISYDK=Mul(SGS(ME_iState)%ISm(L1),ISYDI)
! and E_rs takes J to K.
if (ISYUJ /= Mul(SGS(ME_iState)%ISm(U2),ISYUK)) return
if (ISYDJ /= Mul(SGS(ME_iState)%ISm(L2),ISYDK)) return

do MV1 = 1, 2
   ! Intermediate midvertex m_K reached from m_I by E_pq.
   MVK = merge(EXS(ME_iState)%MVR(MVI,MV1),EXS(ME_iState)%MVL(MVI,MV1),Swap1)
   if (MVK == 0) cycle
   ! Band of E_rs from m_K that reaches m_J.
   MV2 = 0
   do iK = 1, 2
      if (merge(EXS(ME_iState)%MVR(MVK,iK),EXS(ME_iState)%MVL(MVK,iK),Swap2) == MVJ) MV2 = iK
   end do
   if (MV2 == 0) cycle

   ! Upper half: XU = Sum_Ku A_pq(Iu,Ku) A_rs(Ku,Ju)
   XU = OO_Half(0,U1,U2,ISYUI,ISYUK,ISYUJ,MVI,MVK,MVJ,MV1,MV2,Swap1,Swap2,IUI,IUJ)
   if (XU == Zero) cycle
   ! Lower half: XL = Sum_Kl A_pq(Il,Kl) A_rs(Kl,Jl)
   XL = OO_Half(1,L1,L2,ISYDI,ISYDK,ISYDJ,MVI,MVK,MVJ,MV1,MV2,Swap1,Swap2,IDI,IDJ)

   Value = Value+XU*XL
end do

end function MatEl_EpqErs_OO

!
! One half of MatEl_EpqErs_OO, X = Sum_K A_pq(I,K) A_rs(K,J) over the half
! walks K (upper half iHalf = 0 with the upper levels A1 = max(p,q),
! A2 = max(r,s), lower half iHalf = 1 with the lower levels), for the bra and
! ket half walks IWI, IWJ of symmetries ISYI, ISYJ (ISYK for K), the
! midvertices MVI, MVK, MVJ reached through the bands MV1 (E_pq) and MV2
! (E_rs), and the orientations Swap1, Swap2 of the operators. X depends on
! the graph only and is kept in the cache OO_Key/OO_Val.
!
function OO_Half(iHalf,A1,A2,ISYI,ISYK,ISYJ,MVI,MVK,MVJ,MV1,MV2,Swap1,Swap2,IWI,IWJ) result(X)

use stdalloc, only: mma_allocate
use, intrinsic :: iso_fortran_env, only: int64

integer(kind=iwp), intent(in) :: iHalf, A1, A2, ISYI, ISYK, ISYJ, MVI, MVK, MVJ, MV1, MV2, IWI, IWJ
logical(kind=iwp), intent(in) :: Swap1, Swap2
real(kind=wp) :: X
integer(kind=iwp) :: nLev, INDEO, N1, Off1, N2, Off2, Pack, h

X = Zero
nLev = SGS(ME_iState)%nLev
INDEO = (MV2-1)*nLev+A2
N2 = merge(EXS(ME_iState)%NOCP(INDEO,ISYK,MVK),EXS(ME_iState)%NOCP(INDEO,ISYJ,MVJ),Swap2)
if (N2 <= 0) return
INDEO = (MV1-1)*nLev+A1
N1 = merge(EXS(ME_iState)%NOCP(INDEO,ISYI,MVI),EXS(ME_iState)%NOCP(INDEO,ISYK,MVK),Swap1)
if (N1 <= 0) return
Off2 = merge(EXS(ME_iState)%IOCP((MV2-1)*nLev+A2,ISYK,MVK),EXS(ME_iState)%IOCP((MV2-1)*nLev+A2,ISYJ,MVJ),Swap2)
Off1 = merge(EXS(ME_iState)%IOCP(INDEO,ISYI,MVI),EXS(ME_iState)%IOCP(INDEO,ISYK,MVK),Swap1)

! For small coupling blocks (at most OO_MinN transitions each) Path_Value
! costs about as much as a cache lookup.
if ((max(N1,N2) <= OO_MinN) .or. (nLev > OO_MaxLev) .or. (ME_iState > OO_MaxState)) then
  ! Walk index (Path_Value): only the transitions of E_rs starting at the
  ! ket walk IWJ and of E_pq starting at K.
  X = Path_Value(EXS(ME_iState),Off2,N2,merge(2,1,Swap2),IWJ,Off1,N1,merge(2,1,Swap1),IWI)
  return
end if

Pack = A1+128*(A2+128*(merge(1,0,Swap1)+2*(merge(1,0,Swap2)+2*((MV1-1)+2*((MV2-1)+2*((ISYI-1)+8*(iHalf+2*ME_iState)))))))
! hash in 64-bit arithmetic (no overflow: each term < 2**54)
h = int(mod(int(IWI,kind=int64)*40503_int64+int(IWJ,kind=int64)*65599_int64+int(Pack,kind=int64)*2654435_int64+ &
            int(MVI,kind=int64)*97_int64,int(nOO_Cache,kind=int64)),kind=iwp)+1
if (.not. allocated(OO_Key)) then
  call mma_allocate(OO_Key,4,nOO_Cache,Label='OO_Key')
  call mma_allocate(OO_Val,nOO_Cache,Label='OO_Val')
  OO_Key(:,:) = 0
  OO_Val(:) = Zero
else if ((OO_Key(4,h) == MVI) .and. (OO_Key(1,h) == IWI) .and. (OO_Key(2,h) == IWJ) .and. (OO_Key(3,h) == Pack)) then
  X = OO_Val(h)
  return
end if
X = Path_Value(EXS(ME_iState),Off2,N2,merge(2,1,Swap2),IWJ,Off1,N1,merge(2,1,Swap1),IWI)
OO_Key(:,h) = [IWI,IWJ,Pack,MVI]
OO_Val(h) = X

end function OO_Half

!
! <I|E_pq E_rs|J> for I and J of symmetry LSM, by topology: each operator is
! U (both levels upper), L (both lower) or O (open). Screened with the exact
! quick zero test first, from the occupations OccI and OccJ of I and J.
!
function MatEl_EpqErs(I,J,p,q,r,s,OccI,OccJ) result(X)

  integer(kind=iwp), intent(in) :: I, J, p, q, r, s, OccI(:), OccJ(:)
  real(kind=wp) :: X

  X = Zero
  if (Quick_Zero(OccI,OccJ,p,q,r,s)) return
  X = MatEl_EpqErs_Top(I,J,p,q,r,s)

end function MatEl_EpqErs

!> <I|E_pq E_rs|J> by the topology of the two operators, without the
!> occupation test of MatEl_EpqErs: for callers whose operators are already
!> allowed by the occupations of I and J (MatEl_H_Direct, H_Double).
function MatEl_EpqErs_Top(I,J,p,q,r,s) result(X)

  integer(kind=iwp), intent(in) :: I, J, p, q, r, s
  real(kind=wp) :: X
  character(len=2) :: Top

  X = Zero
  Top = Op_Class(p,q)//Op_Class(r,s)
  select case (Top)
    case ('UU')
      X = MatEl_EpqErs_UU(ME_LSM,I,J,p,q,r,s)
    case ('LL')
      X = MatEl_EpqErs_LL(ME_LSM,I,J,p,q,r,s)
    case ('UL')
      X = MatEl_EpqErs_UL(ME_LSM,I,J,p,q,r,s)
    case ('LU')
      X = MatEl_EpqErs_LU(ME_LSM,I,J,p,q,r,s)
    case ('UO')
      X = MatEl_EpqErs_UO(ME_LSM,I,J,p,q,r,s)
    case ('OU')
      X = MatEl_EpqErs_OU(ME_LSM,I,J,p,q,r,s)
    case ('LO')
      X = MatEl_EpqErs_LO(ME_LSM,I,J,p,q,r,s)
    case ('OL')
      X = MatEl_EpqErs_OL(ME_LSM,I,J,p,q,r,s)
    case ('OO')
      X = MatEl_EpqErs_OO(ME_LSM,I,J,p,q,r,s)
  end select

end function MatEl_EpqErs_Top

!> 'U', 'L' or 'O' for an operator E_pq with both levels in the upper half,
!> both in the lower half, or crossing the midlevel.
function Op_Class(p,q) result(C)
  integer(kind=iwp), intent(in) :: p, q
  character(len=1) :: C
  if ((p > SGS(ME_iState)%MidLev) .and. (q > SGS(ME_iState)%MidLev)) then
    C = 'U'
  else if ((p <= SGS(ME_iState)%MidLev) .and. (q <= SGS(ME_iState)%MidLev)) then
    C = 'L'
  else
    C = 'O'
  end if
end function Op_Class

!
! <I|H|J> as MatEl_H, but only from the operator strings that can connect J
! to I, as given by the occupation difference: with the creations P and the
! annihilations M (multisets, I <- J) the term E_pq E_rs contributes only if
! {p,r} = P + X and {q,s} = M + X, X of size 2-|P|:
!
!   |P| > 2 : zero
!   |P| = 2 : the orderings of P and M, at most 2 loop walks (H_Double)
!   |P| = 1 : a <- b: <I|E_ab|J> times k_ab and the analytic part of the
!             strings with X = {x} (with [E_ax,E_xb] = E_ab):
!               x = a, b : E_ab E_xx, E_xx E_ab = n_x(J), n_x(I) E_ab
!               x other  : (ab|xx) n_x E_ab, and E_ax E_xb + E_xb E_ax =
!                          E_ab (x empty), -E_ab (x doubly occupied),
!                          2 E_ax E_xb - E_ab (x open: one loop walk)
!   |P| = 0 : the diagonal Sum_t k_tt n_t + 1/2 Sum_xy (xx|yy) n_x n_y for
!             I = J, and the exchange terms E_xy E_yx for x /= y: by loop
!             walks only for pairs of open shells, otherwise (I = J only)
!             analytically
!
! with H = Sum_pq k_pq E_pq + 1/2 Sum_pqrs (pq|rs) E_pq E_rs (the form of
! MatEl_H). The cost per element is O(o) MatEl_EpqErs calls for a single
! and O(o**2) (o open shells) for equal occupations, instead of O(n**4).
!
function MatEl_H_Direct(I,J,OccI,OccJ,TU,nTU,TUVX,nTUVX) result(H)

  use Constants, only: Half

  integer(kind=iwp), intent(in) :: I, J, OccI(:), OccJ(:), nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp) :: H
  integer(kind=iwp) :: nL, t, d, k, nP, nM, x, y, a, b
  real(kind=wp) :: E1, W, G
  integer(kind=iwp) :: Pl(2), Mi(2)
  logical(kind=iwp) :: Hit

  H = Zero
  nL = SGS(ME_iState)%nLev

  ! creations (Pl) and annihilations (Mi) taking J to I
  nP = 0
  nM = 0
  do t = 1, nL
    d = OccI(t)-OccJ(t)
    do k = 1, d
      nP = nP+1
      if (nP > 2) return
      Pl(nP) = t
    end do
    do k = 1, -d
      nM = nM+1
      if (nM > 2) return
      Mi(nM) = t
    end do
  end do

  select case (nP)
    case (2)
      H = H_Double(I,J,Pl(1),Pl(2),Mi(1),Mi(2),TUVX,nTUVX)
    case (1)
      a = Pl(1)
      b = Mi(1)
      E1 = MatEl_Epq(SGS(ME_iState),CIS(ME_iState),EXS(ME_iState),a,b,ME_LSM,I,J,.False.)
      ! coefficient W of <I|E_ab|J>: it depends only on the occupations, so it
      ! is kept for the CSF pairs of the same configuration pair (ME_HW).
      if (ME_HW_Valid .and. (ME_HW_ab(1) == a) .and. (ME_HW_ab(2) == b)) then
        Hit = all(ME_HW_Occ(1:nL) == OccJ(1:nL))
      else
        Hit = .false.
      end if
      if (Hit) then
        W = ME_HW
      else
        W = H_Int1(a,b,TU,nTU,TUVX,nTUVX)+Half*H_Int2(a,b,a,a,TUVX,nTUVX)*(OccI(a)+OccJ(a))+ &
            Half*H_Int2(a,b,b,b,TUVX,nTUVX)*(OccI(b)+OccJ(b))
        do x = 1, nL
          if ((x == a) .or. (x == b)) cycle
          W = W+H_Int2(a,b,x,x,TUVX,nTUVX)*OccJ(x)+merge(Half,-Half,OccJ(x) == 0)*H_Int2(a,x,x,b,TUVX,nTUVX)
        end do
        call Keep_Occ(ME_HW_Occ,OccJ,nL)
        ME_HW_ab(:) = [a,b]
        ME_HW = W
        ME_HW_Valid = .true.
      end if
      ! open shells x: the spin coupling needs the loop walk. In a RAS graph
      ! the walk is a product over the graph CSFs K, exact only if the first
      ! operator cannot leave the space (see Top_Comm); otherwise use
      ! E_ax E_xb = E_xb E_ax + E_ab.
      do x = 1, nL
        if ((OccJ(x) /= 1) .or. (x == a) .or. (x == b)) cycle
        G = H_Int2(a,x,x,b,TUVX,nTUVX)
        if ((x <= b) .or. (a >= x)) then
          H = H+G*MatEl_EpqErs_Top(I,J,a,x,x,b)
        else
          H = H+G*(MatEl_EpqErs_Top(I,J,x,b,a,x)+E1)
        end if
      end do
      H = H+W*E1
    case (0)
      ! the part of <I|H|I> that depends only on the occupations (ME_HD):
      ! one-electron and Coulomb terms, and the exchange of the pairs that
      ! are not open-open (see below)
      if (I == J) then
        if (ME_HD_Valid) then
          Hit = all(ME_HD_Occ(1:nL) == OccJ(1:nL))
        else
          Hit = .false.
        end if
        if (Hit) then
          H = ME_HD
        else
          do x = 1, nL
            if (OccJ(x) == 0) cycle
            H = H+H_Int1(x,x,TU,nTU,TUVX,nTUVX)*OccJ(x)
            do y = 1, nL
              H = H+Half*H_Int2(x,x,y,y,TUVX,nTUVX)*OccJ(x)*OccJ(y)
            end do
          end do
          do x = 1, nL
            do y = x+1, nL
              if ((OccJ(x) == 1) .and. (OccJ(y) == 1)) cycle
              H = H+Half*H_Int2(x,y,y,x,TUVX,nTUVX)*(OccJ(x)+OccJ(y)-OccJ(x)*OccJ(y))
            end do
          end do
          call Keep_Occ(ME_HD_Occ,OccJ,nL)
          ME_HD = H
          ME_HD_Valid = .true.
        end if
      end if
      ! Exchange terms. With E_xy E_yx = n_x - 1/2 n_x n_y - 2 S_x.S_y (x /= y)
      ! and S_t = 0 for a closed or empty level, a pair that is not open-open
      ! contributes only to I = J, analytically:
      ! E_xy E_yx + E_yx E_xy = n_x + n_y - n_x n_y.
      ! For an open-open pair E_yx E_xy = E_xy E_yx - (E_xx - E_yy) has the
      ! same matrix element (n_x = n_y), so one loop walk suffices; with
      ! x < y the walk of E_yx E_xy (E_xy first) is also exact in a RAS graph.
      do x = 1, nL
        if (OccJ(x) /= 1) cycle
        do y = x+1, nL
          if (OccJ(y) /= 1) cycle
          H = H+H_Int2(x,y,y,x,TUVX,nTUVX)*MatEl_EpqErs_Top(I,J,y,x,x,y)
        end do
      end do
  end select

end function MatEl_H_Direct

! Double excitation a1,a2 <- b1,b2 (a1 <= a2, b1 <= b2, no a equal to a b):
! the sum of 1/2 (pq|rs) <I|E_pq E_rs|J> over the distinct orderings {p,r}
! of {a1,a2} and {q,s} of {b1,b2}. The two operators of each string have no
! common index, so they commute ([E_pq,E_rs] = delta_qr E_ps - delta_ps E_rq
! = 0) and the strings with the order reversed give the same value: at most
! two loop walks.
function H_Double(I,J,a1,a2,b1,b2,TUVX,nTUVX) result(V)

  use Constants, only: Half

  integer(kind=iwp), intent(in) :: I, J, a1, a2, b1, b2, nTUVX
  real(kind=wp), intent(in) :: TUVX(nTUVX)
  real(kind=wp) :: V

  if ((a1 == a2) .and. (b1 == b2)) then
    ! the single string E_ab E_ab
    V = Half*H_Int2(a1,b1,a1,b1,TUVX,nTUVX)*MatEl_EpqErs_Top(I,J,a1,b1,a1,b1)
  else if ((a1 == a2) .or. (b1 == b2)) then
    ! E_a1b1 E_a2b2 and its reverse
    V = H_Int2(a1,b1,a2,b2,TUVX,nTUVX)*Top_Comm(I,J,a1,b1,a2,b2)
  else
    ! E_a1b1 E_a2b2, E_a1b2 E_a2b1 and their reverses
    V = H_Int2(a1,b1,a2,b2,TUVX,nTUVX)*Top_Comm(I,J,a1,b1,a2,b2)+ &
        H_Int2(a1,b2,a2,b1,TUVX,nTUVX)*Top_Comm(I,J,a1,b2,a2,b1)
  end if

end function H_Double

! <I|E_pq E_rs|J> for commuting E_pq and E_rs, as the loop walk of the order
! that is exact in a RAS graph. The walk is a product over the graph CSFs K,
! Sum_K <I|E_pq|K><K|E_rs|J>; with the levels ordered RAS1 < RAS2 < RAS3 an
! intermediate E_rs|J> can leave the space only if r > s, and <I|E_pq| only
! if p < q, so E_pq E_rs is exact unless r > s and p < q, and then E_rs E_pq
! is. In CAS both orders are exact.
function Top_Comm(I,J,p,q,r,s) result(X)

  integer(kind=iwp), intent(in) :: I, J, p, q, r, s
  real(kind=wp) :: X

  if ((r <= s) .or. (p >= q)) then
    X = MatEl_EpqErs_Top(I,J,p,q,r,s)
  else
    X = MatEl_EpqErs_Top(I,J,r,s,p,q)
  end if

end function Top_Comm

! (pq|rs) from the triangular packing of TUVX.
function H_Int2(p,q,r,s,TUVX,nTUVX) result(V)

  integer(kind=iwp), intent(in) :: p, q, r, s, nTUVX
  real(kind=wp), intent(in) :: TUVX(nTUVX)
  real(kind=wp) :: V

  V = TUVX(Tri_SG(Tri_SG(p,q),Tri_SG(r,s)))

end function H_Int2

! k_pq = h_pq - 1/2 Sum_k (pk|kq).
function H_Int1(p,q,TU,nTU,TUVX,nTUVX) result(V)

  use Constants, only: Half

  integer(kind=iwp), intent(in) :: p, q, nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp) :: V
  integer(kind=iwp) :: kO

  V = Zero
  do kO = 1, SGS(ME_iState)%nLev
    V = V+H_Int2(p,kO,kO,q,TUVX,nTUVX)
  end do
  V = TU(Tri_SG(p,q))-Half*V

end function H_Int1

function MatEl_EpqErs_LU( &
           ISYCI, &
           I,J, &
           p,q,r,s ) result(Value)

use Constants,   only: Zero
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J
integer(kind=iwp), intent(in) :: p,q,r,s

integer(kind=iwp) :: MVI
integer(kind=iwp) :: MVJ
integer(kind=iwp) :: ISYUPI
integer(kind=iwp) :: ISYUPJ
integer(kind=iwp) :: IUPI
integer(kind=iwp) :: IUPJ
integer(kind=iwp) :: IDWNI
integer(kind=iwp) :: IDWNJ

real(kind=wp) :: XL
real(kind=wp) :: XU
real(kind=wp) :: Value

Value = Zero

!
! Left operator is lower, right operator is upper.
!
if (p < 1 .or. q < 1) return
if (p > SGS(ME_iState)%MidLev) return
if (q > SGS(ME_iState)%MidLev) return

if (r <= SGS(ME_iState)%MidLev) return
if (s <= SGS(ME_iState)%MidLev) return
if (r >  SGS(ME_iState)%nLev) return
if (s >  SGS(ME_iState)%nLev) return

call Decode_ME( &
     ISYCI,I, &
     MVI,ISYUPI,IUPI,IDWNI )

call Decode_ME( &
     ISYCI,J, &
     MVJ,ISYUPJ,IUPJ,IDWNJ )

!
! The upper operator changes the upper-walk symmetry by Gamma(r)Gamma(s)
! (and the lower one the lower-walk symmetry by Gamma(p)Gamma(q), the
! same for a totally symmetric product). A coupling block is indexed by
! the symmetry of its ICoup(1) walk: the source (J) in canonical
! orientation, the destination (I) if reversed.
!
if (MVI    /= MVJ)    return
if (ISYUPI /= Mul(Mul(SGS(ME_iState)%ISm(r),SGS(ME_iState)%ISm(s)),ISYUPJ)) return

XL = MatEl_Lower_Half( &
        ISYCI,merge(ISYUPI,ISYUPJ,p < q),MVI, &
        IDWNI,IDWNJ, &
        p,q )

if (abs(XL) <= 1.0e-15_wp) return

XU = MatEl_Upper_Half( &
        merge(ISYUPI,ISYUPJ,r < s),MVI, &
        IUPI,IUPJ, &
        r,s )

Value = XL*XU

end function MatEl_EpqErs_LU

function MatEl_EpqErs_UL( &
           ISYCI, &
           I,J, &
           p,q,r,s ) result(Value)

!
! <I|E_pq E_rs|J> with E_pq in the upper and E_rs in the lower half.
!
! For real CI vectors <I|E_pq E_rs|J> = <J|(E_pq E_rs)^+|I> = <J|E_sr E_qp|I>,
! and E_sr E_qp is a lower x upper product, so this is MatEl_EpqErs_LU with
! bra and ket exchanged and the indices (p,q,r,s) -> (s,r,q,p). (Since the
! two operators share no level they also commute, E_pq E_rs = E_rs E_pq.)
!
! ISYCI is the symmetry of both bra and ket, as in MatEl_EpqErs_LU: the
! present implementation only treats operators that individually preserve
! spatial symmetry. Once that is generalised, the symmetry passed here must
! be that of I, the ket of the exchanged matrix element.
!
use Definitions, only: wp, iwp
implicit none

integer(kind=iwp), intent(in) :: ISYCI
integer(kind=iwp), intent(in) :: I
integer(kind=iwp), intent(in) :: J
integer(kind=iwp), intent(in) :: p
integer(kind=iwp), intent(in) :: q
integer(kind=iwp), intent(in) :: r
integer(kind=iwp), intent(in) :: s

real(kind=wp) :: Value

Value = MatEl_EpqErs_LU(ISYCI,J,I,s,r,q,p)

end function MatEl_EpqErs_UL

function MatEl_Epq(SGS,CIS,EXS,IP,IQ,ISYJ,I_CSF,J_CSF,Verbose) result(EIJ)
  use Symmetry_Info, only: Mul
  use Constants,     only: Zero
  use Definitions,   only: wp, iwp
#include "intent.fh"
  type(SGStruct), intent(in) :: SGS
  type(CIStruct), intent(in) :: CIS
  type(EXStruct), intent(in) :: EXS
  integer(kind=iwp), intent(in) :: IP, IQ, ISYJ, I_CSF, J_CSF
  logical(kind=iwp) :: Verbose
  real(kind=wp) :: EIJ

  integer(kind=iwp) :: ISYP, ISYQ, ISYPQ, ISYI
  integer(kind=iwp) :: MVI, MVJ, ISYUI, ISYUJ
  integer(kind=iwp) :: IUP, IDN, JUP, JDN
  logical(kind=iwp) :: Swap, Lower, Upper
  integer(kind=iwp) :: INDEO, NCP, LICP
  integer(kind=iwp) :: INDEO_P, INDEO_Q
  integer(kind=iwp) :: MV, MVX, NCP1, NCP2, LICP1, LICP2, iBand
  integer(kind=iwp) :: ISYDI, ISYUC, ISYDC, JP, JQ, ISYDJ
  integer(kind=iwp) :: IOUW, IOLW, IPSHFT, IPPOW, IC, ICS
  real(kind=wp) :: Upper_Factor, Lower_Factor

  EIJ  = Zero
  ISYP = SGS%ISm(IP)      ! Gamma(p)
  ISYQ = SGS%ISm(IQ)      ! Gamma(q)
  ISYPQ= Mul(ISYP,ISYQ)  ! Gamma(pq)
  ISYI = Mul(ISYPQ,ISYJ)  ! Gamma(I)=Gamma(pq) x Gamma(J)
  If (Verbose) Then
     Write (u6,*)
     Write (u6,*) 'Enter: MatEl_epq'
     Write (u6,*) 'Testing: I,p,q,J=',i_CSF,ip,iq,j_CSF
     Write (u6,*) '        ISYI, ISYJ  =', ISYI, ISYJ
     Write (u6,*) '        I_CSF,j_CSF=', I_CSF,j_CSF
     Write (u6,*) '        nCSFi,nCSFj=',CIS%NCSF(ISYI),CIS%NCSF(ISYJ)
  End If

  if (I_CSF < 1 .or. I_CSF > CIS%NCSF(ISYI)) return
  if (J_CSF < 1 .or. J_CSF > CIS%NCSF(ISYJ)) return

  call Decode_CSF(CIS,ISYI,I_CSF,MVI,ISYUI,IUP,IDN)
  call Decode_CSF(CIS,ISYJ,J_CSF,MVJ,ISYUJ,JUP,JDN)
  ISYDI = Mul(ISYI,ISYUI)   ! Gamma(Il) = Gamma(I) x Gamma(Iu)
  ISYDJ = Mul(ISYJ,ISYUJ)   ! Gamma(Jl) = Gamma(J) x Gamma(Ju)

  If (Verbose) Then
     Write (u6,*) '        ISYP, ISYQ  =', ISYP, ISYQ
     Write (u6,*) '        MVI, MVJ  =', MVI, MVJ
     Write (u6,*) '        IUP, JUP  =', IUP, JUP
     Write (u6,*) '        IDN, JDN  =', IDN, JDN
     Write (u6,*) '        ISYI, ISYJ=', ISYI, ISYJ
  End If

  if (IP == IQ) then
    If (Verbose) Write (u6,*) 'Diagonal case'
    ! E_pp is diagonal in this representation.

    if (ISYI /= ISYJ .or. I_CSF /= J_CSF) return

    if (IP > SGS%MidLev) then
      IOUW = CIS%IOW(1,ISYUI,MVI)
      IPSHFT = 2*(IP-1-SGS%MidLev)
      IC = CIS%ICase(1+IOUW-CIS%nIpWlk+IPSHFT/(2*nPack) &
                     +IUP*CIS%nIpWlk)
    else
      IOLW = CIS%IOW(2,ISYDI,MVI)
      IPSHFT = 2*(IP-1)
      IC = CIS%ICase(1+IOLW-CIS%nIpWlk+IPSHFT/(2*nPack) &
                     +IDN*CIS%nIpWlk)
    end if
    IPSHFT = mod(IPSHFT,2*nPack)
    IPPOW = 2**IPSHFT
    ICS = mod(IC/IPPOW,4)
    if (ICS /= 0) EIJ = real((1+ICS)/2,kind=wp)
    If (Verbose) Then
    Upper_Factor = 1.0_wp
    Lower_Factor = 1.0_wp
    if (IP > SGS%MidLev) then
      Upper_Factor = EIJ
    else
      Lower_Factor = EIJ
    end if
    write(u6,'(/,A)') 'MatEl_Epq: CLOSED LOOP (p=q)'
    write(u6,'(A,5I8)') '  I,J,p,q,ICS       =', &
                        I_CSF,J_CSF,IP,IQ,ICS
    write(u6,'(A,ES24.14)') '  upper factor      =',Upper_Factor
    write(u6,'(A,ES24.14)') '  lower factor      =',Lower_Factor
    write(u6,'(A,ES24.14)') '  contribution      =',EIJ
    Write (u6,*) 'Exit: MatEl_epq'
    Write (u6,*)
    End If
    return
  end if

  Swap = IQ > IP
  If (Verbose) Write (u6,*) 'SWAP=',SWAP
  Lower = merge(IQ <= SGS%MidLev,IP <= SGS%MidLev,Swap)
  Upper = merge(SGS%MidLev < IP,SGS%MidLev < IQ,Swap)

  if (Lower) then
    If (Verbose) Write (u6,*) 'Non-diagonal lower closed-loop case'
    ! Both orbital levels are in the lower half.  The upper walk is unchanged.
    if (MVI /= MVJ .or. ISYUI /= ISYUJ .or. IUP /= JUP) return
    INDEO = merge(2*SGS%nLev+(IQ*(IQ-1))/2+IP, &
                  2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
    If (Verbose) Write (u6,*) ' INDEO_PQ =',INDEO

    ! Gamma(pq) x Gamma(Il)
    ISYDC  = Mul(ISYPQ,ISYDI)
    If (ISYDC /= ISYDJ) Call Abend()

    NCP = merge(EXS%NOCP(INDEO,ISYDI,MVI), &
                EXS%NOCP(INDEO,ISYDJ,MVI),Swap)
    LICP = merge(EXS%IOCP(INDEO,ISYDI,MVI), &
                 EXS%IOCP(INDEO,ISYDJ,MVI),Swap)
    If (Verbose) Write (u6,*) ' NCP, LICP=',NCP,LICP
    Upper_Factor = 1.0_wp
    ! the transition JDN -> IDN (walk index)
    Lower_Factor = Coup_Value(EXS,LICP,NCP,merge(2,1,Swap),JDN,IDN)
    EIJ = Lower_Factor*Upper_Factor
    If (Verbose) Then
    write(u6,'(/,A)') 'MatEl_Epq: CLOSED LOOP (p/=q)'
    write(u6,'(A)') '  case              = lower-half excitation'
    write(u6,'(A,4I8)') '  I,J,p,q           =',I_CSF,J_CSF,IP,IQ
    Write(u6,'(A,L1)') 'Swap=',Swap
    write(u6,'(A,ES24.14)') '  upper factor      =',Upper_Factor
    write(u6,'(A,ES24.14)') '  lower factor      =',Lower_Factor
    write(u6,'(A,ES24.14)') '  EIJ               =',EIJ
    Write (u6,*) 'Exit: MatEl_epq'
    Write (u6,*)
    End If

  else if (Upper) then
    If (Verbose) Write (u6,*) 'Non-diagonal upper closed-loop case'
    ! Both orbital levels are in the upper half.  The lower walk is unchanged.
    ! Walk indices are local to their symmetry block, so the symmetry of the
    ! lower walk must be compared as well as its index.
    if (MVI /= MVJ .or. ISYDI /= ISYDJ .or. IDN /= JDN) return
    INDEO = merge(2*SGS%nLev+(IQ*(IQ-1))/2+IP, &
                  2*SGS%nLev+(IP*(IP-1))/2+IQ,Swap)
    If (Verbose) Write (u6,*) ' INDEO_PQ =',INDEO
    NCP = merge(EXS%NOCP(INDEO,ISYUI,MVI), &
                EXS%NOCP(INDEO,ISYUJ,MVI),Swap)
    LICP = merge(EXS%IOCP(INDEO,ISYUI,MVI), &
                 EXS%IOCP(INDEO,ISYUJ,MVI),Swap)
    If (Verbose) Write (u6,*) ' NCP, LICP=',NCP,LICP
    Lower_Factor = 1.0_wp
    ! the transition JUP -> IUP (walk index)
    Upper_Factor = Coup_Value(EXS,LICP,NCP,merge(2,1,Swap),JUP,IUP)
    EIJ = Lower_Factor*Upper_Factor
    If (Verbose) Then
    write(u6,'(/,A)') 'MatEl_Epq: CLOSED LOOP (p/=q)'
    write(u6,'(A)') '  case              = upper-half excitation'
    write(u6,'(A,4I8)') '  I,J,p,q           =',I_CSF,J_CSF,IP,IQ
    Write(u6,'(A,L1)') 'Swap=',Swap
    write(u6,'(A,ES24.14)') '  upper factor      =',Upper_Factor
    write(u6,'(A,ES24.14)') '  lower factor      =',Lower_Factor
    write(u6,'(A,ES24.14)') '  EIJ               =',EIJ
    Write (u6,*) 'Exit: MatEl_epq'
    Write (u6,*)
    End If

  else
    If (Verbose) Write (u6,*) 'Non-diagonal open-loop case'
    ! The excitation crosses MIDLEV.  Contract the two coupling tables
    ! directly, retaining only paths from J to I.

    ! Evaluate <I(IUP,IDN,MVI)|Epq|J(JUP,JDN,MVJ)> =
    !          A_p(IUP,MVI,JUP,MVJ) x A_q(IDN,MVI,JDN,MVJ)

    JP=merge(IQ,IP,Swap)
    JQ=merge(IP,IQ,Swap)

    ! Find if the midvertex MVI connects with MVJ for delta b = -1 or +1
    iBand=0
    do MV=1,2
      ! Find the Midvertex to the left of MVI
      MVX=merge(EXS%MVR(MVI,MV),EXS%MVL(MVI,MV),Swap)
      If (Verbose) Then
         Write (u6,*) 'MVX=',MVX
         Write (u6,*) 'MVX/=MVJ',MVX /= MVJ
         Write (u6,*) 'Skip= ',MVX == 0 .or. MVX /= MVJ
      End If
      ! branch out if there are no such vertex or if the index don't match
      if (MVX == 0 .or. MVX /= MVJ) cycle
      iBand=MV
    end do

    If (iBand==0) Return
    MV=iBand
    If (Verbose) Write (u6,*) 'MV=',MV

    !Gamma(Ju) = Gamma(p) x Gamma(Iu)
    ISYUC=merge(Mul(ISYQ,ISYUI),Mul(ISYP,ISYUI),Swap)
    if (ISYUC /= ISYUJ) Return

    ! Process A^p(Iu(MVI),Ju(MVJ))
    ! Accessed as (INDEO_P,ISYUJ,MVJ)
    INDEO_P=(MV-1)*SGS%nLev+JP
    If (Verbose) Write (u6,*) 'INDEO_P=',INDEO_P
    NCP1= merge(EXS%NOCP(INDEO_P,ISYUI,MVI), &
                EXS%NOCP(INDEO_P,ISYUJ,MVJ),Swap)
    LICP1=merge(EXS%IOCP(INDEO_P,ISYUI,MVI), &
                EXS%IOCP(INDEO_P,ISYUJ,MVJ),Swap)
    ! the transition JUP -> IUP (walk index)
    Upper_Factor=Coup_Value(EXS,LICP1,NCP1,merge(2,1,Swap),JUP,IUP)

    If (Upper_Factor==Zero) Return

    !Gamma(Jl) = Gamma(q) x Gamma(Il)
    ISYDC=merge(Mul(ISYP,ISYDI),Mul(ISYQ,ISYDI),Swap)
    if (iSYDC /= ISYDJ) Return

    ! Process A^q(Il(MVI),Jl(MVJ))
    ! Accessed as (INDEO_Q,ISYDJ,MVJ)
    INDEO_Q=(MV-1)*SGS%nLev+JQ
    If (Verbose) Write (u6,*) 'INDEO_Q=',INDEO_Q
    NCP2= merge(EXS%NOCP(INDEO_Q,ISYDI,MVI), &
                EXS%NOCP(INDEO_Q,ISYDJ,MVJ),Swap)
    LICP2=merge(EXS%IOCP(INDEO_Q,ISYDI,MVI), &
                EXS%IOCP(INDEO_Q,ISYDJ,MVJ),Swap)
    ! the transition JDN -> IDN (walk index)
    Lower_Factor=Coup_Value(EXS,LICP2,NCP2,merge(2,1,Swap),JDN,IDN)

    If (Lower_Factor==Zero) Return

    EIJ= Upper_Factor*Lower_Factor

    If (Verbose) Then
       write(u6,'(/,A)') 'MatEl_Epq: OPEN LOOP (p/=q)'
       write(u6,'(A)') '  case              = crossing excitation'
       write(u6,'(A,4I8)') '  I,J,p,q           =',I_CSF,J_CSF,IP,IQ
       Write(u6,'(A,L1)') 'Swap=',Swap
       write(u6,'(A,2I8)') '  MVI, MVX          =',MVI,MVX
       write(u6,'(A,I8)') '  crossing branch   =',MV
       write(u6,'(A,ES24.14)') '  upper factor      =',Upper_Factor
       write(u6,'(A,ES24.14)') '  lower factor      =',Lower_Factor
       write(u6,'(A,ES24.14)') '  EIJ               =',EIJ
       Write (u6,*) 'Exit: MatEl_epq'
       Write (u6,*)
    EndIf

  end if


end function MatEl_Epq

!> <I|H|J> between the GUGA CSFs I and J with the occupations OccI and OccJ
!> (SG_CSF_Occ), H = Sum_pq h_pq E_pq + 1/2 Sum_pqrs (pq|rs) (E_pq E_rs -
!> delta_qr E_ps) with the triangularly packed TU and TUVX, without core
!> energy (MatEl_H_Direct, after SG_MatEl_Init); I and J of the same
!> symmetry.
function SG_MatEl_H(I,J,OccI,OccJ,TU,nTU,TUVX,nTUVX) result(H)
  integer(kind=iwp), intent(in) :: I, J, OccI(:), OccJ(:), nTU, nTUVX
  real(kind=wp), intent(in) :: TU(nTU), TUVX(nTUVX)
  real(kind=wp) :: H
  H = MatEl_H_Direct(I,J,OccI,OccJ,TU,nTU,TUVX,nTUVX)
end function SG_MatEl_H

end module SGUGA_MatEl
