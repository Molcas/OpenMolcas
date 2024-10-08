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

subroutine Cho_VecBuf_Init(Frac,lVec)
!
! Purpose: allocate and initialize vector buffer.
!          RUN_MODE=RUN_INTERNAL: buffer used during decomposition.
!          RUN_MODE=RUN_EXTERNAL: buffer used after decomposition,
!                                 i.e. vectors are available on disk.

use Cholesky, only: ip_CHVBFI_SYM, l_CHVBFI_SYM, nSym, RUN_INTERNAL, RUN_EXTERNAL, RUN_MODE
#ifdef _DEBUGPRINT_
use Cholesky, only: LuPri
use stdalloc, only: mma_maxDBLE
#endif
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: Frac
integer(kind=iwp), intent(in) :: lVec(*)
#ifdef _CHO_DEBUGPRINT_
#define _DEBUGPRINT_
#endif
#ifdef _DEBUGPRINT_
#define _DBG_ .true.
#else
#define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: LocDbg = _DBG_
character(len=*), parameter :: SecNam = 'Cho_VecBuf_Init'
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: l_Max, MF
real(kind=wp) :: xMF
character(len=2) :: Unt
#endif

#ifdef _DEBUGPRINT_
call mma_maxDBLE(l_max)
write(Lupri,*) '>>>>> Enter ',SecNam,' <<<<<'
write(Lupri,*) 'Memory fraction requested for buffer: ',Frac
call Cho_Word2Byte(l_Max,8,xMF,Unt)
write(Lupri,*) 'Memory available: ',l_Max,' = ',xMF,Unt
MF = int(Frac*dble(l_Max))
call Cho_Word2Byte(MF,8,xMF,Unt)
write(Lupri,*) 'Memory fraction : ',MF,' = ',xMF,Unt
call XFlush(Lupri)
#endif

l_ChVBfI_Sym(1:nSym) = 0
ip_ChVBfI_Sym(1:nSym) = 0

if (RUN_MODE == RUN_INTERNAL) then
  call Cho_VecBuf_Init_I(Frac,lVec,LocDbg)
else if (RUN_MODE == RUN_EXTERNAL) then
  call Cho_VecBuf_Init_X(Frac,LocDbg)
else
  call Cho_Quit('RUN_MODE error in '//SecNam,103)
end if

#ifdef _DEBUGPRINT_
write(Lupri,*) '>>>>> Exit  ',SecNam,' <<<<<'
call XFlush(Lupri)
#endif

end subroutine Cho_VecBuf_Init
