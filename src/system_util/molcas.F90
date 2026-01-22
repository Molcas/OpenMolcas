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

module Molcas

use Definitions, only: iwp

implicit none
private

#ifdef _DEMO_
# define _MAXBFN_ 100
# define _MXATOM_ 6
# define _MXROOT_ 2
# define _MXNEMOATOM_ 6
# define _MXDBSC_ 10
#else
# define _MAXBFN_ 10000
# define _MXATOM_ 5000
# define _MXROOT_ 600
# define _MXNEMOATOM_ 200
# define _MXDBSC_ 1000
#endif
! We normally expect the number of auxiliary functions to be
! in the range of 3-5 times the normal basis set.
#define _MAXBFN_AUX_ 7*_MAXBFN_
#define _LENIN_ 6

integer(kind=iwp), parameter :: lCache = 64*1024/8, LenIn = _LENIN_, LenIn1 = _LENIN_+1, LenIn2 = _LENIN_+2, LenIn3 = _LENIN_+3, &
                                LenIn4 = _LENIN_+4, LenIn5 = _LENIN_+5, LenIn6 = _LENIN_+6, LenIn8 = _LENIN_+8, MaxBfn = _MAXBFN_, &
                                MaxBfn_Aux = _MAXBFN_AUX_, MxAct = 100, MxAO = _MAXBFN_+_MAXBFN_AUX_, MxAtom = _MXATOM_, &
                                MxBas = _MAXBFN_, Mxdbsc = _MXDBSC_, MxGAS = 16, MxIna = _MAXBFN_, MxNemoAtom = _MXNEMOATOM_, &
                                MxOrb = _MAXBFN_, MxRoot = _MXROOT_, MxSym = 8

public :: lCache, LenIn, LenIn1, LenIn2, LenIn3, LenIn4, LenIn5, LenIn6, LenIn8, MaxBfn, MaxBfn_Aux, MxAct, MxAO, MxAtom, MxBas, &
          Mxdbsc, MxGAS, MxIna, MxNemoAtom, MxOrb, MxRoot, mxSym

end module Molcas
