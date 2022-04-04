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

module AMFI_global

use Definitions, only: wp, iwp

implicit none
private

!bs Lmax            : max. angular momentum of basis functions
!bs                   DO NOT INCREASE TO MORE THAN SIX   !!!!!
!bs                   if you do, you will have to edit the ixyzpow array by hand...
!bs Lmax_occ        : highest L-value for occupied orbitals
!bs MxprimL         : max. of primitives per angular momentum
!bs MxcontL         : max. of contracted functions per angular momentum
!bs MxCart          : max. number of contracted functions in the atom
!bs ndfmx           : dimension of precomputed double factorials
!bs normovlp        : overlap of normalized functions
!bs contrarray      : big array with enough space for all modified contraction coefficients
!bs                   for each l-value there are five blocks of size MxprimL*MxcontL
!bs                   the original contraction coefficients (for normalized functions)
!bs                   and four modified blocks depending on different kinematic factors and included exponents
!bs exponents       : the exponents
!bs nprimt,ncontrac : the numbers of primitive and contracted functions for each l-value
!bs Lfirst(i)       : gives the first L-value, for which radial integrals are calculated
!bs                   for type i and l1,l2,l3,l4 - Integral block.
!bs Llast(i)        : gives the last L-value
!bs Lblocks(i)      : gives the number of L-values
!bs Lstarter(i)     : gives the adress of each integral block on cont4
!bs AOcoeffs        : express AOs in contracted functions
!bs                   first index: number of contracted function
!bs                   second index: number of AO
!bs                   third index: L-value
!bs occup           : occupation numbers
!bs                   first index: number of AO
!bs                   second index: L-value
!bs numbofsym       : number of symmetries
!bs ipow2ired       : gives IR by checking powers
!bs incrLM          : shift of orbitalnumber in IR for L,M
!bs shiftIRED       : shift to get to absolute number from relative number in IR
!bs iredLM          : IR for L and M
!bs shiftIRIR       : shift for (IR1,IR2)-block (IR1<=IR2)
!bs Loffunction     : gives L value of cartesian function
!bs Moffunction     : gives M value of cartesian function
!bs Iredoffunctnew  : give IRED of cartesian function incl. add. functions
!bs itotalperIR     : total number of functions per IR
!bs df, dffrac      : some double factorials and their fractions, initialized by inidf
!bs ipowxyz         : array that includes information about
!bs                   odd powers of x y z in the real harmonics
!bs                   this is used to check whether integrals in the
!bs                   cartesian representation appear

integer(kind=iwp), parameter :: Lmax = 6, Lmax_occ = 3, MxCart = 300, MxcontL = 40, MxprimL = 40, ndfmx = 4*(Lmax+1)

integer(kind=iwp) :: icore(0:Lmax), ikeeplist(Mxcart), ikeeporb, incrLM(-Lmax:Lmax,0:Lmax), ipow2ired(0:1,0:1,0:1), &
                     ipowxyz(3,-Lmax:Lmax,0:Lmax), iredLM(-Lmax:Lmax,0:Lmax), Iredoffunctnew(Mxcart), itotalperIR(8), Lblocks(4), &
                     Lfirst(4), Llast(4), Loffunction(Mxcart), Lstarter(4), Lvalues(4), Moffunction(Mxcart), nblock, &
                     ncontrac(0:Lmax), ncontrac_keep, noccorb(0:Lmax), nprimit(0:Lmax), nrtofiperIR(8), numbofsym, shiftIRED(8), &
                     shiftIRIR(8*(8+1)/2)
real(kind=wp) :: AOcoeffs(MxcontL,MxcontL,0:Lmax), charge, cntscrtch(MxPrimL,MxcontL,0:Lmax), &
                 contrarray(MxcontL*MxprimL,0:4,0:Lmax), df(0:ndfmx), dffrac(0:ndfmx,0:ndfmx), Exp_Finite, &
                 exponents(MxprimL,0:Lmax), normovlp(MxprimL,MxprimL,0:Lmax), occup(MxcontL,0:Lmax), &
                 OVLPinv(MxprimL,MxprimL,0:Lmax), rootOVLP(MxprimL,MxprimL,0:Lmax), rootOVLPinv(MxprimL,MxprimL,0:Lmax)

public :: AOcoeffs, charge, cntscrtch, contrarray, df, dffrac, Exp_finite, exponents, icore, ikeeplist, ikeeporb, incrLM, &
          ipow2ired, ipowxyz, iredLM, iredoffunctnew, itotalperIR, Lblocks, Lfirst, Llast, Lmax, Lmax_occ, Loffunction, Lstarter, &
          Lvalues, Moffunction, MxCart, MxcontL, MxprimL, nblock, ncontrac, ncontrac_keep, noccorb, normovlp, nprimit, &
          nrtofiperIR, numbofsym, occup, OVLPinv, rootOVLP, rootOVLPinv, shiftIRED, shiftIRIR

end module AMFI_global
