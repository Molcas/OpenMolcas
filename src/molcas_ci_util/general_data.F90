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

! This is just an encapsulation of the common blocks in
! src/Include/rasdim.fh
! src/Include/general.fh
! into a data module

module general_data

use Definitions, only: wp, iwp

implicit none

#include "rasdim.fh"
Private :: MaxBfn,MaxBfn_Aux,MxAO,mxAtom,mxroot,mxNemoAtom,Mxdbsc,lCache,mxact,mxina,mxbas,mxOrb,mxSym,mxGAS, &
           LENIN,LENIN1,LENIN2,LENIN3,LENIN4,LENIN5,LENIN6,LENIN8
Private :: mxRef,mxIter,mxCiIt,mxSxIt,mxTit

REAL(kind=wp)     SXDAMP
INTEGER(kind=iwp), PARAMETER:: MAXALTER=16
INTEGER(kind=iwp) NSYM,                                             &
                  NBAS(mxSym),NFRO(mxSym),NORB(mxSym),NDEL(mxSym),  &
                  NISH(mxSym),NASH(mxSym),NSSH(mxSym),              &
                  NRS1(mxSym),NRS2(mxSym),NRS3(mxSym),nSkipX(MxSym),&
                  NTOT,NTOT1,NTOT2,NFROT,NDELT,NRS1T,NRS2T,NRS3T,   &
                  NACTEL,ISPIN,STSYM,NCONF,NHOLE1,NELEC3,NSEL,      &
                  NTOTSP,INVEC,NALTER,MALTER(MAXALTER,3),           &
                  NCRVEC,NCRPROJ
!
!     common logical unit numbers
!
!     LUStartOrb : MO-coefficients and occupation numbers
!                   (formatted ASCI file, input)
!     JOBIPH : MO-coefficients and occupation numbers etc.
!                   (binary, output)
!     JOBOLD : MO-coefficients and occupation numbers etc.
!                   (binary, input)
!     LUONEL : one-electron integrals in AO basis
!                   (binary, input)
!     LUINTA : two-electron integrals in AO basis
!                   (binary, input)
!     LUINTM : two-electron integrals in MO basis
!                   (binary, temporary)
!     LURLX  : geometries, gradients, hessians etc.
!                   (binary, input/output)
!     LUCOM  : seward's info block, reaction fields etc.
!                   (binary, input/output)
!     LUEXT  : summary of results
!                   (formatted ASCI file, input)
!     LUQUNE : orbital gradients
!                   (binary, temporary)
!     LUDAVID: Intermediate results of the diagonalization
!                   (binary, temporary)
!
Character(LEN=256)  StartOrbFile
INTEGER(kind=iwp)LUStartOrb,JOBIPH,JOBOLD,LUDAVID,                 &
                 LUONEL,LUINTA,LUINTM,                             &
                 LURLX,LUCOM,LUEXT,LUQUNE,ITERFILE

Logical Lowdin_ON

integer(kind=iwp), allocatable :: CleanMask(:)
real(kind=wp), allocatable :: CRPROJ(:), CRVEC(:)

end module general_data
