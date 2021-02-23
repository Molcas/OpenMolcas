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
!----------------------------------------------------------------------*
!                                                                      *
!     Define the vocabulary (ComTab) :                                 *
!                                                                      *
!     Command        1st Subcom.  2nd Subcom.  Parameters              *
!     TITL              --           --        Text (10*72 Characters) *
!     FFPT             DIPO         COMP       X=xx,Y=yy,Z=zz          *
!                      EFLD         COMP       X=xx,Y=yy,Z=zz          *
!                                   ORIG       N=i,X=xx,Y=yy,Z=zz      *
!                      QUDR         COMP       XX=xx,YY=yy,ZZ=zz       *
!                                              XY=xy,XZ=xz,YZ=yz,RR=rr *
!                                   ORIG       N=i,X=xx,Y=yy,Z=zz      *
!                      OCTU         COMP       XXX=xxx,XYY=xyy,XZZ=xzz *
!                                              XXY=xxy,YYY=yyy,YZZ=yzz *
!                                              XXZ=xxz,YYZ=yyz,ZZZ=zzz *
!                                              XYZ=xyz                 *
!                                   ORIG       N=i,X=xx,Y=yy,Z=zz      *
!                      EFGR         COMP       XX=xx,YY=yy,ZZ=zz       *
!                                              XY=xy,XZ=xz,YZ=yz,RR=rr *
!                                   ORIG       N=i,X=xx,Y=yy,Z=zz      *
!                      RELA          --        W=ww                    *
!     END               --           --        --                      *
!                                                                      *
!     Define also the command read control tables:                     *
!     ComCtl : Count the number of entries for each hierarchy          *
!              level of the vocabulary                                 *
!     ComStk : flag for each command which has been entered            *
!     ComVal : parameter values read in                                *
!                                                                      *
!----------------------------------------------------------------------*

module FFPT_global

use Definitions, only: wp, iwp

implicit none
private

#include "Molcas.fh"

integer(kind=iwp), parameter :: nCom = 5, MxSub1 = 6, MxSub2 = 2, MxParm = 10
character(len=4) :: ComTab(nCom,0:MxSub1,0:MxSub2,0:MxParm)
integer(kind=iwp) :: ComCtl(nCom,0:MxSub1,0:MxSub2)
logical(kind=iwp) :: ComStk(nCom,0:MxSub1,0:MxSub2,0:MxParm)
real(kind=wp) :: ComVal(nCom,0:MxSub1,0:MxSub2,0:MxParm)

!----------------------------------------------------------------------*
!                                                                      *
!     Allocate space to store general perturbation labels, components  *
!     and weights                                                      *
!     gLblN : general label name                                       *
!     gLblW : general label weight                                     *
!     gLblC : general label component                                  *
!                                                                      *
!----------------------------------------------------------------------*

integer(kind=iwp), parameter :: mxLbl = 20
character(len=8) :: gLblN(mxLbl)
integer(kind=iwp) :: gLblC(mxLbl), mLbl
real(kind=wp) :: gLblW(mxLbl)

!----------------------------------------------------------------------*
!                                                                      *
!     Define the length of the recognition area for:                   *
!     Commands, 1st level - and 2nd level subcommands                  *
!                                                                      *
!----------------------------------------------------------------------*

integer(kind=iwp), parameter :: lCom = 4, lSub = 4, lParm = 3

!----------------------------------------------------------------------*
!                                                                      *
!     Allocate space to store the title                                *
!                                                                      *
!----------------------------------------------------------------------*

integer(kind=iwp), parameter :: MxTitL = 10
character(len=72) :: Title(MxTitL)
integer(kind=iwp) :: mTit

!----------------------------------------------------------------------*
!                                                                      *
!     Allocate space to store the header of the one-electron           *
!     integral file.                                                   *
!                                                                      *
!----------------------------------------------------------------------*

integer(kind=iwp) :: nSym, nBas(MxSym), nAtoms
real(kind=wp), allocatable :: Coor(:,:)
character :: Header(144)

!----------------------------------------------------------------------*
!     An input vector for the SELEctive keyword.                       *
!----------------------------------------------------------------------*

real(kind=wp) :: TranCoo(3)
integer(kind=iwp) :: nSets
logical(kind=iwp) :: LCumulate
integer(kind=iwp), allocatable :: iSelection(:,:)
logical(kind=iwp), allocatable :: Atoms(:), Bonds(:,:)

public :: Atoms, Bonds, ComCtl, ComStk, ComTab, ComVal, Coor, Header, LCumulate, MxLbl, MxTitL, Title, TranCoo, gLblC, gLblN, &
          gLblW, iSelection, mLbl, mTit, nAtoms, nBas, nCom, nSets, nSym, Cleanup

contains

subroutine Cleanup()
  use stdalloc, only: mma_deallocate
  if (allocated(Coor)) call mma_deallocate(Coor)
  if (allocated(iSelection)) call mma_deallocate(iSelection)
  if (allocated(Atoms)) call mma_deallocate(Atoms)
  if (allocated(Bonds)) call mma_deallocate(Bonds)
end subroutine Cleanup

end module FFPT_global
