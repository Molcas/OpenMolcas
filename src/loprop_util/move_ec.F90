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

subroutine Move_EC(rMP,EC,Scratch_New,Scratch_Org,xrMP,xxrMP,xnrMP,lMax,A,B,nij,nElem,Coor,nAtoms,Q_Nuc,C_o_C,nPert, &
                   Bond_Threshold,iANr,T_Values,iT_Sets,iWarnings,Num_Warnings,Opt_Method,iPlot,iPrint)

use Real_Spherical

implicit real*8(A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
#include "status.fh"
real*8 rMP(nij,0:nElem-1,0:nPert-1), EC(3,nij)
real*8 Scratch_Org(nij*(2+lMax+1))
real*8 Scratch_New(nij*(2+lMax+1)), xrMP(nij,nElem), C_o_C(3)
real*8 xxrMP(nij,nElem), T_Values(nij)
real*8 xnrMP(nij,nElem), Coor(3,nAtoms), CX(3), Q_Nuc(nAtoms)
real*8 A(3,nij), B(3,nij)
logical Bond_OK, Check_Bond
character*12 Opt_Method
integer iAnr(nAtoms), iT_Sets(nij), iWarnings(nij)
!parameter (mxAtoms=500)

call Sphere(lMax)

call dcopy_(3*nij,[Zero],0,B,1)

call dCopy_(nij*nElem,[Zero],0,xnrMP,1)
ij = 0
do iAtom=1,nAtoms
  do jAtom=1,iAtom
    ij = ij+1
    if (iAtom == jAtom) then
      xnrMP(ij,1) = Q_nuc(iAtom)
      CX(1) = Coor(1,iAtom)
      CX(2) = Coor(2,iAtom)
      CX(3) = Coor(3,iAtom)
      call ReExpand(xnrMP,nij,nElem,CX,C_o_C,ij,lMax)
    end if
  end do
end do

do iAtom=1,nAtoms
  ii = iAtom*(iAtom+1)/2
  do jAtom=iAtom,1,-1
    jj = jAtom*(jAtom+1)/2
    ij = iAtom*(iAtom-1)/2+jAtom
    call dcopy_(3*nij,EC,1,A,1)

    if (iAtom == jAtom) then
      ! Atomic domain
      ! Do nothing for now - put appropriate expresion in later.
      iT_sets(ij) = 1
      T_values(ij) = 0.0d0
      B(1,ij) = EC(1,ij)
      B(2,ij) = EC(2,ij)
      B(3,ij) = EC(3,ij)
      !?? end if
    else
      ! Bond domain
      !
      ! Check whether bond is too long - if so skip it.

      Bond_OK = Check_Bond(EC(1,ii),EC(1,jj),iANr(iAtom),iANr(jAtom),Bond_Threshold)
      if (.not. Bond_OK) then
        T_values(ij) = 0.0d0
        B(1,ij) = EC(1,ij)
        B(2,ij) = EC(2,ij)
        B(3,ij) = EC(3,ij)
      else

        ! Bond is real bond - proceed

        iT_sets(ij) = 1
        if (Opt_Method(1:9) == 'Multipole') then
          call Rotate_Dipole(rMP,EC,nij,nElem,nPert,ij,ii,jj,Dipole_Rot_A,Dipole_Rot_B,Dipole_Rot_AB,R_A,R_B)

          call Find_Dipole_Center(rMP(ii,0,0),rMP(jj,0,0),Dipole_Rot_A,Dipole_Rot_B,xnrMP(ii,1),xnrMP(jj,1),R_A,R_B,EC(1,ii), &
                                  EC(1,jj),EC(1,ij),T_Values(ij),iPlot)

          B(1,ij) = EC(1,ij)+T_Values(ij)*(EC(1,jj)-EC(1,ii))
          B(2,ij) = EC(2,ij)+T_Values(ij)*(EC(2,jj)-EC(2,ii))
          B(3,ij) = EC(3,ij)+T_Values(ij)*(EC(3,jj)-EC(3,ii))
        else
          call Min_Mult_Error(EC,A,B,EC(1,ii),EC(1,jj),rMP,xrMP,xxrMP,xnrMP,lMax,nij,nElem,iAtom,jAtom,nAtoms,nPert,C_o_C, &
                              Scratch_New,Scratch_Org,iPlot,T_Values,iWarnings,Num_Warnings)
        end if
      end if

    end if
  end do
end do

do ij=1,nij
  do iPert=0,nPert-1
    call ReExpand(rMP(1,0,iPert),nij,nElem,EC(1,ij),B(1,ij),ij,lMax)
  end do
end do
call dcopy_(3*nij,B,1,EC,1)

call Sphere_Free()

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(iPrint)

end subroutine Move_EC
