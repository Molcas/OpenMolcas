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

subroutine Xprop(short,ifallorb,nirrep,nbas,ntotv,vec,ntoto,occ,ntotd,opel,dout)
!***********************************************************************
!                                                                      *
!     Purpose: the calculation of the average value of an operator,    *
!              whose matrix elements are supplied in the array opel,   *
!              for a set of eigenvectors supplied in vec with occupa-  *
!              tion numbers supplied in occ                            *
!                                                                      *
!     Parameters                                                       *
!                                                                      *
!       short          logical, if (short) then only the total average *
!                      value is calculated and transferred in dout(1)  *
!       ifallorb       logical, if (ifallorb) then print the property  *
!                      of all orbitals (including virtuals) and none is*
!                      weighted by occupation number (S.S.Dong, 2018)  *
!       nirrep         number of irreps                                *
!       nbas(0:nirrep) dimension for each irrep                        *
!       ntotv          the total number of elemnts for all eigenvectors*
!       vec(ntotv)     if (short) then vec stores all lower triangles  *
!                      for all diagonal blocks of the density matrix   *
!                      size: sum(i,i=0,nirrep-1)(nbas(i)*(nbas(i)+1)/2)*
!                      else                                            *
!                      vec stores the eigenvectors                     *
!                      size: sum(i,i=0,nirrep-1)(nbas(i)*nbas(i))      *
!       ntoto          the total number of vectors for all represen-   *
!                      tations=the number of basis functions           *
!       occ(ntoto)     if (short) the occ array is a dummy             *
!                      else                                            *
!                      occ stores the occupation numbers               *
!                      size: sum(i,i=0,nirrep-1)(nbas(i))              *
!       ntotd          the total number of elements in lower triangles *
!                      of all diagonal blocks                          *
!       opel(ntotd)    a storage area for transferring all lower       *
!                      triangles of diagonal blocks of the operator    *
!                      matrix                                          *
!                      size: sum(i,i=0,nirrep-1)(nbas(i)*(nbas(i)+1)/2)*
!       dout           on return if (short) dout(1) contains the total *
!                      average value                                   *
!                      else                                            *
!                      dout(i), i=1,sum(k,k=0,nirrep-1)(nbas(i))       *
!                      contains the orbital contributions (multiplied  *
!                      by the corresponding occupation numbers)        *
!                                                                      *
! Modified by S.S.Dong, 2018, Univ. of Minnesota                       *
! - Enable properties to be printed for all orbitals                   *
! (including virtuals) and not weighted by occupation numbers          *
!***********************************************************************

use Constants, only: Zero, Two
use Definitions, only: wp, iwp

implicit none
logical(kind=iwp), intent(in) :: short, ifallorb
integer(kind=iwp), intent(in) :: nirrep, nbas(0:nirrep-1), ntotv, ntoto, ntotd
real(kind=wp), intent(in) :: vec(ntotv), occ(ntoto), opel(ntotd)
real(kind=wp), intent(inout) :: dout(ntoto)
integer(kind=iwp) :: i, iado, iadout, iadv, icount, iv, iv1, iv2, jCount, nDim2
real(kind=wp) :: dsum
real(kind=wp), external :: DDOT_

if (short) then
  icount = 1
  dsum = DDOT_(ntotd,vec,1,opel,1)
  dout(1) = dsum
else if (.not. ifallorb) then
  ndim2 = 0
  do i=0,nirrep-1
    ndim2 = ndim2+nbas(i)**2
  end do

  iadv = 0
  iado = 0
  iadout = 0
  jCount = 1
  do i=0,nirrep-1
    do iv=1,nbas(i)
      iado = iado+1
      iadout = iadout+1
      dsum = Zero
      icount = jCount
      do iv1=1,nbas(i)
        do iv2=1,iv1-1
          dsum = dsum+Two*vec(iadv+iv1)*vec(iadv+iv2)*opel(icount)
          icount = icount+1
        end do
        dsum = dsum+vec(iadv+iv1)*vec(iadv+iv1)*opel(icount)
        icount = icount+1
      end do
      dout(iadout) = occ(iado)*dsum
      iadv = iadv+nbas(i)
    end do
    jCount = jCount+nBas(i)*(nBas(i)+1)/2
  end do
else if (ifallorb) then
  ndim2 = 0
  do i=0,nirrep-1
    ndim2 = ndim2+nbas(i)**2
  end do

  iadv = 0
  iado = 0
  iadout = 0
  jCount = 1
  do i=0,nirrep-1
    do iv=1,nbas(i)
      iado = iado+1
      iadout = iadout+1
      dsum = Zero
      icount = jCount
      do iv1=1,nbas(i)
        do iv2=1,iv1-1
          dsum = dsum+Two*vec(iadv+iv1)*vec(iadv+iv2)*opel(icount)
          icount = icount+1
        end do
        dsum = dsum+vec(iadv+iv1)*vec(iadv+iv1)*opel(icount)
        icount = icount+1
      end do
      dout(iadout) = dsum
      iadv = iadv+nbas(i)
    end do
    jCount = jCount+nBas(i)*(nBas(i)+1)/2
  end do
end if

!if (short) then
!  write(u6,'(3x,a,f18.10)') 'Total = ',dout(1)
!else
!  ii = 0
!  do i=0,nirrep-1
!    write(u6,'(1x,a,i2)') 'Irrep No.',i
!    write(u6,'(5(3x,i3,i3,f18.10))') (j,dout(ii+j),j=1,nbas(i))
!    ii = ii+nbas(i)
!  end do
!end if

return

end subroutine Xprop
