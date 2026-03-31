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

subroutine AppDerCpl(natom,nST,ChgNuc,Prop,DerCpl,HAM)
! Approximate derivative couplings:         <\Psi_I|\nabla H|\Psi_J>
!                                    f_IJ =  ----------------------
!                                                 E_J - E_I
!
! If the wfn are real-valued: f_II = 0 ; f_JI = - f_IJ -> lower triangular storage

use rassi_aux, only: ipglob
use Cntrl, only: ICOMP, NPROP, NSTATE, PNAME
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: natom, nST
real(kind=wp) :: ChgNuc(natom), Prop(nState,nState,NProp), DerCpl(nST,3,natom), Ham(Nstate,Nstate)
integer(kind=iwp) :: IST, ISTA, JSTA, KPROP, KXYZ, LAT
real(kind=wp) :: EI, EJ, SumX, SumY, SumZ

nST = nState*(nState+1)/2
call FZero(DerCpl,3*natom*nST)
do iSta=1,nState-1
  Ei = Ham(iSta,iSta)
  do jSta=iSta+1,nState
    Ej = Ham(jSta,jSta)
    iST = iSta*(jSta-1)/2+jSta
    write(u6,1000) iSta,jSta,Ej-Ei
    do kProp=1,nProp
      if (PName(kProp)(1:3) == 'EF1') then
        read(PName(kProp)(5:8),'(i4)') lAt
        DerCpl(iST,IComp(kProp),lAt) = Prop(iSta,jSta,kProp)*ChgNuc(lAt)/(Ej-Ei)
      end if
    end do
    SumX = Zero
    SumY = Zero
    SumZ = Zero
    do lAt=1,natom
      write(u6,1100) lAt,(DerCpl(iST,kXYZ,lAt),kXYZ=1,3)
      SumX = SumX+DerCpl(iST,1,lAt)
      SumY = SumY+DerCpl(iST,2,lAt)
      SumZ = SumZ+DerCpl(iST,3,lAt)
    end do
    if (IPGLOB >= 4) write(u6,1200) SumX,SumY,SumZ
  end do
end do

1000 format(/,' Approximate derivative couplings for states ',2i3,/,' Energy difference = ',F15.8,/, &
            '   Atom          X              Y              Z')
1100 format(i7,3f15.8)
1200 format('   Sum:',3f15.8)

end subroutine AppDerCpl
