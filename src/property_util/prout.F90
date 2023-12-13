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

subroutine PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,PrTot,iPL,iPoint,ifallorb)
!***********************************************************************
!                                                                      *
!     purpose: printing of tables for different tensor properties      *
!                                                                      *
!     Short           logical option for either Short (Short           *
!                     =.true., the total electronic contribution)      *
!                     long (Short=.false., orbital contributions)      *
!                     output                                           *
!     sig             the sign factor for the electronic contri-       *
!                     bution, positive for ndiamagnetic shielding      *
!     nIrrep          the number of irreducible representations        *
!     nBas            the number of functions in each representa-      *
!     (0:nIrrep-1)    tion                                             *
!     nTot            the total number of elements supplied for        *
!                     each component of the property tensor; equal     *
!                     either to 1 (total electronic and total nuc-     *
!                     lear contributions) or to the dimension of       *
!                     the basis set.                                   *
!     Occ(1:nTot)     Occupation numbers for all eigenvectors,         *
!                     a dummy for Short outputs                        *
!     ThrSV           threshold for Occupation numbers; if             *
!                     Occ(i).le.ThrSV the contribution will not        *
!                     be printed                                       *
!     PrEl(1:nTot,    matrix elements for all components 1,2,...,      *
!          1:maxlab)  maxlab, nTot entries for each component          *
!                                                                      *
!     PrNu(1:maxlab)  nuclear contributions for each component         *
!     maxlab          total number of cartesian components             *
!     labs(1:maxlab)  labels for each component                        *
!     TotEl(1:6)      auxiliary storage area                           *
!     PrTot(1:maxlab) Total value for each component                   *
!     iPL             Print level                                      *
!     iPoint          The number of the center                         *
!     ifallorb        logical option for whether the property of       *
!                     all orbitals are printed (and not weighted by    *
!                     occupation number)in property calculation when   *
!                     short=.false. (S.S.Dong, 2018)                   *
!                                                                      *
! 2000 Dept. of Chem. Phys., Univ. of Lund, Sweden                     *
! Modified by S.S.Dong, 2018, Univ. of Minnesota                       *
! - Enable properties to be printed for all orbitals                   *
! (including virtuals) and not weighted by occupation numbers          *
!***********************************************************************

use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: Short, ifallorb
integer(kind=iwp), intent(in) :: nIrrep, nBas(0:nIrrep-1), nTot, maxlab, iPL, iPoint
real(kind=wp), intent(in) :: sig, Occ(nTot), ThrSV, PrEl(nTot,maxlab), PrNu(maxlab)
character(len=16), intent(in) :: labs(maxlab)
real(kind=wp), intent(out) :: PrTot(maxlab)
integer(kind=iwp) :: i, icount, ii, j, jcount, jj, nDec, nLog10
real(kind=wp) :: TotEl(6), val
character(len=132) :: outlab
character(len=18) :: Format1
character(len=16) :: lab
character(len=14) :: Format3
character(len=13) :: Format2
character(len=11) :: Format4
character(len=4) :: lbb4
character(len=5) :: lab5

Format1 = '(2i5,f14.8,6f16.8)'
Format2 = '(1x,a,6f16.8)'
Format3 = '(1x,a,6f16.8/)'
Format4 = '(a5,6f16.8)'
val = Zero
do i=1,maxlab
  val = max(val,abs(PrNu(i)))
  do j=1,nTot
    val = max(val,abs(PrEl(j,i)))
  end do
end do
val = max(val,One)
nLog10 = max(1,int(log10(val)+One)+1)
nDec = min(14-nLog10,8)
write(Format1(17:17),'(I1)') nDec
write(Format2(12:12),'(I1)') nDec
write(Format3(12:12),'(I1)') nDec
write(Format4(10:10),'(I1)') nDec

if ((.not. Short) .and. (.not. ifallorb)) then
  write(u6,'(A,ES9.2/)') ' orbital contributions printed for occupation numbers >',ThrSV
else if ((.not. Short) .and. ifallorb) then
  write(u6,'(A)') ' orbital properties printed for all occupation numbers'
end if

do i=1,maxlab,6

  do j=1,132
    outlab(j:j) = ' '
  end do

  if (Short) then
    if (iPL >= 3) then
      if ((maxlab > 1) .or. (labs(1) /= ' ')) then
        write(outlab,'(1x,a,6a16)') 'Component              ',(labs(j),j=i,min(i+5,maxlab))
        write(u6,'(a)') trim(outlab)
      end if
    end if
    jcount = 0
    do j=i,min(i+5,maxlab)
      jcount = jcount+1
      TotEl(jcount) = PrEl(1,j)
    end do
  else
    do j=1,6
      TotEl(j) = Zero
    end do
    write(outlab,'(1x,a,6a16)') 'Irrep  Orb   Occupation',(labs(j),j=i,min(i+5,maxlab))
    write(u6,'(a)') trim(outlab)
    lbb4 = '    '
    write(outlab,'(33a4)') (lbb4,j=1,33)
    lbb4 = '----'
    lab = '----------------'
    write(outlab,'(6a4,6a16)') (lbb4,j=1,6),(lab,j=i,min(i+5,maxlab))
    write(u6,'(1x,a)') trim(outlab)
    icount = 0
    do ii=0,nIrrep-1
      do jj=1,nBas(ii)
        icount = icount+1
        jcount = 0
        do j=i,min(i+5,maxlab)
          jcount = jcount+1
          if (ifallorb) then
            TotEl(jcount) = TotEl(jcount)+PrEl(icount,j)*Occ(icount)
          else
            TotEl(jcount) = TotEl(jcount)+PrEl(icount,j)
          end if
        end do
        if (ifallorb .or. (Occ(icount) > ThrSV)) then
          write(u6,Format1) ii+1,jj,Occ(icount),(sig*PrEl(icount,j),j=i,min(i+5,maxlab))
        end if
      end do
    end do
    write(u6,'(1x,a)') trim(outlab)
  end if

  do j=i,min(i+5,maxlab)
    PrTot(j) = sig*TotEl(j-i+1)+PrNu(j)
  end do
  if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
    write(u6,Format2) 'Total electronic       ',(sig*TotEl(j-i+1),j=i,min(i+5,maxlab))
    write(u6,Format2) 'Total nuclear          ',(PrNu(j),j=i,min(i+5,maxlab))
    write(u6,Format3) 'Total                  ',(PrTot(j),j=i,min(i+5,maxlab))
  else
    lab5 = '     '
    if ((iPoint > 0) .and. (i == 1)) write(lab5,'(I5)') iPoint
    write(u6,Format4) lab5,(PrTot(j),j=i,min(i+5,maxlab))
  end if
end do

return

end subroutine PrOut
