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

subroutine InsReqo2v4(NaGrp,NbeGrp,NaSGrp)
! this routine does:
! Inspect W3 and W4 files requirements of o2v4 procedure
! on this node. It checks which of the W3 and W4 files
! are used on this node

use Index_Functions, only: nTri_Elem
use chcc_global, only: ABID, DimGrpa, DimGrpbe, DimSGrpbe, GrpaLow, GrpaUp, GrpbeLow, GrpbeUp, InqW3, InqW4
use Para_Info, only: MyRank
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: NaGrp, NbeGrp, NaSGrp
integer(kind=iwp) :: adda, addb, addbe, addbepp, addga, addgapp, aGrp, aSGrp, beGrp, beSGrp, bGrp, bSGrp, bSGrpUp, dima, dimb, &
                     dimbe, dimga, gaGrp, gaSGrp, gaSGrpUp, i, LenW3, LenW4, NSGrp

!* Initial Set InqW3,InqW4 = F
NSGrp = NaGrp*NaSGrp
InqW3(1:nTri_Elem(NSGrp),1:NSGrp) = .false.
InqW4(1:nTri_Elem(NSGrp),1:nTri_Elem(NSGrp)) = .false.
LenW3 = 0
LenW4 = 0

!* cycle over all groups of (a>=b)
adda = 0
do aGrp=1,NaGrp
  dima = DimGrpa(aGrp)

  !## test, if on this node at least one combination with this
  !   aGrp is scheduled. Skip if no
  i = sum(ABID(myRank,aGrp,1:NaGrp))

  if (i /= 0) then
    addb = 0
    do bGrp=1,aGrp
      dimb = DimGrpa(bGrp)

      !## test, if this a'b' combination is planed to be run on this node
      if (ABID(myRank,aGrp,bGrp) /= 0) then

        ! cycle over all groups of (be>=ga)
        addbe = 0
        do beGrp=1,NbeGrp
          dimbe = DimGrpbe(beGrp)
          addga = 0
          do gaGrp=1,beGrp
            dimga = DimGrpbe(gaGrp)

            ! cycle over all subgroups of (be>=ga)'
            addbepp = addbe
            do beSGrp=GrpbeLow(beGrp),GrpbeUp(beGrp)

              if (beGrp == gaGrp) then
                gaSGrpUp = beSGrp
              else
                gaSGrpUp = GrpbeUp(gaGrp)
              end if

              addgapp = addga
              do gaSGrp=GrpbeLow(gaGrp),gaSGrpUp

                ! cycle over all subgroups of (a>=b)'
                do aSGrp=GrpaLow(aGrp),GrpaUp(aGrp)
                  if (aGrp == bGrp) then
                    bSGrpUp = aSGrp
                  else
                    bSGrpUp = GrpaUp(bGrp)
                  end if
                  do bSGrp=GrpaLow(bGrp),bSGrpUp

                    ! clasical reading or generating of (VV|VV) integrals

                    ! term W1(a",be",b",ga") <- -(b"ga"|a"i) . T1(i,be")
                    ! Get Ww(b",ga",a",i) (Wx used as Aux)
                    call InsReaW3(bSGrp,gaSGrp,aSGrp,LenW3)

                    ! Upgrade W1(a",be",b",ga") <<- (a",be"|b",ga")
                    call InsReaW4(aSGrp,beSGrp,bSGrp,gaSGrp,LenW4)

                    if (aSGrp /= bSGrp) then
                      ! term W2(b",be",a",ga") <- -(a"ga"|b"i) . T1(i,ga")
                      ! Get Ww(a",ga",b",i) (Wx used as Aux)
                      call InsReaW3(aSGrp,gaSGrp,bSGrp,LenW3)

                      ! term W2(b",be",a",ga") <<- -(b"be"|a"i) . T1(i,ga")
                      ! Get Ww(b",be",a",i) (Wx used as Aux)
                      call InsReaW3(bSGrp,beSGrp,aSGrp,LenW3)

                      ! Add W2(b",be",a",ga") <<- (b",be"Ia",ga")
                      ! Upgrade:W2, destroy:Ww
                      call InsReaW4(bSGrp,beSGrp,aSGrp,gaSGrp,LenW4)
                    end if

                  ! end cycle over all subgroups of (a>=b)'
                  end do
                end do

                ! end cycle over all subgroups of (be>=ga)'
                addgapp = addgapp+DimSGrpbe(gaSGrp)
              end do
              addbepp = addbepp+DimSGrpbe(beSGrp)
            end do

            ! end cycle over all groups of (be>=ga)
            addga = addga+dimga
          end do
          addbe = addbe+dimbe
        end do

      end if
      ! end cycle over all groups of (a>=b)
      addb = addb+dimb
    end do
  end if
  adda = adda+dima
end do

return

end subroutine InsReqo2v4
