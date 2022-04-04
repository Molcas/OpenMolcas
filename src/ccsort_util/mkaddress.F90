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

subroutine mkaddress(NOIPSB)
!FUE It is strictly forbidden to use any construct to fetch disk
!FUE addresses which do not make use of subroutines provided in
!FUE the MOLCAS utilities as otherwise transparency and portability
!FUE is lost.
!
! this routine prepares information arrays
! for integral blocks per symmetry in primitive form
! it defines:
! noipsb(nijkl) - number of integrals per symmetry block
! idispsb(nijkl)- initial addresses for each symmetry block
! typ(pa,qa,ra) - type of symmetry block
! types of (ij|kl):
! 1 - si=sk, si=sj, sk=sl
! 2 - si=sk, si=sj, sk>sl
! 3 - si=sk, si>sj, sk=sl
! 4 - si=sk, si>sj, sk>sl
! 5 - si>sk, si=sj, sk=sl
! 6 - si>sk, si=sj, sk>sl
! 7 - si>sk, si>sj, sk=sl
! 8 - si>sk, si>sj, sk>sl
!
! idis(pa,qa,ra)- initial addresses for given symmetry block
! np(pa,qa,ra)  - position of p index in original block (ij|kl) (on tape)
! nq(pa,qa,ra)  - position of q index in original block (ij|kl) (on tape)
! nr(pa,qa,ra)  - position of r index in original block (ij|kl) (on tape)
! ns(pa,qa,ra)  - position of s index in original block (ij|kl) (on tape)
!
! N.B. typ,idis,np,nq,nr,ns are imported from ccsort_global

use ccsort_global, only: fullprint, idis, LUINTM, NORB, np, nq, nr, ns, NSYM, typ
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: NOIPSB(106)
#include "tratoc.fh"
integer(kind=iwp) :: idishelp, idistemp, ilow, IND, INDT, iold, ISPQRS, iup, jlow, jold, jup, kold, kup, lold, lup, N_INT, norbp, &
                     nsi, nsij, nsijk, nsj, nsk, nsl, NSLM, p, pa, q, qa, r, ra, s, sense !, idispsb(106)
real(kind=wp) :: dum(1)

!FUE - pick the start addresses of each symmetry allowed integral
!FUE   block from the tranformed two electron integral file
idistemp = 0
call iDaFile(LUINTM,0,iTraToc,nTraToc,idistemp)
!FUE - the option 0 in the call to dafile does not make any I/O
!FUE   but just returns the updated disk address
!FUE - at this point idistemp points to the first record

IND = 0
INDT = 0
!FUE idistemp = 1
!FUE idisadd = 150

typ(1:NSYM,1:NSYM,1:NSYM) = 0

if (fullprint > 0) then
  write(u6,'(6X,A)') 'Transformed integral blocks:'
  write(u6,'(6X,A)') '----------------------------'
  write(u6,*)
  write(u6,'(6X,A)') 'block  symmetry      no. of        no. of '
  write(u6,'(6X,A)') ' no.    spec.        orbitals     integrals'
  write(u6,'(6X,A)') '-------------------------------------------'
end if

ISPQRS = 0
do NSI=1,NSYM
  do NSJ=1,NSI
    NSIJ = MUL(NSI,NSJ)
    do NSK=1,NSI
      NSIJK = MUL(NSK,NSIJ)
      NSLM = NSK
      if (NSK == NSI) NSLM = NSJ
      do NSL=1,NSLM
        if (NSIJK /= NSL) cycle
        NORBP = NORB(NSI)*NORB(NSJ)*NORB(NSK)*NORB(NSL)
        if (NORBP == 0) cycle
        ISPQRS = ISPQRS+1

        ! def

        ! redefine indices from Parr to Dirac notation <p,q|r,s> <-> (i,j|k,l)

        p = nsi
        q = nsk
        r = nsj
        s = nsl

        ! def. sense

        ! type (ij|kl)
        ! 1 - si=sk, si=sj, sk=sl
        ! 2 - si=sk, si=sj, sk>sl
        ! 3 - si=sk, si>sj, sk=sl
        ! 4 - si=sk, si>sj, sk>sl
        ! 5 - si>sk, si=sj, sk=sl
        ! 6 - si>sk, si=sj, sk>sl
        ! 7 - si>sk, si>sj, sk=sl
        ! 8 - si>sk, si>sj, sk>sl

        if (nsi == nsk) then
          if (nsi == nsj) then
            if (nsk == nsl) then
              sense = 1
            else
              sense = 2
            end if
          else
            if (nsk == nsl) then
              sense = 3
            else
              sense = 4
            end if
          end if
        else
          if (nsi == nsj) then
            if (nsk == nsl) then
              sense = 5
            else
              sense = 6
            end if
          else
            if (nsk == nsl) then
              sense = 7
            else
              sense = 8
            end if
          end if
        end if

        if (nsijk /= nsl) then
          sense = 0
        else if (NORB(NSI)*NORB(NSJ)*NORB(NSK)*NORB(NSL) == 0) then
          sense = 0
        end if

        !1: perm <pq|rs> -> <pq|rs>

        pa = p
        qa = q
        ra = r
        typ(pa,qa,ra) = sense
        idis(pa,qa,ra) = idistemp
        np(pa,qa,ra) = 1
        nq(pa,qa,ra) = 3
        nr(pa,qa,ra) = 2
        ns(pa,qa,ra) = 4

        !2: perm <pq|rs> -> <rq|ps> 1-3

        pa = r
        qa = q
        ra = p
        typ(pa,qa,ra) = sense
        idis(pa,qa,ra) = idistemp
        np(pa,qa,ra) = 2
        nq(pa,qa,ra) = 3
        nr(pa,qa,ra) = 1
        ns(pa,qa,ra) = 4

        !3: perm <pq|rs> -> <ps|rq> 2-4

        pa = p
        qa = s
        ra = r
        typ(pa,qa,ra) = sense
        idis(pa,qa,ra) = idistemp
        np(pa,qa,ra) = 1
        nq(pa,qa,ra) = 4
        nr(pa,qa,ra) = 2
        ns(pa,qa,ra) = 3

        !4: perm <pq|rs> -> <pq|rs> 1-3,2-4

        pa = r
        qa = s
        ra = p
        typ(pa,qa,ra) = sense
        idis(pa,qa,ra) = idistemp
        np(pa,qa,ra) = 2
        nq(pa,qa,ra) = 4
        nr(pa,qa,ra) = 1
        ns(pa,qa,ra) = 3

        !5: perm <pq|rs> -> <qp|sr> 1-2,3-4

        pa = q
        qa = p
        ra = s
        typ(pa,qa,ra) = sense
        idis(pa,qa,ra) = idistemp
        np(pa,qa,ra) = 3
        nq(pa,qa,ra) = 1
        nr(pa,qa,ra) = 4
        ns(pa,qa,ra) = 2

        !6: perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <sp|qr> 1-3

        pa = s
        qa = p
        ra = q
        typ(pa,qa,ra) = sense
        idis(pa,qa,ra) = idistemp
        np(pa,qa,ra) = 4
        nq(pa,qa,ra) = 1
        nr(pa,qa,ra) = 3
        ns(pa,qa,ra) = 2

        !7: perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <qr|sp> 2-4

        pa = q
        qa = r
        ra = s
        typ(pa,qa,ra) = sense
        idis(pa,qa,ra) = idistemp
        np(pa,qa,ra) = 3
        nq(pa,qa,ra) = 2
        nr(pa,qa,ra) = 4
        ns(pa,qa,ra) = 1

        !8: perm <pq|rs> -> <qp|sr> 1-2,3-4 -> <sr|qp> 1-3,2-4

        pa = s
        qa = r
        ra = q
        typ(pa,qa,ra) = sense
        idis(pa,qa,ra) = idistemp
        np(pa,qa,ra) = 4
        nq(pa,qa,ra) = 2
        nr(pa,qa,ra) = 3
        ns(pa,qa,ra) = 1

        !idispsb(ispqrs) = idistemp
        idishelp = 0

        !***************************************************************

        ! LOOP OVER THE USED ORBITALS OF EACH SYMMETRY BLOCK
        ! THIS LOOPING IS COPIED FROM THE MANUAL OF THE 4-INDEX
        ! TRANSFORMATION PROGRAM
        !
        ! N_INT COUNTS INTEGRAL LABELS IN THE GIVEN SYMMETRY BLOCK

        N_INT = 0
        KUP = NORB(NSK)
        do KOLD=1,KUP

          LUP = NORB(NSL)
          if (NSK == NSL) LUP = KOLD
          do LOLD=1,LUP

            ILOW = 1
            if (NSI == NSK) ILOW = KOLD
            IUP = NORB(NSI)
            do IOLD=ILOW,IUP

              JLOW = 1
              if ((NSI == NSK) .and. (IOLD == KOLD)) JLOW = LOLD
              JUP = NORB(NSJ)
              if (NSI == NSJ) JUP = IOLD
              do JOLD=JLOW,JUP

                IND = IND+1
                INDT = INDT+1
                N_INT = N_INT+1
                idishelp = idishelp+1
                if (idishelp > nTraBuf) then
                  !FUE - all integrals in the reord were processed, hence
                  !FUE   update the disk address by a dummy I/O
                  dum(1) = Zero
                  call dDaFile(LUINTM,0,dum,nTraBuf,idistemp)
                  !FUE idistemp = idistemp+idisadd
                  idishelp = 1
                end if

                if (IND >= nTraBuf) IND = 0

              end do
            end do
          end do
        end do
        if (idishelp > 0) then
          !FUE - all integrals in the reord were processed, hence
          !FUE   update the disk address by a dummy I/O
          dum(1) = Zero
          call dDaFile(LUINTM,0,dum,nTraBuf,idistemp)
          !FUE idistemp = idistemp+idisadd
        end if

        ! WRITING THE LAST RECORD OF LABELS IN THE GIVEN SYMMETRY BLOCK
        ! RECORDS ON LUPACK ARE FORMATTED TO 28KB LENGTH

        if (IND /= 0) IND = 0

        NOIPSB(ISPQRS) = N_INT

        !***************************************************************

        if (fullprint > 0) write(u6,'(6X,I5,2X,4I2,2X,4I4,2X,I8)') ISPQRS,NSI,NSJ,NSK,NSL,IUP,JUP,KUP,LUP,N_INT

        !***************************************************************

      end do
    end do
  end do
end do
if (fullprint > 0) write(u6,'(6X,A)') '-------------------------------------------'

return

end subroutine mkaddress
