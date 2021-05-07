************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
************************************************************************
*                                                                      *
* This routine dots the UHF spin density with the property integrals   *
* and trace the resulting matrix for each of the 9 components of the   *
* hyperfine magnetic integrals to obtain the 3 by 3 HFC tensor matrix  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Parameters:                                                          *
* nb     -  Number of total basis functions, input.                    *
* nat    -  Number of atoms, input.                                    *
*                                                                      *
************************************************************************
      Subroutine cmp_hfc(nb, nat)
      Implicit none
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
      Integer nb, nat, nb2, irc, iopt, icomp, toper, lu_one, kdir
      Integer nbtri, ip, jp, isd, itd, stri
      Integer ita, isa, isr, idir, jdir, iat
      Real*8 hfc(3,3), trace, amat(3,3)
      character*8 label
      character*16 nome

* Sizes of the matrices
      nb2 = nb * nb
      nbtri = nb * (nb + 1) / 2
      irc = -1
      iopt = 0
      lu_one = 2
      toper = 255

      call allocate_work(itd, nbtri+4)
      call allocate_work(isd, nb2+4)
      call get_d1sao(Work(itd),nbtri)
      call square(Work(itd), Work(isd), nb, 1, nb)

      stri = 0
      do ip = 0, nb - 1
        stri = isd + ip * nb
        do jp = 0, nb - 1
          if (.not. ip.eq.jp) then
            Work(stri + jp) = 0.5d0 * Work(stri + jp)
          endif
        enddo
      enddo

      call allocate_work(ita, nbtri+4)
      call allocate_work(isa, nb2+4)
      call allocate_work(isr, nb2+4)

      irc = -1
      write(label,'(A,I3)') 'DEBUG', 0
      call opnone(irc,iopt,'ONEINT',lu_one)
      if (irc.ne.0) goto 999

      do iat = 1, nat
        icomp = 0
        do idir = 1, 3
          do jdir = 1, 3
            icomp = icomp + 1
            trace = 0
            write(label,'(A,I3)') 'MAGXP', iat
            irc = -1
            call rdone(irc,iopt,label,icomp,work(ita),toper)
            if (irc.ne.0) goto 999
            call square(Work(ita), Work(isa), nb, 1, nb)
            call dgemm_('N', 'N', nb, nb, nb, 1.0d0, Work(isd),
     &                  nb, Work(isa), nb, 0.0d0, Work(isr), nb)
            do kdir = 1, nb
              trace = trace + Work(isr + nb * (kdir - 1) + kdir - 1)
            enddo
            hfc(idir,jdir) = trace
          enddo
        enddo

!  spind field reduction to get the amat elements match
        amat(1,1) = hfc(2,2) + hfc(3,3)
        amat(2,2) = hfc(1,1) + hfc(3,3)
        amat(3,3) = hfc(1,1) + hfc(2,2)
        do idir = 1, 3
            do jdir = 1, 3
                if(idir.ne.jdir) amat(idir,jdir) = -hfc(jdir,idir)
            enddo
        enddo

        write(6,*) ''
        write(6,'(A,I3)') 'Hyperfine coupling tensor matrix for atom
     &                    :', iat
        write(6,*) ''
        write(6,'(A,A)')    '   --------------------------------------',
     &                      '-------------------'
        do idir = 1, 3
          write(6,'(3E20.10)') (-amat(idir,jdir),jdir=1,3)
        enddo
        write(6,'(A,A)')    '   --------------------------------------',
     &                      '-------------------'
        write(6,'(A)') ''
      enddo
      Call Add_Info('AMAT',AMAT,9,5)

      call free_work(itd)
      call free_work(isd)
      call free_work(ita)
      call free_work(isa)
      call free_work(isr)
      call clsone(irc,iopt)

      return

 999  continue
      write(6,*) ' *** Error in subroutine cmp_hfc ***'
      write (6,'(A,A)') '     Label = ', Label
      Call Abend
      call free_work(itd)
      call free_work(isd)
      call free_work(ita)
      call free_work(isa)
      call free_work(isr)
      call clsone(irc,iopt)

      End subroutine cmp_hfc
