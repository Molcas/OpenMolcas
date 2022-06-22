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
        subroutine autoSegmentation(Nprocs, maxspace,
     & Jal1, Jal2,
     & NvGrp, NvSGrp, NchBlk, wrksize, maxdim)
c
c           Issues:
c
c           1) o2v4: n'(n'+1)/2 + n' "poltaskov"
c           2) odhadnut overhead? on-the-fly vs precalculate?
c
        implicit none

c MaxGrp, MaxSGrp, ...
#include "chcc1.fh"

        integer Nprocs, NvGrp, NvSGrp, NchBlk
        real*8 eff_thrs, eff
        parameter (eff_thrs = 80.0d0)

        logical requireEfficiency
        parameter (requireEfficiency=.true.)
c
        integer Jal1,Jal2
        integer wrksize, maxspace, maxdim
c
        write (6,*)
        write (6,*) "==============================="
        write (6,*) 'Autogenerating segmentation'
        write (6,'(A,i8)') ' Nprocs: ',Nprocs
        write (6,*)
c
c reset small segmentation
        NvSGrp = 1
c
c get rough estimated of N', e.g. Nprocs == nntasks
        NvGrp = int(sqrt(2.0d0*Nprocs))
        if (printkey.ge.2) then
          write (6,'(A,i4)')
     & " 1st Np estimate: ",
     & NvGrp
        end if

c if no such Np, get the smallest Np for which nntasks > Nprocs
        if ((NvGrp*NvGrp/2.0d0).lt.Nprocs) then
          NvGrp = NvGrp + 1
          if (printkey.ge.2) then
            write (6,'(A,i4)')
     & " Corrected 1st Np estimate: ",
     & NvGrp
          end if
        end if

c make correction for efficiency, if required
        if (requireEfficiency) then
           call findNextEffSeg(NvGrp, eff,
     & Nprocs, eff_thrs, maxGrp, printkey)
           if (printkey.ge.2) then
             write (6,'(A,i4,A,f6.2)')
     & " (maybe) further correction for efficiency: ",
     & NvGrp,", efficiency: ",eff*100
           end if
        endif

        write (6,*)

c does the segmentation fit into memory (start with Npp = 1)?
        call checkMem(NvGrp, NvSGrp, NchBlk,
     & Jal1, Jal2, wrksize, maxdim)
c
12      continue
        if (wrksize.gt.maxspace) then
c
           if (printkey.ge.10) then
             write (6,'(A,i13)') " Not enough memory. Max: ",
     & maxspace,", Current: ",wrksize
           end if
c
c          increase small segmentation, if possible
c
           if ((NvSGrp.lt.8).and.
     & ((NvGrp*(NvSGrp+1)).le.maxSGrp)) then
             NvSGrp = NvSGrp + 1
             if (printkey.ge.10) then
               write (6,'(A,i4)') ' Npp increased: ',NvSGrp
             end if
c
c          if not, increase large segmentation and reset small
c
           else
c
c            reset small segmentation in any case
             NvSGrp = 1
c
c            increment large segmentation
             NvGrp = NvGrp + 1
             if (printkey.ge.10) then
               write (6,'(A,i4)') " Np increased: ",NvGrp
             end if
c
c            make correction for efficiency, if required
             if (requireEfficiency) then
               call findNextEffSeg(NvGrp, eff,
     & Nprocs, eff_thrs, maxGrp, printkey)
                if (printkey.ge.10) then
                  write (6,'(A,i4)')
     & " Np increased (corrected for efficiency): ",NvGrp
                end if
             endif
c
c            check new large segmentation
c
             if (NvGrp.gt.maxGrp) then
               write (6,*)
               write (6,*) " No suitable segmentation found, quiting"
               call abend
             endif
c
           endif ! increase small or large
c
c          print final results and
c          recompute memory requirements
           if (printkey.ge.10) then
             write (6,'(2(A,i4))') " Current Np: ",NvGrp,', Npp: ',
     & NvSgrp
           end if
c
           call checkMem(NvGrp, NvSGrp, NchBlk,
     & Jal1, Jal2, wrksize, maxdim)
           goto 12

        endif ! we fit to memory

        write (6,*)
        write (6,'(A,2(i4),i8)')
     & " Final segmentation (Large/Small/Cholesky): ",
     & NvGrp,NvSGrp,NchBlk
        write (6,*) "==============================="
        write (6,*)

        return
        end
