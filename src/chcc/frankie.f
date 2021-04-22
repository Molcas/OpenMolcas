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
        subroutine frankie(nfro,no,nv,printkey)
c
        use Data_Structures, only: DSBA_Type
        use Data_Structures, only: Allocate_DSBA, Deallocate_DSBA
        implicit none
c
        integer nbas,norb,nocc,nfro,ndel
        integer no,nv
        integer printkey
c
        Type (DSBA_Type) CMO
        integer rc
c
        real*8  FracMem
#include "chotime.fh"
        integer idum(1)
c
c.1 - get the info on  nBas, nOrb, nOcc. Use nFro from input
c
c# nbas = nfro + nocc + nvirt + ndel
c
         Call Get_iArray('nBas',idum,1)
         nBas=idum(1)
         Call Get_iArray('nOrb',idum,1)
         nOrb=idum(1)
         Call Get_iArray('nIsh',idum,1) ! in general > no
         nOcc=idum(1)

         ndel=nbas-no-nv-nfro

        if (printkey.ge.10) then
            write (6,*) 'nbas = ',nbas
            write (6,*) 'norb = ',norb
            write (6,*) 'nocc = ',nocc
            write (6,*) 'nfro = ',nfro
            write (6,*) 'no   = ',no,' (nocc-nfro)'
            write (6,*)
            write (6,*) 'ndel = ',ndel
        end if

        if ( (no+nfro+nv+ndel).ne.nbas ) then
          write (6,*) 'Problem '
          write (6,*) 'nbas from Runfile : ',nbas
          write (6,*) 'nbas control      : ',nfro+no+nv+ndel
          call abend()
        end if
c
        timings=.False.
        if (printkey.gt.1) timings=.True.
c
c.2 - allocate space for CMO with removed SCF deleted and frozen orbitals
c     final ordering of indexes : (o+v,nbas)
c
        Call Allocate_DSBA(CMO,[no+nv],[nbas],1)
        if (printkey.ge.10) then
        write (6,*) 'Dopice 1 - Allo'
        end if
c
c.3 - read CMO
        call read_mo(Cmo,nfro,no,nv,ndel,nbas,nOrb)
c.3 - invert the CMO matrix
        FracMem=0.0d0 ! in a parallel run set it to a sensible value
        rc=0
        Call Cho_X_init(rc,FracMem) ! initialize cholesky info
        if (printkey.ge.10) then
        write (6,*) 'Dopice 2 ',rc
        end if

        call CHO_CC_drv(rc,CMO)
        if (printkey.ge.10) then
        write (6,*) 'Dopice 3 '
        end if

        Call Cho_X_final(rc)
        if (printkey.ge.10) then
        write (6,*) 'Dopice 4 '
        end if
c
        if (rc.ne.0) then
          write (6,*) 'cho_cc_drv failed'
          call abend()
        end if
c
c.  -  deallocate CMO
        Call Deallocate_DSBA(CMO)
c
        return
        end
