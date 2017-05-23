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
#ifdef _HDF5_

      subroutine one2h5_ovlmat(fileid, nsym, nbas)
*     SVC: read basic molecular information from the RunFile
*     and 1-electron integral file and write it to the HDF5
*     file specified with fileid.
*     This routine does nothing if HDF5 is not supported.
*
*     Attributes:
*       NSYM, IRREPS, POTNUC, NBAS
*     Datasets:
*       BASIS_LABELS, OVLMAT

      implicit none
      integer :: fileid
      integer :: nsym, nbas(*)
#  include "Molcas.fh"
#  include "stdalloc.fh"
#  include "mh5.fh"

      integer :: isym
      integer :: nb, nbast, nbast1, nbast2

      real*8, allocatable :: SAO(:), Scr(:)

      integer :: iRc, iOpt, iComp, iSyLbl
      character(8) :: Label

      integer :: iOff1, iOff2

      integer :: dsetid

      nbast=0
      nbast1=0
      nbast2=0
      do isym=1,nsym
        nb=nbas(isym)
        nbast=nbast+nb
        nbast1=nbast1+(nb*(nb+1))/2
        nbast2=nbast2+nb**2
      end do

*     atomic orbital overlap matrix
      dsetid = mh5_create_dset_real(fileid,
     $        'AO_OVERLAP_MATRIX', 1, [NBAST2])
      call mh5_init_attr(dsetid, 'description',
     $        'Overlap matrix of the atomic orbitals, '//
     $        'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')

      call mma_allocate(SAO,NBAST1)
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,SAO,iSyLbl)
      iOff1 = 0
      iOff2 = 0
      Do iSym = 1,nSym
        nB = nBas(iSym)
        If ( nb.gt.0 ) then
          call mma_allocate(scr,nb*nb)
          Call Square(SAO(1+iOff1),Scr,1,nb,nb)
          call mh5_put_dset(dsetid,
     $            Scr,[nb*nb],[iOff2])
          call mma_deallocate(scr)
        end if
        iOff1 = iOff1 + (nb*nb+nb)/2
        iOff2 = iOff2 + (nb*nb)
      end do
      call mma_deallocate(SAO)

      call mh5_close_dset(dsetid)

      end

      subroutine one2h5_fckint(fileid, nsym, nbas)
*     SVC: read basic molecular information from the RunFile
*     and 1-electron integral file and write it to the HDF5
*     file specified with fileid.
*     This routine does nothing if HDF5 is not supported.
*
*     Attributes:
*       NSYM, IRREPS, POTNUC, NBAS
*     Datasets:
*       BASIS_LABELS, OVLMAT

      implicit none
      integer :: fileid
      integer :: nsym, nbas(*)
#  include "Molcas.fh"
#  include "stdalloc.fh"
#  include "mh5.fh"

      integer :: isym
      integer :: nb, nbast, nbast1, nbast2

      real*8, allocatable :: SAO(:), Scr(:)

      integer :: iRc, iOpt, iComp, iSyLbl
      character(8) :: Label

      integer :: iOff1, iOff2

      integer :: dsetid

      nbast=0
      nbast1=0
      nbast2=0
      do isym=1,nsym
        nb=nbas(isym)
        nbast=nbast+nb
        nbast1=nbast1+(nb*(nb+1))/2
        nbast2=nbast2+nb**2
      end do

*     atomic orbital overlap matrix
      dsetid = mh5_create_dset_real(fileid,
     $        'AO_FOCKINT_MATRIX', 1, [NBAST2])
      call mh5_init_attr(dsetid, 'description',
     $        'Overlap matrix of the atomic orbitals, '//
     $        'arranged as blocks of size [NBAS(i)**2], i=1,#irreps')

      call mma_allocate(SAO,NBAST1)
      iRc=-1
      iOpt=6
      iComp=1
      iSyLbl=1
      Label='FckInt  '
      Call RdOne(iRc,iOpt,Label,iComp,SAO,iSyLbl)
      iOff1 = 0
      iOff2 = 0
      Do iSym = 1,nSym
        nB = nBas(iSym)
        If ( nb.gt.0 ) then
          call mma_allocate(scr,nb*nb)
          Call Square(SAO(1+iOff1),Scr,1,nb,nb)
          call mh5_put_dset(dsetid,
     $            Scr,[nb*nb],[iOff2])
          call mma_deallocate(scr)
        end if
        iOff1 = iOff1 + (nb*nb+nb)/2
        iOff2 = iOff2 + (nb*nb)
      end do
      call mma_deallocate(SAO)

      call mh5_close_dset(dsetid)

      end

      subroutine one2h5_crtmom(fileid, nsym, nbas)
*     SVC: read cartesian moments from the 1-electron integral file
*     and write it to the HDF5 file specified with fileid.
*     This routine does nothing if HDF5 is not supported.
*     FP: also include the origins used for the operators
*
*     Datasets:
*       MLTPL_X, MLTPL_Y, MLTPL_Z
*       MLTPL_XX, MLTPL_YY, MLTPL_ZZ, MLTPL_XY, MLTPL_YZ, MLTPL_XZ
*       MLTPL_ORIG

      implicit none
      integer :: fileid
      integer :: nsym, nbas(*)
#  include "Molcas.fh"
#  include "stdalloc.fh"
#  include "mh5.fh"

      integer :: isym, jsym, msym
      integer :: nb, nbast, nB1, nB2

      real*8, allocatable :: MLTPL(:,:), Scratch(:)
      real*8, dimension(3,3) :: mp_orig

      integer :: iRc, iOpt, iComp, iSyMsk
      character(8) :: Label

      character(1) :: mltpl1_comp(3) = ['X','Y','Z']
      character(2) :: mltpl2_comp(6) = ['XX','XY','XZ','YY','YZ','ZZ']

      integer :: i, j, iOff, jOff, iScrOff, iBas, jBas

      integer :: dsetid

      integer, external :: symmetry

      nbast=0
      do isym=1,nsym
        nb=nbas(isym)
        nbast=nbast+nb
      end do

      mp_orig(:,:) = 0.

      call mma_allocate(MLTPL,NBAST,NBAST)
      call mma_allocate(Scratch,NBAST**2+3)

      do icomp=1,3
      MLTPL=0.0D0
      iRc=-1
      iOpt=4
      iSyMsk=0
      Label='Mltpl  1'
      Call RdOne(iRc,iOpt,Label,iComp,Scratch,iSyMsk)
* iSyMsk tells us which symmetry combination is valid
      iScrOff = 0
      iOff = 0
      Do iSym = 1,nSym
        jOff = 0
        nB1 = nBas(iSym)
        Do jSym = 1,iSym
          mSym = symmetry(iSym,jSym)
          nB2 = nBas(jSym)
          If (IAND(2**(mSym-1),iSyMsk).ne.0) Then
            If (iSym.eq.jSym) Then
              Do j=1,nB2
                jBas = jOff + j
                Do i=1,j
                  iBas = iOff + i
                  MLTPL(iBas,jBas) = Scratch(1+iScrOff)
                  iScrOff = iScrOff + 1
                End Do
              End Do
            Else
              Do j=1,nB2
                jBas = jOff + j
                Do i=1,nB1
                  iBas = iOff + i
                  MLTPL(iBas,jBas) = Scratch(1+iScrOff)
                  iScrOff = iScrOff + 1
                End Do
              End Do
            End If
          End If
          jOff = jOff + nB2
        End Do
        iOff = iOff + nB1
      End Do
      Do j=1,nBasT
        Do i=1,j-1
          MLTPL(j,i)=MLTPL(i,j)
        End Do
      End Do
      dsetid = mh5_create_dset_real(fileid,
     $        'AO_MLTPL_'//mltpl1_comp(icomp), 2, [NBAST,NBAST])
      call mh5_init_attr(dsetid, 'description',
     $        '1st-order multipole matrix of the atomic orbitals, '//
     $        'arranged as matrix of size [NBAST,NBAST]')
      call mh5_put_dset_array_real(dsetid,MLTPL)
      call mh5_close_dset(dsetid)
      end do

      mp_orig(1:3,2) = Scratch(iScrOff+1:iScrOff+3)

      do icomp=1,6
      MLTPL=0.0D0
      iRc=-1
      iOpt=4
      iSyMsk=0
      Label='Mltpl  2'
      Call RdOne(iRc,iOpt,Label,iComp,Scratch,iSyMsk)
* iSyMsk tells us which symmetry combination is valid
      iScrOff = 0
      iOff = 0
      Do iSym = 1,nSym
        jOff = 0
        nB1 = nBas(iSym)
        Do jSym = 1,iSym
          mSym = symmetry(iSym,jSym)
          nB2 = nBas(jSym)
          If (IAND(2**(mSym-1),iSyMsk).ne.0) Then
            If (iSym.eq.jSym) Then
              Do j=1,nB2
                jBas = jOff + j
                Do i=1,j
                  iBas = iOff + i
                  MLTPL(iBas,jBas) = Scratch(1+iScrOff)
                  iScrOff = iScrOff + 1
                End Do
              End Do
            Else
              Do j=1,nB2
                jBas = jOff + j
                Do i=1,nB1
                  iBas = iOff + i
                  MLTPL(iBas,jBas) = Scratch(1+iScrOff)
                  iScrOff = iScrOff + 1
                End Do
              End Do
            End If
          End If
          jOff = jOff + nB2
        End Do
        iOff = iOff + nB1
      End Do
      Do j=1,nBasT
        Do i=1,j-1
          MLTPL(j,i)=MLTPL(i,j)
        End Do
      End Do
      dsetid = mh5_create_dset_real(fileid,
     $        'AO_MLTPL_'//mltpl2_comp(icomp), 2, [NBAST,NBAST])
      call mh5_init_attr(dsetid, 'description',
     $        '2nd-order multipole matrix of the atomic orbitals, '//
     $        'arranged as matrix of size [NBAST,NBAST]')
      call mh5_put_dset_array_real(dsetid,MLTPL)
      call mh5_close_dset(dsetid)
      end do

      mp_orig(1:3,3) = Scratch(iScrOff+1:iScrOff+3)

      call mma_deallocate(MLTPL)
      call mma_deallocate(Scratch)

      dsetid = mh5_create_dset_real(fileid,
     $        'MLTPL_ORIG', 2, [3,3])
      call mh5_init_attr(dsetid, 'description',
     $        'Origin used for the multipole moment operators: '//
     $        'arranged as overlap, dipole, quadrupole')
      call mh5_put_dset_array_real(dsetid,mp_orig,[3,3],[0,0])
      call mh5_close_dset(dsetid)

      end

#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_one2h5_ovlmat()
      End
#endif
