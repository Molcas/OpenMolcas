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
      subroutine getAOs2(lhigh)
      implicit real*8(a-h,o-z)
cbs   get expansions of atomic orbitals in contracted functions
#include "para.fh"
#include "param.fh"
      common /nucleus/ charge,Exp_finite
      integer closedshells(0:LMAX),openshells(0:LMAX)
      call getocc_ao(int(charge),closedshells,openshells)
      do lrun=0,lhigh
      do irun=1,MxcontL
      do jrun=1,MxcontL
      AOcoeffs(jrun,irun,lrun)=0d0
      enddo
      enddo
      enddo
CBS   write(6,*) 'Orbitals for mean-field'
      do lrun=0,lhigh
CBS   write(6,'(A3,I3)') 'L= ',lrun
      do i=1,closedshells(lrun)
      occup(i,lrun)=2.0
      AOcoeffs(i,i,lrun)=1d0
      enddo
      noccorb(lrun)=closedshells(lrun)
      if (openshells(lrun).gt.0) then
      i=closedshells(lrun)+1
      occup(i,lrun)=1d0*DBLE(openshells(lrun))/DBLE(lrun+lrun+1)
      AOcoeffs(i,i,lrun)=1d0
      noccorb(lrun)=i
      endif
      if (noccorb(lrun).gt.0) then
CBS   write(6,'(A,I3)') 'number of orbitals ',noccorb(lrun)
CBS   do iorbital=1,noccorb(lrun)
CBS   write(6,'(A,8F8.4)') 'OCCUPATION: ',(occup(iorbital,lrun),
CBS  *iorbital=1,noccorb(lrun))
CBS   enddo
      endif
      enddo
      return
      end
cbs
      subroutine getocc_ao(icharge,iclosed,iopen)
      implicit real*8(a-h,o-z)
#include "para.fh"
      parameter (ichargemax=103)
      dimension iclocc(0:Lmax_occ,0:ichargemax)
      dimension iopocc(0:Lmax_occ,0:ichargemax)
      character*30 occtxt(0:ichargemax)
      character*35 txt
      data txt/'  SO-integrals were calculated for '/
      dimension iclosed(0:LMAX),iopen(0:LMAX)
*
      data (occtxt(i),i= 0,10) /
     *'dummy atom (no integrals)     ',
     *' H: no mean-field             ',
     *'He: 1s^2                      ',
     *'Li: [He]2s^1                  ',
     *'Be: [He]2s^2                  ',
     *' B: [He]2s^2 2p^1             ',
     *' C: [He]2s^2 2p^2             ',
     *' N: [He]2s^2 2p^3             ',
     *' O: [He]2s^2 2p^4             ',
     *' F: [He]2s^2 2p^5             ',
     *'Ne: [He]2s^2 2p^6             '/
      data (occtxt(i),i= 11,20) /
     *'Na: [Ne]3s^1                  ',
     *'Mg: [Ne]3s^2                  ',
     *'Al: [Ne]3s^2 3p^1             ',
     *'Si: [Ne]3s^2 3p^2             ',
     *' P: [Ne]3s^2 3p^3             ',
     *' S: [Ne]3s^2 3p^4             ',
     *'Cl: [Ne]3s^2 3p^5             ',
     *'Ar: [Ne]3s^2 3p^6             ',
     *' K: [Ar]4s^1                  ',
     *'Ca: [Ar]4s^2                  '/
      data (occtxt(i),i= 21,30) /
     *'Sc: [Ar]4s^2 3d^1             ',
     *'Ti: [Ar]4s^2 3d^2             ',
     *' V: [Ar]4s^2 3d^3             ',
     *'Cr: [Ar]4s^2 3d^4             ',
     *'Mn: [Ar]4s^2 3d^5             ',
     *'Fe: [Ar]4s^2 3d^6             ',
     *'Co: [Ar]4s^2 3d^7             ',
     *'Ni: [Ar]4s^2 3d^8             ',
     *'Cu: [Ar]4s^1 3d^10            ',
     *'Zn: [Ar]4s^2 3d^10            '/
      data (occtxt(i),i= 31,40) /
     *'Ga: [Ar]4s^2 3d^10 4p^1       ',
     *'Ge: [Ar]4s^2 3d^10 4p^2       ',
     *'As: [Ar]4s^2 3d^10 4p^3       ',
     *'Se: [Ar]4s^2 3d^10 4p^4       ',
     *'Br: [Ar]4s^2 3d^10 4p^5       ',
     *'Kr: [Ar]4s^2 3d^10 4p^6       ',
     *'Rb: [Kr]5s^1                  ',
     *'Sr: [Kr]5s^2                  ',
     *' Y: [Kr]5s^2 4d^1             ',
     *'Zr: [Kr]5s^2 4d^2             '/
      data (occtxt(i),i= 41,50) /
     *'Nb: [Kr]5s^2 4d^3             ',
     *'Mo: [Kr]5s^2 4d^4             ',
     *'Tc: [Kr]5s^2 4d^5             ',
     *'Ru: [Kr]5s^2 4d^6             ',
     *'Rh: [Kr]5s^2 4d^7             ',
     *'Pd: [Kr]5s^2 4d^8             ',
     *'Ag: [Kr]5s^1 4d^10            ',
     *'Cd: [Kr]5s^2 4d^10            ',
     *'In: [Kr]5s^2 4d^10 5p^1       ',
     *'Sn: [Kr]5s^2 4d^10 5p^2       '/
      data (occtxt(i),i= 51,60) /
     *'Sb: [Kr]5s^2 4d^10 5p^3       ',
     *'Te: [Kr]5s^2 4d^10 5p^4       ',
     *' I: [Kr]5s^2 4d^10 5p^5       ',
     *'Xe: [Kr]5s^2 4d^10 5p^6       ',
     *'Cs: [Xe]6s^1                  ',
     *'Ba: [Xe]6s^2                  ',
     *'La: [Xe]6s^2 5d^1             ',
     *'Ce: [Xe]6s^2 4f^2             ',
     *'Pr: [Xe]6s^2 4f^3             ',
     *'Nd: [Xe]6s^2 4f^4             '/
      data (occtxt(i),i= 61,70) /
     *'Pm: [Xe]6s^2 4f^5             ',
     *'Sm: [Xe]6s^2 4f^6             ',
     *'Eu: [Xe]6s^2 4f^7             ',
     *'Gd: [Xe]6s^2 4f^8             ',
     *'Tb: [Xe]6s^2 4f^9             ',
     *'Dy: [Xe]6s^2 4f^10            ',
     *'Ho: [Xe]6s^2 4f^11            ',
     *'Er: [Xe]6s^2 4f^12            ',
     *'Tm: [Xe]6s^2 4f^13            ',
     *'Yb: [Xe]6s^2 4f^14            '/
      data (occtxt(i),i= 71,80) /
     *'Lu: [Xe+4f^14]6s^2 5d^1       ',
     *'Hf: [Xe+4f^14]6s^2 5d^2       ',
     *'Ta: [Xe+4f^14]6s^2 5d^3       ',
     *' W: [Xe+4f^14]6s^2 5d^4       ',
     *'Re: [Xe+4f^14]6s^2 5d^5       ',
     *'Os: [Xe+4f^14]6s^2 5d^6       ',
     *'Ir: [Xe+4f^14]6s^2 5d^7       ',
     *'Pt: [Xe+4f^14]6s^1 5d^9       ',
     *'Au: [Xe+4f^14]6s^1 5d^10      ',
     *'Hg: [Xe+4f^14]6s^2 5d^10      '/
      data (occtxt(i),i= 81,90) /
     *'Tl: [Xe+4f^14+5d^10]6s^2 6p^1 ',
     *'Pb: [Xe+4f^14+5d^10]6s^2 6p^2 ',
     *'Bi: [Xe+4f^14+5d^10]6s^2 6p^3 ',
     *'Po: [Xe+4f^14+5d^10]6s^2 6p^4 ',
     *'At: [Xe+4f^14+5d^10]6s^2 6p^5 ',
     *'Rn: [Xe+4f^14+5d^10]6s^2 6p^6 ',
     *'Fr: [Rn]7s^1                  ',
     *'Ra: [Rn]7s^2                  ',
     *'Ac: [Rn]7s^2 6d^1             ',
     *'Th: [Rn]7s^2 6d^2             '/
      data (occtxt(i),i= 91,iChargeMax) /
     *'Pa: [Rn]7s^2 6d^1 5f^2        ',
     *' U: [Rn]7s^2 6d^1 5f^3        ',
     *'Np: [Rn]7s^2 6d^1 5f^4        ',
     *'Pu: [Rn]7s^2 6d^0 5f^6        ',
     *'Am: [Rn]7s^2 6d^0 5f^7        ',
     *'Cm: [Rn]7s^2 6d^0 5f^8        ',
     *'Bk: [Rn]7s^2 6d^0 5f^9        ',
     *'Cf: [Rn]7s^2 6d^0 5f^10       ',
     *'Es: [Rn]7s^2 6d^0 5f^11       ',
     *'Fm: [Rn]7s^2 6d^0 5f^12       ',
     *'Md: [Rn]7s^2 6d^0 5f^13       ',
     *'No: [Rn]7s^2 6d^0 5f^14       ',
     *'Lr: [Rn]7s^2 6d^1 5f^14       '/
*
      data ((iclocc(i,j),i=0,LMAX_occ),j=0,10) /
     & 0 , 0, 0, 0,       !0
     & 0 , 0, 0, 0,       !1
     & 1 , 0, 0, 0,       !2
     & 1 , 0, 0, 0,       !3
     & 2 , 0, 0, 0,       !4
     & 2 , 0, 0, 0,       !5
     & 2 , 0, 0, 0,       !6
     & 2 , 0, 0, 0,       !7
     & 2 , 0, 0, 0,       !8
     & 2 , 0, 0, 0,       !9
     & 2 , 1, 0, 0/       !10
      data ((iclocc(i,j),i=0,LMAX_occ),j=11,20) /
     & 2 , 1, 0, 0,       !11
     & 3 , 1, 0, 0,       !12
     & 3 , 1, 0, 0,       !13
     & 3 , 1, 0, 0,       !14
     & 3 , 1, 0, 0,       !15
     & 3 , 1, 0, 0,       !16
     & 3 , 1, 0, 0,       !17
     & 3 , 2, 0, 0,       !18
     & 3 , 2, 0, 0,       !19
     & 4 , 2, 0, 0/       !20
      data ((iclocc(i,j),i=0,LMAX_occ),j=21,30) /
     & 4 , 2, 0, 0,       !21
     & 4 , 2, 0, 0,       !22
     & 4 , 2, 0, 0,       !23
     & 4 , 2, 0, 0,       !24
     & 4 , 2, 0, 0,       !25
     & 4 , 2, 0, 0,       !26
     & 4 , 2, 0, 0,       !27
     & 4 , 2, 0, 0,       !28
     & 3 , 2, 1, 0,       !29
     & 4 , 2, 1, 0/       !30
      data ((iclocc(i,j),i=0,LMAX_occ),j=31,40) /
     & 4 , 2, 1, 0,       !31
     & 4 , 2, 1, 0,       !32
     & 4 , 2, 1, 0,       !33
     & 4 , 2, 1, 0,       !34
     & 4 , 2, 1, 0,       !35
     & 4 , 3, 1, 0,       !36
     & 4 , 3, 1, 0,       !37
     & 5 , 3, 1, 0,       !38
     & 5 , 3, 1, 0,       !39
     & 5 , 3, 1, 0/       !40
      data ((iclocc(i,j),i=0,LMAX_occ),j=41,50) /
     & 5 , 3, 1, 0,       !41
     & 5 , 3, 1, 0,       !42
     & 5 , 3, 1, 0,       !43
     & 5 , 3, 1, 0,       !44
     & 5 , 3, 1, 0,       !45
     & 5 , 3, 1, 0,       !46
     & 4 , 3, 2, 0,       !47
     & 5 , 3, 2, 0,       !48
     & 5 , 3, 2, 0,       !49
     & 5 , 3, 2, 0/       !50
      data ((iclocc(i,j),i=0,LMAX_occ),j=51,60) /
     & 5 , 3, 2, 0,       !51
     & 5 , 3, 2, 0,       !52
     & 5 , 3, 2, 0,       !53
     & 5 , 4, 2, 0,       !54
     & 5 , 4, 2, 0,       !55
     & 6 , 4, 2, 0,       !56
     & 6 , 4, 2, 0,       !57
     & 6 , 4, 2, 0,       !58
     & 6 , 4, 2, 0,       !59
     & 6 , 4, 2, 0/       !60
      data ((iclocc(i,j),i=0,LMAX_occ),j=61,70) /
     & 6 , 4, 2, 0,       !61
     & 6 , 4, 2, 0,       !62
     & 6 , 4, 2, 0,       !63
     & 6 , 4, 2, 0,       !64
     & 6 , 4, 2, 0,       !65
     & 6 , 4, 2, 0,       !66
     & 6 , 4, 2, 0,       !67
     & 6 , 4, 2, 0,       !68
     & 6 , 4, 2, 0,       !69
     & 6 , 4, 2, 1/       !70
      data ((iclocc(i,j),i=0,LMAX_occ),j=71,80) /
     & 6 , 4, 2, 1,       !71
     & 6 , 4, 2, 1,       !72
     & 6 , 4, 2, 1,       !73
     & 6 , 4, 2, 1,       !74
     & 6 , 4, 2, 1,       !75
     & 6 , 4, 2, 1,       !76
     & 6 , 4, 2, 1,       !77
     & 5 , 4, 2, 1,       !78
     & 5 , 4, 3, 1,       !79
     & 6 , 4, 3, 1/       !80
      data ((iclocc(i,j),i=0,LMAX_occ),j=81,90) /
     & 6 , 4, 3, 1,       !81
     & 6 , 4, 3, 1,       !82
     & 6 , 4, 3, 1,       !83
     & 6 , 4, 3, 1,       !84
     & 6 , 4, 3, 1,       !85
     & 6 , 5, 3, 1,       !86
     & 6 , 5, 3, 1,       !87
     & 7 , 5, 3, 1,       !88
     & 7 , 5, 3, 1,       !89
     & 7 , 5, 3, 1/       !90
      data ((iclocc(i,j),i=0,LMAX_occ),j=91,ichargemax) /
     & 7 , 5, 3, 1,       !91
     & 7 , 5, 3, 1,       !92
     & 7 , 5, 3, 1,       !93
     & 7 , 5, 3, 1,       !94
     & 7 , 5, 3, 1,       !95
     & 7 , 5, 3, 1,       !96
     & 7 , 5, 3, 1,       !97
     & 7 , 5, 3, 1,       !98
     & 7 , 5, 3, 1,       !99
     & 7 , 5, 3, 1,       !100
     & 7 , 5, 3, 1,       !101
     & 7 , 5, 3, 2,       !102
     & 7 , 5, 3, 2/       !103
cbs
      data ((iopocc(i,j),i=0,LMAX_occ),j=0,10) /
     & 0 , 0, 0, 0,    !0
     & 0 , 0, 0, 0,    ! 1
     & 0 , 0, 0, 0,    ! 2
     & 1 , 0, 0, 0,    ! 3
     & 0 , 0, 0, 0,    ! 4
     & 0 , 1, 0, 0,    ! 5
     & 0 , 2, 0, 0,    ! 6
     & 0 , 3, 0, 0,    ! 7
     & 0 , 4, 0, 0,    ! 8
     & 0 , 5, 0, 0,    ! 9
     & 0 , 0, 0, 0/    ! 10
      data ((iopocc(i,j),i=0,LMAX_occ),j=11,20) /
     & 1 , 0, 0, 0,    ! 11
     & 0 , 0, 0, 0,    ! 12
     & 0 , 1, 0, 0,    ! 13
     & 0 , 2, 0, 0,    ! 14
     & 0 , 3, 0, 0,    ! 15
     & 0 , 4, 0, 0,    ! 16
     & 0 , 5, 0, 0,    ! 17
     & 0 , 0, 0, 0,    ! 18
     & 1 , 0, 0, 0,    ! 19
     & 0 , 0, 0, 0/    ! 20
      data ((iopocc(i,j),i=0,LMAX_occ),j=21,30) /
     & 0 , 0, 1, 0,    ! 21
     & 0 , 0, 2, 0,    ! 22
     & 0 , 0, 3, 0,    ! 23
     & 0 , 0, 4, 0,    ! 24
     & 0 , 0, 5, 0,    ! 25
     & 0 , 0, 6, 0,    ! 26
     & 0 , 0, 7, 0,    ! 27
     & 0 , 0, 8, 0,    ! 28
     & 1 , 0, 0, 0,    ! 29
     & 0 , 0, 0, 0/    ! 30
      data ((iopocc(i,j),i=0,LMAX_occ),j=31,40) /
     & 0 , 1, 0, 0,    ! 31
     & 0 , 2, 0, 0,    ! 32
     & 0 , 3, 0, 0,    ! 33
     & 0 , 4, 0, 0,    ! 34
     & 0 , 5, 0, 0,    ! 35
     & 0 , 0, 0, 0,    ! 36
     & 1 , 0, 0, 0,    ! 37
     & 0 , 0, 0, 0,    ! 38
     & 0 , 0, 1, 0,    ! 39
     & 0 , 0, 2, 0/    ! 40
      data ((iopocc(i,j),i=0,LMAX_occ),j=41,50) /
     & 0 , 0, 3, 0,    ! 41
     & 0 , 0, 4, 0,    ! 42
     & 0 , 0, 5, 0,    ! 43
     & 0 , 0, 6, 0,    ! 44
     & 0 , 0, 7, 0,    ! 45
     & 0 , 0, 8, 0,    ! 46
     & 1 , 0, 0, 0,    ! 47
     & 0 , 0, 0, 0,    ! 48
     & 0 , 1, 0, 0,    ! 49
     & 0 , 2, 0, 0/    ! 50
      data ((iopocc(i,j),i=0,LMAX_occ),j=51,60) /
     & 0 , 3, 0, 0,    ! 51
     & 0 , 4, 0, 0,    ! 52
     & 0 , 5, 0, 0,    ! 53
     & 0 , 0, 0, 0,    ! 54
     & 1 , 0, 0, 0,    ! 55
     & 0 , 0, 0, 0,    ! 56
     & 0 , 0, 1, 0,    ! 57
     & 0 , 0, 0, 2,    ! 58
     & 0 , 0, 0, 3,    ! 59
     & 0 , 0, 0, 4/    ! 60
      data ((iopocc(i,j),i=0,LMAX_occ),j=61,70) /
     & 0 , 0, 0, 5,    ! 61
     & 0 , 0, 0, 6,    ! 62
     & 0 , 0, 0, 7,    ! 63
     & 0 , 0, 0, 8,    ! 64
     & 0 , 0, 0, 9,    ! 65
     & 0 , 0, 0, 10,    ! 66
     & 0 , 0, 0, 11,    ! 67
     & 0 , 0, 0, 12,    ! 68
     & 0 , 0, 0, 13,    ! 69
     & 0 , 0, 0,  0/    ! 70
      data ((iopocc(i,j),i=0,LMAX_occ),j=71,80) /
     & 0 , 0, 1, 0,    ! 71
     & 0 , 0, 2, 0,    ! 72
     & 0 , 0, 3, 0,    ! 73
     & 0 , 0, 4, 0,    ! 74
     & 0 , 0, 5, 0,    ! 75
     & 0 , 0, 6, 0,    ! 76
     & 0 , 0, 7, 0,    ! 77
     & 1 , 0, 9, 0,    ! 78
     & 1 , 0, 0, 0,    ! 79
     & 0 , 0, 0, 0/    ! 80
      data ((iopocc(i,j),i=0,LMAX_occ),j=81,90) /
     & 0 , 1, 0, 0,    ! 81
     & 0 , 2, 0, 0,    ! 82
     & 0 , 3, 0, 0,    ! 83
     & 0 , 4, 0, 0,    ! 84
     & 0 , 5, 0, 0,    ! 85
     & 0 , 0, 0, 0,    ! 86
     & 1 , 0, 0, 0,    ! 87
     & 0 , 0, 0, 0,    ! 88
     & 0 , 0, 1, 0,    ! 89
     & 0 , 0, 2, 0/    ! 90
      data ((iopocc(i,j),i=0,LMAX_occ),j=91,ichargemax) /
     & 0 , 0, 1, 2,    ! 91
     & 0 , 0, 1, 3,    ! 92
     & 0 , 0, 1, 4,    ! 93
     & 0 , 0, 0, 6,    ! 94
     & 0 , 0, 0, 7,    ! 95
     & 0 , 0, 0, 8,    ! 96
     & 0 , 0, 0, 9,    ! 97
     & 0 , 0, 0, 10,   ! 98
     & 0 , 0, 0, 11,   ! 99
     & 0 , 0, 0, 12,   ! 100
     & 0 , 0, 0, 13,   ! 101
     & 0 , 0, 0, 0,    ! 102
     & 0 , 0, 1, 0/    ! 103
cbs
      if (icharge.gt.ichargemax) then
         write(6,*) 'occupations not implemented'
         Call Abend()
      endif
*
      iPL=iPrintLevel(-1)
      If (iPL.ge.3) write(6,'(A35,A30)') txt,occtxt(icharge)
*
      do irun=0,min(lmax,lmax_occ)
         iclosed(irun)=iclocc(irun,icharge)
         iopen(irun)=iopocc(irun,icharge)
      end do
      do irun=min(lmax,lmax_occ)+1,lmax
         iclosed(irun)=0
         iopen(irun)=0
      end do
      return
      end
