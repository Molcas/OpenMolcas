!!  dmrg-interface-utils: interface to the Maquis DMRG program for various
!!                        quantum-chemistry program packages.
!!  Copyright 2013-2018 Leon Freitag, Erik Hedegaard, Sebastian Keller,
!!                      Stefan Knecht, Yingjin Ma, Christopher Stein
!!                      and Markus Reiher
!!                      Laboratory for Physical Chemistry, ETH Zurich
!!
!!  dmrg-interface-utils is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!!
!!  dmrg-interface-utils is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with dmrg-interface-utils. If not, see <http://www.gnu.org/licenses/>.

module qcmaquis_interface_main

 !> stefan: interface to DMRG tasks - provide functionalities based on Fortran-->python interface to
 !>         exchange data (1-e,2-e integrals, X-particle reduced (transition) density matrices)

 use qcmaquis_interface_cfg
 use qcmaquis_interface_utility_routines
 use qcmaquis_interface_environment, only: set_dmrg_runtime_environment
 use qcmaquis_interface_measurements

implicit none

  !> public
  public qcmaquis_interface_ctl
  public file_name_generator

contains

  subroutine qcmaquis_interface_ctl(                                  &
                                    task,                             &
                                    x1,                               &
                                    x2,                               &
                                    x3,                               &
                                    x4,                               &
                                    energy,                           &
                                    ndim,                             &
                                    mdim,                             &
                                    odim,                             &
                                    pdim,                             &
                                    state,                            &
                                    stateL,                           &
                                    msproj,                           &
                                    msprojL,                          &
                                    multiplet,                        &
                                    multipletL,                       &
                                    rdm1,                             &
                                    rdm2,                             &
                                    rdm3,                             &
                                    rdm4,                             &
                                    Key_CION,                         &
                                    IterSCF,                          &
                                    checkpoint1,                      &
                                    checkpoint2                       &
                                   )

      character(len=8),                        intent(in)    :: task
      real*8 , optional, dimension(*),         intent(inout) :: x1
      real*8 , optional, dimension(*),         intent(inout) :: x2
      real*8 , optional, dimension(*),         intent(inout) :: x3
      real*8 , optional, dimension(*),         intent(inout) :: x4
      real*8 , optional,                       intent(inout) :: energy
      integer(kind=8), optional,               intent(in)    :: ndim
      integer(kind=8), optional,               intent(in)    :: mdim
      integer(kind=8), optional,               intent(in)    :: odim
      integer(kind=8), optional,               intent(in)    :: pdim
      integer(kind=8), optional,               intent(in)    :: state
      integer(kind=8), optional,               intent(in)    :: stateL
      integer(kind=8), optional,               intent(in)    :: msproj
      integer(kind=8), optional,               intent(in)    :: msprojL
      integer(kind=8), optional,               intent(in)    :: multiplet
      integer(kind=8), optional,               intent(in)    :: multipletL
      logical, optional,                       intent(in)    :: rdm1
      logical, optional,                       intent(in)    :: rdm2
      logical, optional,                       intent(in)    :: rdm3
      logical, optional,                       intent(in)    :: rdm4
      logical, optional,                       intent(in)    :: Key_CION
      integer(kind=8), optional,               intent(in)    :: iterSCF
      character(len=*),optional,               intent(in)    :: checkpoint1
      character(len=*),optional,               intent(in)    :: checkpoint2

        ! select task
        select case(trim(task))
          case('fci dump')
            call dmrg_task_fcidump(                                 &
                                   oneint    = x1,                  &
                                   twoint    = x2,                  &
                                   corenergy = energy               &
                                  )
          case('tra dump')
            call dmrg_task_tradump(                                 &
                                   tra       = x1,                  &
                                   no        = ndim,                &
                                   ni        = mdim,                &
                                   jorb      = state,               &
                                   tirrep    = stateL               &
                                  )
          case('overlap ')
            call dmrg_task_overlap(                                 &
                                   iroot     = state,               &
                                   jroot     = stateL,              &
                                   overlap   = energy,              &
                                 checkpoint1 = checkpoint1,         &
                                 checkpoint2 = checkpoint2          &
                          )
          case('overlapR')
            call dmrg_task_overlapR(                                &
                                    jroot     = stateL,             &
                                    overlap   = energy              &
                          )
          case('overlapU')
          ! Leon: 2U1 overlap. Only with checkpoint names so far.
          ! TODO: extend for non-checkpoint name usage
          ! (with iroot/jroot etc.)
            call dmrg_task_overlapU(                                &
                                   energy,                          &
                                   checkpoint1,                     &
                                   checkpoint2                      &
                          )
          case('imp spdX')
            call dmrg_task_import_spd(                              &
                                      spd1  = x1,                   &
                                      nrdm  = ndim,                 &
                                      iroot = state,                &
                                      sp1   = rdm1                  &
                                     )
          case('imp rdmY')
            call dmrg_task_import_rdmY(                             &
                                       dv         = x1,             &
                                       pv         = x2,             &
                                       tv         = x3,             &
                                       nrdm1      = ndim,           &
                                       nrdm2      = mdim,           &
                                       nrdm3      = pdim,           &
                                       iroot      = state,          &
                                       jroot      = stateL,         &
                                       msproj     = msproj,         &
                                       msprojL    = msprojL,        &
                                       multiplet  = multiplet,      &
                                       multipletL = multipletL,     &
                                       rdm1       = rdm1,           &
                                       rdm2       = rdm2,           &
                                       rdm3       = rdm3,           &
                                      checkpoint1 = checkpoint1,    &
                                      checkpoint2 = checkpoint2     &
                                     )
          case('imp rdmX')
            call dmrg_task_import_rdm(                              &
                                      dv    = x1,                   &
                                      pv    = x2,                   &
                                      tv    = x3,                   &
                                      fv    = x4,                   &
                                      nrdm  = ndim,                 &
                                      mrdm  = mdim,                 &
                                      ordm  = odim,                 &
                                      prdm  = pdim,                 &
                                      iroot = state,                &
                                      rdm1  = rdm1,                 &
                                      rdm2  = rdm2,                 &
                                      rdm3  = rdm3,                 &
                                      rdm4  = rdm4                  &
                                     )
          case('run DMRG')
            call dmrg_task_run_dmrg  (                              &
                               Key_DMRGonly = Key_CION,             &
                                       Iter = iterSCF    )
          case('MPS crot')
            call dmrg_task_MPScrot(                                 &
                                   iroot     = state,               &
                                 checkpoint  = checkpoint1          &
                          )
          case('MPS back')
            call dmrg_task_MPSbackup(                               &
                                     iroot     = state,             &
                                     switch    = stateL,            &
                                   checkpoint  = checkpoint1        &
                          )
! write the QCMaquis 4rdm or 3tdm evaluation to disk
          case('w 4rdmin')
            call dmrg_task_prepare_hirdm_template(                  &
                                     task      = '4rdm',            &
                                     iroot     = state,             &
                                     run_qcmaquis = rdm4,           &
                                    compress_Mmax = msproj          &
                          )
          case('w 3tdmin')
            call dmrg_task_prepare_hirdm_template(                  &
                                     task      = '3tdm',            &
                                     iroot     = state,             &
                                     jroot     = stateL,            &
                                     run_qcmaquis = rdm3            &
                          )
! Compress an mps
          case('MPS comp')
            call dmrg_task_mps_compress(                            &
                                     iroot     = state,             &
                                     Mmax      = stateL             &
                          )
          case default
            print *, 'Unknown task --> ',trim(task),' <-- in qcmaquis_interface_ctl ...'
            stop -911
        end select

  end subroutine qcmaquis_interface_ctl
!**********************************************************************

      subroutine dmrg_task_fcidump(oneint,twoint,corenergy)
!     -----------------------------------------------------------------
!
!     purpose: write one- and two-electron integrals for a given active
!              space to disk in a formatted file.
!
!     filename: FCIDUMP
!     -----------------------------------------------------------------
!
      real*8 , intent(in), dimension(*) :: oneint
      real*8 , intent(in), dimension(*) :: twoint
      real*8 , intent(in)               :: corenergy
!     -----------------------------------------------------------------
      integer                           :: isym, ksym, lsym, jsym, ijsym, klsym
      integer                           :: i, j, k, l, ij, kl, ijkl
      integer                           :: ndummy,lenorbstring
      integer                           :: norbtot, noccend, noccendi, noccendj
      integer                           :: offset, offseti, offsetj, offsetk, offsetl
      character(len=30)                 :: form1
      character(len=30)                 :: form2
      character(len=30)                 :: form3
      character(len=5000)               :: orbstring
      integer, parameter                :: fcidump = 99
      real*8 , parameter                :: threshold = 1.0d-16
!     -----------------------------------------------------------------

!
!     define printing format | Using the same G21.12 as molpro
      form1="(G21.12, 4X, I6, I6, I6, I6)"
      form2="(A11,I3,A7,I2,A5,I2,A1)"
      form3="(G21.12, 4X, I6, I6, I6, I6)"

!     calculate total number of active orbitals
      norbtot = 0
      do ksym = 1, dmrg_symmetry%nirrep
        norbtot = norbtot + dmrg_orbital_space%nash(ksym)
      end do

!     Print header of FCIDUMP file using the MOLPRO format
      open(fcidump,file='FCIDUMP',status='replace',form='formatted',    &
           action='readwrite',position='rewind')

      write(fcidump,form2) ' &FCI NORB=', norbtot , ',NELEC=',          &
      dmrg_state%nactel, ',MS2=', dmrg_state%ms2, ','

      write(fcidump,"(A)",advance='no') '  ORBSYM='
      do isym = 1, dmrg_symmetry%nirrep
        if(dmrg_orbital_space%nash(isym) /= 0)then
          do i = 1, dmrg_orbital_space%nash(isym)
            write(fcidump,"(I1,A1)",advance='no') isym,','
          end do
        end if
      end do
      !> remove trailing ',' from FCIDUMP orbital symmetry string
      backspace(fcidump); read(fcidump,'(a)') orbstring; backspace(fcidump)
      lenorbstring = len_trim(orbstring); write(fcidump,'(a)',advance='no') orbstring(1:lenorbstring-1)

      write(fcidump,*)
      write(fcidump,"(A7,I1)") '  ISYM=',dmrg_state%irefsm
      write(fcidump,"(A5)") ' &END'

      if(dmrg_state%nactel > 1)then

        select case(dmrg_host_program_settings%dmrg_host_program)

        case ('dirac  ')

          if(dmrg_symmetry%nirrep > 1) stop 'symmetry handling not implemented for Dirac interface'

          print *, 'dump integrals in Dirac interface mode'
          !> next two-electron integrals
          !> integrals are sorted in (IJ|KL) order
          offset = 0

!         IJ KL
!         sum over all irreducible representations
          do isym = 1, dmrg_symmetry%nirrep
            if(dmrg_orbital_space%nash(isym) == 0) cycle
            do i = 1, dmrg_orbital_space%nash(isym)
              do j = 1, i
                do k = 1, i
                  do l = 1, k

                    !> check for redundant integrals
                    if(k == i .and. l > j) cycle
                    offset = (l-1)*(dmrg_orbital_space%nash(isym)**3)+(k-1)*(dmrg_orbital_space%nash(isym)**2)+&
                             (j-1)*dmrg_orbital_space%nash(isym)+i
                    !print *, 'offset for (ij|kl) ',i,j,k,l,' is ==> ',offset
                    if(dabs(twoint(offset)) < threshold) cycle
                    write(fcidump,form1) twoint(offset), i, j, k, l

                  end do ! do l
                end do ! do k
              end do ! do j
            end do ! do i
          end do ! do isym

        case default
          !> next two-electron integrals
          !> integrals are sorted in (KL|IJ) order
          offset = 0

!         KL IJ
!         define orbital offset for K
          offsetk = 0

!         sum over all irreducible representations
          do ksym = 1, dmrg_symmetry%nirrep

            if(dmrg_orbital_space%nash(ksym) == 0) cycle

            do k = 1, dmrg_orbital_space%nash(ksym)

              offsetl = 0 !define orbital offset for L

              do lsym = 1, ksym !restrict summation to prevent double counting

                if(dmrg_orbital_space%nash(lsym) == 0)cycle

!               set upper summation bound for orbital index (prevent double counting):
!               if not the same irrep l goes from 1 to number of orbitals
                if(ksym == lsym)then
                  noccend = k
                else
                  noccend = dmrg_orbital_space%nash(lsym)
                end if

                do l = 1, noccend

!                 orbital offset for I
                  offseti = 0

!                 restrict summation to prevent double counting for both irrep ISYM and orbital indices i
                  do isym = 1, ksym

                    if(dmrg_orbital_space%nash(isym) == 0)cycle

                    if(isym == ksym)then
                      noccendi = k
                    else
                      noccendi = dmrg_orbital_space%nash(isym)
                    end if

                    do i = 1, noccendi

!                     set orbital offset J
                      offsetj = 0

!                     double counting issue: irrep of J must be smaller or equal to irrep of I
                      do jsym = 1, isym
!                       fetch integrals which are nonzero by symmetry
!                       two cases have to be distinguished: IJ|KL  and IK|JL
                        if(dmrg_orbital_space%nash(jsym) == 0)cycle

                        if(isym == jsym .and. ksym == lsym)then ! first case
                          ijsym   = dmrg_symmetry%multiplication_table(isym,jsym)
                          klsym   = dmrg_symmetry%multiplication_table(ksym,lsym)
                        else                                    ! second case
                          ijsym   = dmrg_symmetry%multiplication_table(isym,ksym)
                          klsym   = dmrg_symmetry%multiplication_table(jsym,lsym)
                        end if

                        if(dmrg_symmetry%multiplication_table(ijsym,klsym) == 1)then
                          offset = (offsetk+k)*(offsetk+k-1)/2*(norbtot*(norbtot+1)/2)+ &
                                   (offsetl+l-1)*(norbtot*(norbtot+1)/2)+1

!                         prevent double counting of symmetry redundant indices: set upper summation index
!                         if IJKL in same irrep, restrict j to at most i
                          if(jsym == isym .and. ksym == lsym)then !.and.ISYM == KSYM
                            noccendj = i
!                           if LJ in irrep1 and IK in irrep2
                          else if(lsym == jsym .and. ksym == isym .and. lsym /= isym)then
!                           second restriction to prevent double counting, J<=L in KL IJ
                            if(k == i)then
                              noccendj = l
                            else ! otherwise all J are needed
                              noccendj = dmrg_orbital_space%nash(jsym)
                            end if
                          else
                            noccendj = dmrg_orbital_space%nash(jsym)
                          end if
                          offset = offset + (offseti+i)*(offseti+i-1)/2+offsetj
                          do j=1,noccendj,1
!                           check for redundant integrals
                            if(JSYM == ISYM .and. KSYM == lsym .and. ISYM == KSYM)then
                              if(dmrg_host_program_settings%dmrg_host_program == 'dalton ')then
                                if(k == i .and. l == i .and. j < i)then
                                   offset = offset + 1
                                   cycle
                                end if
                                if(k == i .and. j < l)then
                                   offset = offset + 1
                                   cycle
                                end if
                              else
                                if(k == i .and. l == i .and. j < i) cycle
                                if(k == i .and. j < l) cycle
                              end if
                            end if
                            if(dmrg_host_program_settings%dmrg_host_program == 'dalton ')then
                              ijkl   = offset
                              offset = offset + 1
                             !write(6,*) 'org. offset',ijkl
                             !write(6,*) 'indices... ',i+offseti,j+offsetj,k+offsetk,l+offsetl
                            else
                              ij   = max(i+offseti,j+offsetj)*(max(i+offseti,j+offsetj)-1)/2+min(i+offseti,j+offsetj)
                              kl   = max(k+offsetk,l+offsetl)*(max(k+offsetk,l+offsetl)-1)/2+min(k+offsetk,l+offsetl)
                              ijkl = max(ij,kl)*(max(ij,kl)-1)/2+min(ij,kl)
                            end if

!                           write(6,*) 'indices... ',i+offseti,j+offsetj,k+offsetk,l+offsetl
!                           write(6,*) 'offset... ij, kl, ijkl',ij,kl,ijkl
                            if(dabs(twoint(ijkl)) < threshold)then
!                              write(fcidump,form1) 0.000000000000E-15, i+offseti, j+offsetj, k+offsetk, l+offsetl
                              cycle
                            else
                              write(fcidump,form1) twoint(ijkl), i+offseti, j+offsetj, k+offsetk, l+offsetl
                            end if
                          end do ! do j
                        end if ! if dmrg_symmetry%multiplication_table(ijsym,klsym) == 1

                        offsetj = offsetj + dmrg_orbital_space%nash(jsym) !update orbital offset J
                      end do ! do jsym
                    end do ! do i
                    offseti = offseti + dmrg_orbital_space%nash(isym) !update orbital offset I
                  end do ! do isym
                end do ! do l
                offsetl = offsetl + dmrg_orbital_space%nash(lsym) !update orbital offset L
              end do ! do lsym
            end do ! do k
            offsetk = offsetk + dmrg_orbital_space%nash(ksym) !update orbital offset K
          end do ! do ksym

        end select

      end if ! dmrg_state%nactel > 1

!     next step: one-electron integrals
      offset = 0
!     keep track of dummy indices to be ignored in one-electron
!     integrals because the symmetry ISYM < actual ISYM
      ndummy = 0
!
      do isym = 1, dmrg_symmetry%nirrep

        if(dmrg_orbital_space%nash(isym) == 0) cycle

        offset = offset + ndummy

        do i = 1, dmrg_orbital_space%nash(isym) ! only loop through same irrep

          do j = 1, i

            offset = offset + 1

            if(j == i)then
              if (dmrg_host_program_settings%dmrg_host_program(1:7) == 'molcas ')then
                 if(dabs(oneint(offset)-(corenergy/dble(dmrg_state%nactel))) < threshold)then
                   cycle
                 else
                   write(fcidump,form1) oneint(offset)-(corenergy/dble(dmrg_state%nactel)), & ! subtract scaled inactive energy from diagonal elements
                                                        i+ndummy, j+ndummy,0, 0
                 end if
              else
                 if(dabs(oneint(offset)) < threshold)then
                   cycle
                 else
                   write(fcidump,form1) oneint(offset), i+ndummy, j+ndummy,0, 0
                 end if
              end if
            else
              if(dabs(oneint(offset)) < threshold)then
                cycle
              else
                write(fcidump,form1) oneint(offset), i+ndummy, j+ndummy,0, 0
              end if
            end if
          end do
!         add orbital offset if irrep changes
          if(i < dmrg_orbital_space%nash(isym)) offset = offset + ndummy
        end do

!       update orbital offset
        ndummy = ndummy + dmrg_orbital_space%nash(isym)
      end do

!     last step: core energy
      write(fcidump,form3) corenergy , 0,0 ,0,0
!
      close(unit=fcidump,status='KEEP')

      end subroutine dmrg_task_fcidump
!**********************************************************************

      subroutine dmrg_task_tradump(tra,no,ni,jorb,tirrep)
!     -----------------------------------------------------------------
!
!     purpose: write one-electron integrals for a given orbital
!              to perform a sequential MPS rotation (to biorthonormal basis).
!
!     filename: FCIDUMP.tramps.orb#JORB
!     -----------------------------------------------------------------
!
      real*8 , intent(in), dimension(no,no)  :: tra
      integer, intent(in)                    :: jorb
      integer, intent(in)                    :: tirrep
      integer, intent(in)                    :: no
      integer, intent(in)                    :: ni
!     -----------------------------------------------------------------
      integer                                :: isym
      integer                                :: i, lfname, lsuffix
      integer                                :: ndummy
      integer                                :: norbtot
      character(len=30)                      :: form1
      character(len=30)                      :: form2
      character(len=100)                     :: filename
      character(len=4)                       :: suffix
      integer, parameter                     :: mpsdump   = 99
      real*8 , parameter                     :: threshold = 1.0d-50
      real*8                                 :: tjj
!     -----------------------------------------------------------------

      if(jorb > 999)then
        write(6,*) " Error: too many orbitals (>999) for MPS rotation!"
        stop
      end if

      ndummy = 0
      do isym = 1, tirrep - 1
        ndummy = ndummy + dmrg_orbital_space%nash(isym)
      end do
!
      suffix   = ""
      filename = ""
      lfname   = 0
      lsuffix  = 0

      if(jorb < 10)then
        write(suffix,'(i1)') jorb; lsuffix = 1
      else if(jorb < 100)then
        write(suffix,'(i2)') jorb; lsuffix = 2
      else if(jorb < 1000)then
        write(suffix,'(i3)') jorb; lsuffix = 3
      end if

      filename = "tjj.tramps.orb."//trim(suffix(1:lsuffix))
      lfname   = len_trim(filename)

      open(mpsdump,file=filename(1:lfname),status='replace',form='formatted',    &
           action='readwrite',position='rewind')

      !> jorb == 0: scaling factor for the rotation wrt the inactive orbital space
      if(jorb == 0)then
        tjj = tra(1,1)
      else
        tjj = tra(ni+jorb-ndummy,ni+jorb-ndummy)
      end if

      write(mpsdump,'(G21.12)') tjj
      close(unit=mpsdump,status='KEEP')

      if(jorb == 0) return

      filename = "FCIDUMP.tramps.orb."//trim(suffix(1:lsuffix))
      lfname   = len_trim(filename)

      open(mpsdump,file=filename(1:lfname),status='replace',form='formatted',    &
           action='readwrite',position='rewind')
!
!     define printing format | Using the same G21.12 as molpro
      form1="(G21.12, 4X, I6, I6, I6, I6)"
      form2="(A11,I3,A7,I2,A5,I2,A1)"

!     calculate total number of active orbitals
      norbtot = 0
      do isym = 1, dmrg_symmetry%nirrep
        norbtot = norbtot + dmrg_orbital_space%nash(isym)
      end do

!     !> Print the header of the FCIDUMP.tramps.orb.xxx file using the standard (FCIDUMP) MOLPRO format
      write(mpsdump,form2) ' &FCI NORB=', norbtot , ',NELEC=',          &
      dmrg_state%nactel, ',MS2=', dmrg_state%ms2, ','

      write(mpsdump,"(A)",advance='no') '  ORBSYM='
      do isym = 1, dmrg_symmetry%nirrep
        if(dmrg_orbital_space%nash(isym) /= 0)then
          do i = 1, dmrg_orbital_space%nash(isym)
            write(mpsdump,"(I1,A1)",advance='no') isym,','
          end do
        end if
      end do

      write(mpsdump,"(/A7,I1)") '  ISYM=',dmrg_state%irefsm
      write(mpsdump,"(A5)") ' &END'

!       !> punch out the transformation matrix AKA one-electron integrals
      do i = 1, dmrg_orbital_space%nash(tirrep)

        if(i == jorb-ndummy) cycle ! skip element tra_jj

        if(dabs(tra(ni+i,ni+jorb-ndummy)) < threshold)then
          cycle
        else
          write(mpsdump,form1) (tra(ni+i,ni+jorb-ndummy)/tjj), i+ndummy, jorb,0, 0
        end if
      end do

      close(unit=mpsdump,status='KEEP')

      end subroutine dmrg_task_tradump
!**********************************************************************

      subroutine dmrg_task_overlap(iroot,jroot,overlap,checkpoint1,checkpoint2)

      real*8 , optional,               intent(inout) :: overlap
      integer, optional,               intent(in)    :: iroot
      integer, optional,               intent(in)    :: jroot
      ! Leon 02-12-2016: added optional custom checkpoint names
      character(len=*),optional,       intent(in)    :: checkpoint1
      character(len=*),optional,       intent(in)    :: checkpoint2
!----------------------------------------------------------------------
      character(len=2300),              allocatable  :: maquis_name_states(:)
      character(len=3000)                            :: currdir
      character(len=300)                             :: pydriver
      character(len=300)                             :: overlap_exe
      integer                                        :: xroot
      integer                                        :: i
      integer                                        :: j
!----------------------------------------------------------------------

      allocate(maquis_name_states(2))
      maquis_name_states     = ""
      if (present(checkpoint1).and.present(checkpoint2)) then
        !> set prefix for QCMaquis results/checkpoint files
        call getenv("CurrDir",currdir)
        maquis_name_states(1) = trim(currdir)//'/'//trim(checkpoint1)
        maquis_name_states(2) = trim(currdir)//'/'//trim(checkpoint2)
      else if (present(iroot).and.present(jroot)) then
        do j = 1, 2
                  xroot=iroot-1
          if(j == 2) xroot=jroot-1

          if(xroot.lt.1000)then
            call file_name_generator(xroot,"checkpoint_state.",17,".h5",3, maquis_name_states(j))
          else
            write(6,*)"too many states (>999) in DMRG calculation"
            Stop
          end if
        end do
      else
        write (6,*) 'Mandatory parameters (root # or checkpoint name) for dmrg_task_overlap missing'
        stop 99
      end if
      overlap = 0.0d0

      pydriver    = "$MOLCAS/pytools/runDMRG.py "
      overlap_exe = "mps_overlap_su2u1pg "

      !> compute overlap <MPS1 | MPS2>
      call overlap_kernel(pydriver,overlap_exe,maquis_name_states,overlap)

      deallocate(maquis_name_states)

      end subroutine dmrg_task_overlap
!**********************************************************************

      subroutine dmrg_task_overlapR(jroot,overlap)

      real*8 , optional,               intent(inout) :: overlap
      integer, optional,               intent(in)    :: jroot
!----------------------------------------------------------------------
      character(len=2300),              allocatable  :: maquis_name_states(:)
      character(len=300)                             :: pydriver
      character(len=300)                             :: overlap_exe
      integer                                        :: xroot
      integer                                        :: j
!----------------------------------------------------------------------

      allocate(maquis_name_states(2))
      maquis_name_states     = ""

      write(maquis_name_states(1),'(a)') "rf.checkpoint_state.h5"
      if(jroot > 0)then
        do j = 2, 2
          xroot=jroot-1
          if(xroot.lt.1000)then
            call file_name_generator(xroot,"checkpoint_state.",17,".h5",3, maquis_name_states(j))
          else
            write(6,*)"too many states (>999) in DMRG calculation"
            Stop
          end if
        end do
      else
        write(maquis_name_states(2),'(a)') "rf.checkpoint_state.h5"
      end if

      overlap = 0.0d0

      pydriver    = "$MOLCAS/pytools/runDMRG.py "
      overlap_exe = "mps_overlap_su2u1pg "

      !> compute overlap <MPS1 | MPS2>
      call overlap_kernel(pydriver,overlap_exe,maquis_name_states,overlap)

      deallocate(maquis_name_states)

      end subroutine dmrg_task_overlapR

!**********************************************************************
! Leon: 2U1PG overlap. For now, only with checkpoint names
! TODO: maybe allow it to work with iroot/jroot/multiplet
!
      subroutine dmrg_task_overlapU(overlap,checkpoint1,checkpoint2)

      real*8,                 intent(inout) :: overlap
      character(len=*),       intent(in)    :: checkpoint1
      character(len=*),       intent(in)    :: checkpoint2
!----------------------------------------------------------------------
      character(len=2300)                   :: maquis_name_states(2)
      character(len=3000)                   :: currdir
      character(len=300)                    :: pydriver
      character(len=300)                    :: overlap_exe

      overlap = 0.0d0

      pydriver    = "$MOLCAS/pytools/runDMRG.py "
      overlap_exe = "mps_overlap_2u1pg "

      !> set prefix for QCMaquis results/checkpoint files
      call getenv("CurrDir",currdir)
      maquis_name_states(1) = trim(currdir)//'/'//trim(checkpoint1)
      maquis_name_states(2) = trim(currdir)//'/'//trim(checkpoint2)

      print *, 'Calculating overlap of '//trim(checkpoint1)//' and '//trim(checkpoint2)// ' in 2U1 group'

      call overlap_kernel(pydriver,overlap_exe,maquis_name_states,overlap)

      end subroutine dmrg_task_overlapU
!**********************************************************************

      subroutine dmrg_task_MPScrot(iroot,checkpoint)

      integer, optional, intent(in)   :: iroot
      ! Leon 02-12-2016: added optional custom checkpoint names
      character(len=*),optional,intent(in) :: checkpoint
!----------------------------------------------------------------------
      character(len=300)              :: pydriver
      character(len=300)              :: MPScrot_exe
      character(len=300)              :: project
      character(len=600)              :: rotout
      character(len=256)              :: state_tag
      integer                         :: xroot, i, j, jj
      character(len=2300)             :: maquis_name_state
      character(len=4)                :: suffix
      character                       :: cJ,cJJ,cJJJ
      character(len=100)              :: maquis_msproj
      character(len=3000)             :: currdir
!----------------------------------------------------------------------

      call getenv("CurrDir",currdir)
      state_tag(1:5) = ' '
      pydriver       = "$MOLCAS/pytools/runDMRG.py "
      MPScrot_exe    = "mps_rotate_2u1pg "

      maquis_name_state     = ""

      if (present(iroot)) then
        xroot=iroot-1

        !> select target MPS to rotate
        if(xroot < 1000)then
          call file_name_generator(xroot,"checkpoint_state.",17,".h5",3, maquis_name_state)
        else
          write(6,*)"too many states (>999) in DMRG calculation"
          stop
        end if
      else if (present(checkpoint)) then
        maquis_name_state = trim(currdir)//'/'//checkpoint
      else
        write(6,*) "essential parameter (root # or checkpoint name) to dmrg_task_MPScrot missing"
        stop 99
      end if

      !> select Ms value for 2u1(pg) state
      if(abs(dmrg_state%ms2) < 10)then
          J=abs(dmrg_state%ms2)
          cJ=CHAR(J+48)
          if(dmrg_state%ms2 < 0)then
            suffix = "-"//cJ
          else
            suffix = cJ
          end if
      else if(abs(dmrg_state%ms2) < 100)then
          J=mod(abs(dmrg_state%ms2),10)
          JJ=abs(dmrg_state%ms2)/10
          cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          if(dmrg_state%ms2 < 0)then
            suffix = "-"//cJJ//cJ
          else
            suffix = cJJ//cJ
          end if
      end if
      write(maquis_msproj,'(a,a)')   ' --msprojRHS=',trim(suffix)

      print *, 'msRHS == ',suffix

      !> compute counterrotated MPS
      call system(trim(pydriver)//' --rotate '//                          &
                  " --executable="//trim(MPScrot_exe)//" "//              &
                  " --output=mps.rotate.out"//' --tmpfull=$PWD/tmp'//" "//&
                  trim(maquis_msproj)//" "//                              &
                  ' --rhs='//trim(maquis_name_state)                      &
                 )

      !> save output for checking/debugging
      call getenv("Project",Project)
      if (present(iroot)) then
        call get_state_tag(iroot,state_tag,dmrg_file%offset)
      else if(present(checkpoint))then
          state_tag = checkpoint
        ! exception when neither iroot nor checkpoint are present is caught above
      end if
      rotout  = Project(1:index(Project,' ')-1)//'.QCMaquis.state.'//trim(state_tag)//'.mps.rotate.out'
      call system("cp mps.rotate.out "//trim(rotout))

      end subroutine dmrg_task_MPScrot
!**********************************************************************

      subroutine dmrg_task_MPSbackup(iroot,switch,checkpoint)

      integer, optional, intent(in)    :: iroot
      integer, optional, intent(in)    :: switch
      ! Leon 02-12-2016: added optional custom checkpoint names
      character(len=*),optional,intent(in) :: checkpoint
!----------------------------------------------------------------------
      character(len=2300)              :: maquis_name_state
      integer                          :: xroot
      logical                          :: available
      character(len=3000)              :: currdir
!----------------------------------------------------------------------

      call getenv("CurrDir",currdir)

      maquis_name_state     = ""

      if (present(iroot)) then
        xroot=iroot-1
        if(xroot < 1000)then
          call file_name_generator(xroot,"checkpoint_state.",17,".h5",3, maquis_name_state)
        else
          write(6,*)"too many states (>999) in DMRG calculation"
          stop
        end if
      else if (present(checkpoint)) then
        maquis_name_state = trim(currdir)//'/'//checkpoint
      else
        write(6,*) "essential parameter (root # or checkpoint name) to dmrg_task_MPSbackup missing"
        stop 99
      end if

      !> make backup of MPS or restore MPS
      select case(switch)
      case(-1)
        inquire(file=trim(maquis_name_state)//"/props.h5",exist=available)
        if(.not.available) print *, 'cannot back up non-existing MPS: ',trim(maquis_name_state)
        call system("cp -r "//trim(maquis_name_state)//" "//trim(maquis_name_state)//".BACKUP")
      case( 1)
        inquire(file=trim(maquis_name_state)//".BACKUP"//"/props.h5",exist=available)
        if(.not.available) print *, 'cannot restore non-existing MPS: ',trim(maquis_name_state)
        call system("cp -r "//trim(maquis_name_state)//".BACKUP/* "//trim(maquis_name_state))
      case default
         print *, 'MPS backup - unknown case'
      end select

      end subroutine dmrg_task_MPSbackup
!**********************************************************************

      subroutine dmrg_task_import_spd(spd1,nrdm,iroot,sp1)

      real*8 , optional, dimension(*), intent(inout) :: spd1
      integer, optional,               intent(in)    :: nrdm
      integer, optional,               intent(in)    :: iroot
      logical, optional,               intent(in)    :: sp1

      if(sp1) call dmrg_task_import_spd1(spd1=spd1,nrdm1=nrdm,iroot=iroot)

      end subroutine dmrg_task_import_spd

      subroutine dmrg_task_import_spd1(spd1,nrdm1,iroot)

!     -----------------------------------------------------------------
!     !> purpose: import the 1-particle spin-density matrix generated by DMRG
!     -----------------------------------------------------------------
      integer, intent(in)    :: nrdm1
      integer, intent(in)    :: iroot
      real*8,  intent(inout) :: spd1(nrdm1)
!     -----------------------------------------------------------------
      real*8, allocatable    :: tmp(:,:)
      character(len=2300)    :: oneSPDfile
      integer                :: i,j,ij,i1,j1,nact, irootm1
      integer                :: lunit
!     -----------------------------------------------------------------

        lunit      =   140
! Since DMRG start from state - 0
        irootm1    = iroot - 1
        oneSPDfile = ""

        if(irootm1.lt.1000)then
          call file_name_generator(irootm1,"oneparticle.spd.",16,"",0,oneSPDfile)
        else
          write(6,*)"There are too many states (>999) in DMRG calculation"
          Stop
        end if

#ifdef _DMRG_DEBUG_
      write(6,*)oneSPDfile
#endif

      open(unit=lunit,file=trim(oneSPDfile))
        read(lunit,*)nact

        allocate(tmp(nact,nact)); tmp = 0

        do i=1,nact
          do j=1,nact
            read(lunit,*)i1,j1,tmp(i1+1,j1+1)
          end do
        end do
      close(lunit)

      ij=0
      do i=1,nact
        do j=1,i
          ij=ij+1
          spd1(ij)=tmp(i,j)
        end do
      end do

#ifdef _DMRG_DEBUG_
      write(6,*)"SPD1 successfully imported"
#endif

      deallocate(tmp)

      end subroutine dmrg_task_import_spd1
!**********************************************************************

      subroutine dmrg_task_import_rdmY(dv,pv,tv,nrdm1,nrdm2,nrdm3,iroot,jroot,&
                                       msproj,msprojL,multiplet,multipletL,rdm1,rdm2,rdm3,checkpoint1,checkpoint2)

      real*8 , optional, dimension(*), intent(inout) :: dv
      real*8 , optional, dimension(*), intent(inout) :: pv
      real*8 , optional, dimension(*), intent(inout) :: tv
      integer, optional,               intent(in)    :: nrdm1
      integer, optional,               intent(in)    :: nrdm2
      integer, optional,               intent(in)    :: nrdm3
      integer, optional,               intent(in)    :: iroot
      integer, optional,               intent(in)    :: jroot
      integer, optional,               intent(in)    :: msproj
      integer, optional,               intent(in)    :: msprojL
      integer, optional,               intent(in)    :: multiplet
      integer, optional,               intent(in)    :: multipletL
      logical, optional,               intent(in)    :: rdm1
      logical, optional,               intent(in)    :: rdm2
      logical, optional,               intent(in)    :: rdm3
      ! Leon 02-12-2016: added optional custom checkpoint names
      character(len=*),optional,       intent(in)    :: checkpoint1
      character(len=*),optional,       intent(in)    :: checkpoint2

!-------------------------------------------------------------------------------
      integer                                        :: setup_container(6)
!-------------------------------------------------------------------------------

      !> we proceed in two steps:
        !> step 1: compute the TDM(s): <iroot | c+(c+[c+]) c(c[c]) | jroot>
        !> step 2: import the TDM(s)

      !> temporary fix until I understand the problem of why TDMs for <S|op|T> and <T|op|S> are not identical
      !  where T == triplet state and S = singlet state
      !  always calculate <T|op|S> for now...
      !setup_container(1:6) = -100000
      if (present(iroot).and.present(jroot))then
        if(present(multipletL) .and. present(multiplet)) then
          if (multipletL > multiplet)then
            setup_container(1:6) = (/ iroot,jroot,msproj,msprojL,multiplet,multipletL /)
          else
            setup_container(1:6) = (/ jroot,iroot,msprojL,msproj,multipletL,multiplet /)
          end if
        end if

        if(present(rdm1) .and. rdm1)then
          call dmrg_task_compute_rdmY(iroot     = setup_container(1),   jroot    = setup_container(2), &
                                      msproj    = setup_container(3), msprojL    = setup_container(4), &
                                      multiplet = setup_container(5), multipletL = setup_container(6), isrdm=1)
          call dmrg_task_import_rdmY1(dv=dv,nrdm1=nrdm1,iroot=setup_container(1),jroot=setup_container(2))
        else if(present(rdm2) .and. rdm2)then
          call dmrg_task_compute_rdmY(iroot     = setup_container(1),   jroot    = setup_container(2), &
                                      msproj    = setup_container(3), msprojL    = setup_container(4), &
                                      multiplet = setup_container(5), multipletL = setup_container(6), isrdm=2)
          call dmrg_task_import_rdmY2(pv=pv,nrdm2=nrdm2,iroot=setup_container(1),jroot=setup_container(2))
        else if(present(rdm3) .and. rdm3)then
          call dmrg_task_compute_rdmY(iroot=iroot, jroot=jroot, isrdm=3)
          call dmrg_task_import_rdmY3(tv=tv,nrdm3=nrdm3,iroot=iroot,jroot=jroot)
        else
          print *, ' this routine currently supports only 1p-/2p-/3p-transition density matrices'
          stop 99
        end if
      else if (present(checkpoint1).and.present(checkpoint2))then
        if (multipletL > multiplet)then
          setup_container(3:6) = (/ msproj,msprojL,multiplet,multipletL /)
          !write(6,*) 'container... ',setup_container(3:6)
          if(present(rdm1) .and. rdm1)then
            call dmrg_task_compute_rdmY(msproj    = setup_container(3), msprojL    = setup_container(4), &
                                        multiplet = setup_container(5), multipletL = setup_container(6), isrdm=1, &
                                        checkpoint1 = checkpoint2, checkpoint2 = checkpoint1)
            call dmrg_task_import_rdmY1(dv=dv,nrdm1=nrdm1,checkpoint1 = checkpoint2, checkpoint2 = checkpoint1)
          else if(present(rdm2) .and. rdm2)then
            call dmrg_task_compute_rdmY(msproj    = setup_container(3), msprojL    = setup_container(4), &
                                        multiplet = setup_container(5), multipletL = setup_container(6), isrdm=2, &
                                        checkpoint1 = checkpoint2, checkpoint2 = checkpoint1)
            call dmrg_task_import_rdmY2(pv=pv,nrdm2=nrdm2,checkpoint1 = checkpoint2, checkpoint2 = checkpoint1)
          else
            print *,'Only up to 2p-TDMs are supported with checkpoint names for now'
            stop 99
          end if
        else
          setup_container(3:6) = (/ msprojL,msproj,multipletL,multiplet /)
          !write(6,*) 'container... ',setup_container(3:6)
          if(present(rdm1) .and. rdm1)then
            call dmrg_task_compute_rdmY(msproj    = setup_container(3), msprojL    = setup_container(4), &
                                        multiplet = setup_container(5), multipletL = setup_container(6), isrdm=1, &
                                        checkpoint1 = checkpoint1, checkpoint2 = checkpoint2)
            call dmrg_task_import_rdmY1(dv=dv,nrdm1=nrdm1,checkpoint1 = checkpoint1, checkpoint2 = checkpoint2)
          else if(present(rdm2) .and. rdm2)then
            call dmrg_task_compute_rdmY(msproj    = setup_container(3), msprojL    = setup_container(4), &
                                        multiplet = setup_container(5), multipletL = setup_container(6), isrdm=2, &
                                        checkpoint1 = checkpoint1, checkpoint2 = checkpoint2)
            call dmrg_task_import_rdmY2(pv=pv,nrdm2=nrdm2,checkpoint1 = checkpoint1, checkpoint2 = checkpoint2)
          else
            print *,'Only up to 2p-TDMs are supported with checkpoint names for now'
            stop 99
          end if
        end if
      else
        print *,'Mandatory parameters (root # or checkpoint name) for dmrg_task_import_rdmY missing'
        stop 99
      end if
      end subroutine dmrg_task_import_rdmY
!**********************************************************************

      subroutine dmrg_task_compute_rdmY(iroot,jroot,msproj,msprojL,multiplet,multipletL,isrdm,checkpoint1,checkpoint2)

      integer, optional,               intent(in)    :: iroot
      integer, optional,               intent(in)    :: jroot
      integer, optional,               intent(in)    :: msproj
      integer, optional,               intent(in)    :: msprojL
      integer, optional,               intent(in)    :: multiplet
      integer, optional,               intent(in)    :: multipletL
      ! Leon 02-12-2016: added optional custom checkpoint names
      integer,                         intent(in)    :: isrdm
      character(len=*),optional,       intent(in)    :: checkpoint1
      character(len=*),optional,       intent(in)    :: checkpoint2
!----------------------------------------------------------------------
      character(len=300)                             :: input
      character(len=2300),             allocatable   :: maquis_name_states(:)
      character(len=2300),             allocatable   :: maquis_name_results(:)
      integer                                        :: xroot
      integer                                        :: i
      integer                                        :: j,jj
      character(len=256),              allocatable   :: state_tag(:)
      logical                                        :: available
      character(len=3000)    :: currdir
      character(len=4)       :: suffix(4)
      character              :: cJ,cJJ,cJJJ
      !> for the input file
      character(len=100)     :: maquis_norb
      character(len=100)     :: maquis_msproj
      character(len=100)     :: maquis_msprojL
      character(len=2500)    :: maquis_result
      character(len=8000)    :: option
      character(len=100)     :: maquis_nup
      character(len=100)     :: maquis_ndown

!----------------------------------------------------------------------

      call getenv("CurrDir",currdir)

      allocate(maquis_name_states(2), maquis_name_results(2), state_tag(2))

      maquis_name_states  = ""
      maquis_name_results = ""
      state_tag           = ""
      suffix              = ""

      if(present(msprojL))then
        if(abs(msprojL) < 10)then
          J=abs(msprojL)
          cJ=CHAR(J+48)
          if(msprojL < 0)then
            suffix(1) = "-"//cJ
          else
            suffix(1) = cJ
          end if
        else if(abs(msprojL) < 100)then
          J=mod(abs(msprojL),10)
          JJ=abs(msprojL)/10
          cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          if(msprojL < 0)then
            suffix(1) = "-"//cJJ//cJ
          else
            suffix(1) = cJJ//cJ
          end if
        end if
        write(maquis_msprojL,'(a,a)')   ' --msprojLHS=',trim(suffix(1))
      else
        maquis_msprojL = ' '
      end if

      if(present(msproj))then
        if(abs(msproj) < 10)then
          J=abs(msproj)
          cJ=CHAR(J+48)
          if(msproj < 0)then
            suffix(2) = "-"//cJ
          else
            suffix(2) = cJ
          end if
        else if(abs(msproj) < 100)then
          J=mod(abs(msproj),10)
          JJ=abs(msproj)/10
          cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          if(msproj < 0)then
            suffix(2) = "-"//cJJ//cJ
          else
            suffix(2) = cJJ//cJ
          end if
        end if
        write(maquis_msproj,'(a,a)')   ' --msprojRHS=',trim(suffix(2))
      else
        maquis_msproj = ' '
      end if

      if(present(multipletL))then
        if(abs(multipletL) < 10)then
          J=abs(multipletL)
          cJ=CHAR(J+48)
          if(multipletL < 0)then
            suffix(3) = "-"//cJ
          else
            suffix(3) = cJ
          end if
        else if(abs(multipletL) < 100)then
          J=mod(abs(multipletL),10)
          JJ=abs(multipletL)/10
          cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          if(multipletL < 0)then
            suffix(3) = "-"//cJJ//cJ
          else
            suffix(3) = cJJ//cJ
          end if
        end if
      end if

      if(present(multiplet))then
        if(abs(multiplet) < 10)then
          J=abs(multiplet)
          cJ=CHAR(J+48)
          if(multiplet < 0)then
            suffix(4) = "-"//cJ
          else
            suffix(4) = cJ
          end if
        else if(abs(multiplet) < 100)then
          J=mod(abs(multiplet),10)
          JJ=abs(multiplet)/10
          cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          if(multiplet < 0)then
            suffix(4) = "-"//cJJ//cJ
          else
            suffix(4) = cJJ//cJ
          end if
        end if
      end if

      if (present(iroot).and.present(jroot)) then
        do j = 1, 2

                    xroot=jroot-1
          if(j == 2) xroot=iroot-1

          if(xroot.lt.1000)then
            !> MPS rotated - we already have the MPS in 2u1(pg) representation
            if(dmrg_external%MPSrotated)then
              call file_name_generator(xroot,"checkpoint_state.",17,"."//trim(suffix(3+j-1))//"."//trim(suffix(1+j-1))//&
                                      ".h5",len(trim(suffix(3+j-1)))+len(trim(suffix(1+j-1)))+5, maquis_name_states(j))
            else
              call file_name_generator(xroot,"checkpoint_state.",17,".h5",3, maquis_name_states(j))
            end if
            call file_name_generator(xroot,"results_state.",   14,".h5",3, maquis_name_results(j))
            call get_state_tag(xroot+1,state_tag(j),dmrg_file%offset)
          else
            write(6,*)"too many states (>999) in DMRG calculation"
            Stop
          end if
        end do
      else if (present(checkpoint1).and.present(checkpoint2)) then
        if(dmrg_external%MPSrotated)then
          ! in this case we assume to have 2U1 checkpoint names with corresponding naming convention
          ! project_name.checkpoint.state.<state#>.<s>.<ms>.h5
          ! but checkpoint names rasscf provides us are SU2U1 names, so we're going to generate 2U1 names from them
          ! this works only if the SU2U1 checkpoint is named properly (ie as project_name.checkpoint_state.<state#>.h5)

          ! we're cutting off the '.h5' extension and adding the '.<s>.<ms>.h5' suffix

          maquis_name_states(1) = checkpoint1(1:len_trim(checkpoint1)-3)//"."//trim(suffix(3))//"."//    &
                        trim(suffix(1))//".h5"
          maquis_name_states(2) = checkpoint2(1:len_trim(checkpoint2)-3)//"."//trim(suffix(4))//"."//    &
                        trim(suffix(2))//".h5"
        else
          ! otherwise assume we have SU2U1 checkpoints with proper names
          maquis_name_states(1) = trim(checkpoint1)
          maquis_name_states(2) = trim(checkpoint2)
        end if
        !> leon+stefan: use path to currdir to correctly address the checkpoint file
        maquis_name_states(1) = trim(currdir)//'/'//trim(maquis_name_states(1))
        maquis_name_states(2) = trim(currdir)//'/'//trim(maquis_name_states(2))

        ! The result file naming is a bit ugly in this case, since instead of root numbers in the file name we should have full checkpoint names
        maquis_name_results(1) = trim(checkpoint1)//".tdm_results.h5"
        maquis_name_results(2) = trim(checkpoint2)//".tdm_results.h5"
        ! state tags are just full checkpoint names
        state_tag(1) = trim(checkpoint1)
        state_tag(2) = trim(checkpoint2)
      end if
      !> prepare replacements in template input file
      write(maquis_norb,'(a,i3,a)')    ' --replace="orbital_number=',dmrg_external%norb,'"'
      write(maquis_result,'(a,a,a)')   ' --replace="saved_result='//trim(maquis_name_results(2))//'"'

      if(dmrg_external%MPSrotated)then

        write(maquis_nup,'(a,i3,a)')   ' --replace="isup= ', dmrg_external%nalpha,'"'
        write(maquis_ndown,'(a,i3,a)') ' --replace="isdown= ', dmrg_external%nbeta,'"'

        write(option,'(16a)')                                                                     &
           ' --notransform --replace="bra_checkpoint='//trim(maquis_name_states(1))//'"'//' '//   &
           ' --replace="saved_checkpoint='//trim(maquis_name_states(2))//'"'//' '//               &
           ' --replace="isirrep='//"0"//'"'//' '//trim(maquis_nup)//' '//trim(maquis_ndown)//' '
      else
        option = ' '
      end if

      !print *, 'OPTION is for TDMs == ',option

      input = ' $MOLCAS/template-files/template-dmrg-tdm.maquis'

      select case(isrdm)
      case(1)

        inquire(file="oneparticle.tdm."//trim(state_tag(1))//"."//trim(state_tag(2)),exist=available)

        if(available)then
          write(6,*)'import existing 1-TDM for <',trim(state_tag(1)),'|c+ c|',trim(state_tag(2)),'>'
        else
          write(6,*) 'compute 1-TDM for <',trim(state_tag(1)),'|c+ c|',trim(state_tag(2)),'>'

! write(6,*) 'these are my commands: '

! write(6,*) '$MOLCAS/pytools/runDMRG.py '//' --tdm '//' --tdmlevel=1 '//                                   &
!                      " --executable=dmrg_meas"//                                                                    &
!                      " --output=dmrg.TDM_OUT "//' --tmpfull=$PWD/tmp'//                                             &
!                      ' --lhs='//trim(maquis_name_states(1))//' --rhs='//trim(maquis_name_states(2))//               &
!                      trim(maquis_msprojL)//' '//trim(maquis_msproj)//' '//                                          &
!                      trim(maquis_norb)//trim(maquis_result)//' '//trim(option)//' '//trim(input)

          !> off-diagonal part of 1-TDM
          call system('$MOLCAS/pytools/runDMRG.py '//' --tdm '//' --tdmlevel=1 '//                                   &
                      " --executable=dmrg_meas"//                                                                    &
                      " --output=dmrg.TDM_OUT "//' --tmpfull=$PWD/tmp'//                                             &
                      ' --lhs='//trim(maquis_name_states(1))//' --rhs='//trim(maquis_name_states(2))//               &
                      trim(maquis_msprojL)//' '//trim(maquis_msproj)//' '//                                          &
                      trim(maquis_norb)//trim(maquis_result)//' '//trim(option)//' '//trim(input)                    &
                     )

!write(6,*) '$MOLCAS/pytools/runDMRG.py '//' --tdm '//' --tdmlevel=1 '//                                   &
!                      " --executable=onetdm_diag_2u1pg"//                                                            &
!                      " --output=dmrg.diagTDM_OUT "//' --tmpfull=$PWD/tmp'//                                         &
!                      ' --lhs='//trim(maquis_name_states(1))//' --rhs='//trim(maquis_name_states(2))//               &
!                      trim(maquis_msprojL)//' '//trim(maquis_msproj)//' '//                                          &
!                      trim(maquis_norb)//trim(maquis_result)//' '//trim(option)//' '//trim(input)

          !> diagonal part of 1-TDM
          call system('$MOLCAS/pytools/runDMRG.py '//' --tdm '//' --tdmlevel=1 '//                                   &
                      " --executable=onetdm_diag_2u1pg"//                                                            &
                      " --output=dmrg.diagTDM_OUT "//' --tmpfull=$PWD/tmp'//                                         &
                      ' --lhs='//trim(maquis_name_states(1))//' --rhs='//trim(maquis_name_states(2))//               &
                      trim(maquis_msprojL)//' '//trim(maquis_msproj)//' '//                                          &
                      trim(maquis_norb)//trim(maquis_result)//' '//trim(option)//' '//trim(input)                    &
                      )

!write(6,*) "$MOLCAS/pytools/tdmsave.py "//trim(maquis_name_results(2))//" "//trim(state_tag(1))//" "&
!                      //trim(state_tag(2))//" 1 > tdm_save."//trim(state_tag(1))//"."                          &
!                      //trim(state_tag(2))//".out"

          !> save result in oneparticle.tdm.(iroot-1).(jroot-1)
          call system("$MOLCAS/pytools/tdmsave.py "//trim(maquis_name_results(2))//" "//trim(state_tag(1))//" "&
                      //trim(state_tag(2))//" 1 > tdm_save."//trim(state_tag(1))//"."                          &
                      //trim(state_tag(2))//".out"                                                             &
                     )
          print *, 'saving 1-TDM to oneparticle.tdm.'//trim(state_tag(1))//"."//trim(state_tag(2))
        end if

      case(2)

        inquire(file="twoparticle.tdm."//trim(state_tag(1))//"."//trim(state_tag(2)),exist=available)

        if(available)then
          write(6,*)'import existing 2-TDM for <',trim(state_tag(1)),'|c+c+ cc|',trim(state_tag(2)),'>'
        else
          write(6,*) 'compute 2-TDM for <',trim(state_tag(1)),'|c+c+ cc|',trim(state_tag(2)),'>'

          !> 2-TDM
          call system('$MOLCAS/pytools/runDMRG.py '//' --tdm '//' --tdmlevel=2 '//                             &
                      " --executable=dmrg_meas"//                                                              &
                      " --output=dmrg.TDM_OUT "//' --tmpfull=$PWD/tmp'//                                       &
                      ' --lhs='//trim(maquis_name_states(1))//' --rhs='//trim(maquis_name_states(2))//         &
                      trim(maquis_msprojL)//' '//trim(maquis_msproj)//' '//                                    &
                      trim(maquis_norb)//trim(maquis_result)//' '//trim(option)//' '//trim(input)              &
                     )

          !> save result in twoparticle.tdm.(iroot-1).(jroot-1)
          call system("$MOLCAS/pytools/tdmsave.py "//trim(maquis_name_results(2))//" "//trim(state_tag(1))//" "&
                      //trim(state_tag(2))//" 2 > tdm_save."//trim(state_tag(1))//"."                          &
                      //trim(state_tag(2))//".out"                                                             &
                     )

        end if
      case(3)

        !> check both options x.y and y.x as state tags for the 3p-tdms
        inquire(file="threeparticle.tdm."//trim(state_tag(1))//"."//trim(state_tag(2)),exist=available)
        if(.not.available)then
          inquire(file="threeparticle.tdm."//trim(state_tag(2))//"."//trim(state_tag(1)),exist=available)
        end if

        if(available)then
          write(6,*)'import existing 3-TDM for <',trim(state_tag(1)),'|c+c+c+ ccc|',trim(state_tag(2)),'>'
        else
          write(6,*)'compute 3-TDM for <',trim(state_tag(1)),'|c+c+c+ ccc|',trim(state_tag(2)),'>'

          ! 3-TDM
          call system('$MOLCAS/pytools/runDMRG.py '//' --tdm '//' --tdmlevel=3 '//                             &
                      " --executable=dmrg_meas"//" --output=dmrg.TDM_OUT "//' --tmpfull=$PWD/tmp'//            &
                      ' --lhs='//trim(maquis_name_states(1))//' --rhs='//trim(maquis_name_states(2))//         &
                      trim(maquis_msprojL)//' '//trim(maquis_msproj)//' '//                                    &
                      trim(maquis_norb)//trim(maquis_result)//trim(input)                                      &
                     )

          !> save result in threeparticle.tdm.(iroot-1).(jroot-1)
          call system("$MOLCAS/pytools/transition_threeptrdm.py "//trim(maquis_name_results(2))//" "//         &
                      trim(state_tag(1))//" "  &
                      //trim(state_tag(2))//" 3 > tdm_save."//trim(state_tag(1))//"."                          &
                      //trim(state_tag(2))//".out"                                                             &
                     )
        end if
      case default
        write(6,*) 'dmrg_task_compute_rdmY: nothing to do ...'
      end select


      deallocate(maquis_name_states, maquis_name_results, state_tag)

      end subroutine dmrg_task_compute_rdmY
!**********************************************************************

      subroutine dmrg_task_import_rdmY1(dv,nrdm1,iroot,jroot,checkpoint1,checkpoint2)

!     -----------------------------------------------------------------
!     !> purpose: import the special version of the 1p-reduced density
!                 matrix generated by DMRG
!     !> required for RASSI interface in Molcas
!
!     !> TDMs are NOT triangular packed!
!
!     -----------------------------------------------------------------
      real*8,  intent(inout) :: dv(nrdm1**2)
      integer, intent(in)    :: nrdm1
      integer, optional, intent(in)    :: iroot
      integer, optional, intent(in)    :: jroot
      character(len=*),optional,       intent(in)    :: checkpoint1
      character(len=*),optional,       intent(in)    :: checkpoint2
!     -----------------------------------------------------------------
      real*8, allocatable    :: tmp(:,:)
      integer                :: i,ii,iii,j,jj,jjj,ij,i1,j1,nact,k,isorb,jsorb
      integer                :: irootm1, jrootm1, len_suffix
      character(len=4)       :: suffix
      character(len=2300)    :: oneRDMfile
      character              :: cI,cII,CIII,cJ,cJJ,cJJJ
      integer                :: lunit
      character(len=1997)    :: dmrg_file_prefix_save
!     -----------------------------------------------------------------

      lunit   =   140
      oneRDMfile=""

      if (present(iroot).and.present(jroot)) then
        irootm1 = iroot - 1
        jrootm1 = jroot - 1

        if(irootm1 < 10)then
          J=irootm1
          cJ=CHAR(J+48)
          len_suffix = 2
          suffix(1:2) = "."//cJ
        else if(irootm1 < 100)then
          J=mod(irootm1,10)
          JJ=irootm1/10
          cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          len_suffix = 3
          suffix(1:3) = "."//cJJ//cJ
        else if(irootm1 < 1000)then
            J=mod(irootm1,10)
            JJ=mod(irootm1,100)/10
          JJJ=irootm1/100
            cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          cJJJ=CHAR(JJJ+48)
          len_suffix = 4
          suffix(1:4) = "."//cJJJ//cJJ//cJ
        end if

        if(irootm1 < 1000 .and. jrootm1 < 1000)then
          dmrg_file_prefix_save = dmrg_file%prefix
          dmrg_file%prefix = ''
          call file_name_generator(jrootm1,"oneparticle.tdm.",16,suffix,len_suffix,oneRDMfile)
          dmrg_file%prefix = dmrg_file_prefix_save
        else
          print *, "There are too many states (>999) in DMRG calculation"
          stop
        end if
      else if (present(checkpoint1).and.present(checkpoint2)) then
        oneRDMfile="oneparticle.tdm."//trim(checkpoint1)//"."//trim(checkpoint2)
      end if

!#ifdef _DMRG_DEBUG_
      print *, 'reading 1-TDM from file: ', trim(oneRDMfile)
!#endif

      open(unit=lunit,file=oneRDMfile)
        read(lunit,*)nact
        allocate(tmp(4*nact,4*nact))
!       !> import alpha-alpha (k == 1)
!                  alpha-beta (k == 2)
!                  beta-alpha (k == 3)
!                   beta-beta (k == 4)
!          blocks separately
         tmp = 0.0d0
        !> read unsorted TDM
        do k = 1, 4
          do i=1,nact
            do j=1,nact
              read(lunit,*)i1,j1,tmp(((i1+1)+(k-1)*nact),((j1+1)+(k-1)*nact))
            end do
          end do
        end do
      close(lunit)

      !> store TDM
#ifdef _DMRG_DEBUG_
      write(6,*) ' store sorted 1-TDM for lhs|rhs',jroot,iroot
#endif
      do i=1,nact
        isorb = 2*i - 1
        do j=1,nact
          jsorb = 2*j - 1
          dv(1-1+isorb+nrdm1*(jsorb-1)) = tmp(i+(1-1)*nact,j+(1-1)*nact) ! aa
          dv(1-1+isorb+nrdm1*(jsorb  )) = tmp(i+(2-1)*nact,j+(2-1)*nact) ! ab
          dv(1  +isorb+nrdm1*(jsorb-1)) = tmp(i+(3-1)*nact,j+(3-1)*nact) ! ba
          dv(1  +isorb+nrdm1*(jsorb  )) = tmp(i+(4-1)*nact,j+(4-1)*nact) ! bb
#ifdef _DMRG_DEBUG_
          if(abs(dv(1-1+isorb+nrdm1*(jsorb-1))) > 1.0d-14)&
          write(6,*) 'aa: index, value == ',1-1+isorb+nrdm1*(jsorb-1),dv(1-1+isorb+nrdm1*(jsorb-1))
          if(abs(dv(1-1+isorb+nrdm1*(jsorb  ))) > 1.0d-14)&
          write(6,*) 'ab: index, value == ',1-1+isorb+nrdm1*(jsorb  ),dv(1-1+isorb+nrdm1*(jsorb  ))
          if(abs(dv(1  +isorb+nrdm1*(jsorb-1))) > 1.0d-14)&
          write(6,*) 'ba: index, value == ',1  +isorb+nrdm1*(jsorb-1),dv(1  +isorb+nrdm1*(jsorb-1))
          if(abs(dv(1  +isorb+nrdm1*(jsorb  ))) > 1.0d-14)&
          write(6,*) 'bb: index, value == ',1  +isorb+nrdm1*(jsorb  ),dv(1  +isorb+nrdm1*(jsorb  ))
#endif
        end do
      end do

      deallocate(tmp)

      end subroutine dmrg_task_import_rdmY1
!**********************************************************************

      subroutine dmrg_task_import_rdmY2(pv,nrdm2,iroot,jroot,checkpoint1,checkpoint2)

!     -----------------------------------------------------------------
!     purpose: import the 2p-reduced transition density matrix
!
!     -----------------------------------------------------------------
      integer, intent(in)    :: nrdm2
      integer, optional, intent(in)    :: iroot
      integer, optional, intent(in)    :: jroot
      real*8,  intent(inout) :: pv(*)
      character(len=*),optional,       intent(in)    :: checkpoint1
      character(len=*),optional,       intent(in)    :: checkpoint2
!     -----------------------------------------------------------------
      real*8, allocatable    :: tmp(:,:,:,:,:)
      real*8, allocatable    :: tmp2(:,:,:,:)
      real*8                 :: value, sgnjl,sgnik
      character(len=2300)    :: twoRDMfile
      integer                :: irootm1, jrootm1, jj, jjj, io
      integer                :: isorb, jsorb, ij, kl, ijkl
      integer                :: i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,nact,len_suffix
      character(len=1)       :: cJ,cJJ,cJJJ
      character(len=4)       :: suffix
      integer                :: lunit
      integer                :: iaka,iakb,ibia,ibka,ibkb,jala,jalb,jbla,jblb,jbja
      integer                :: iaaaa,ibbbb,ibaab,iabba,iabab,ibaba
      integer                :: iorba,iorbb,jorba,jorbb,korba,korbb,lorba,lorbb
      character(len=1997)    :: dmrg_file_prefix_save
!     -----------------------------------------------------------------

      twoRDMfile=""
      if (present(iroot).and.present(jroot)) then
        irootm1 = iroot - 1
        jrootm1 = jroot - 1

        if(irootm1 < 10)then
          J=irootm1
          cJ=CHAR(J+48)
          len_suffix = 2
          suffix(1:2) = "."//cJ
        else if(irootm1 < 100)then
          J=mod(irootm1,10)
          JJ=irootm1/10
          cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          len_suffix = 3
          suffix(1:3) = "."//cJJ//cJ
        else if(irootm1 < 1000)then
            J=mod(irootm1,10)
            JJ=mod(irootm1,100)/10
          JJJ=irootm1/100
            cJ=CHAR(J+48)
          cJJ=CHAR(JJ+48)
          cJJJ=CHAR(JJJ+48)
          len_suffix = 4
          suffix(1:4) = "."//cJJJ//cJJ//cJ
        end if

        if(irootm1 < 1000 .and. jrootm1 < 1000)then
          dmrg_file_prefix_save = dmrg_file%prefix
          dmrg_file%prefix      = ''
          call file_name_generator(jrootm1,"twoparticle.tdm.",16,suffix,len_suffix,twoRDMfile)
          dmrg_file%prefix      = dmrg_file_prefix_save
        else
          print *, "There are too many states (>999) in DMRG calculation"
          stop
        end if
      else if (present(checkpoint1).and.present(checkpoint2)) then
        twoRDMfile="twoparticle.tdm."//trim(checkpoint1)//"."//trim(checkpoint2)
      end if

#ifdef _DMRG_DEBUG_
      print *, 'reading 2-TDM from file: ', trim(twoRDMfile)
#endif

      lunit   =   140
      open(unit=lunit,file=twoRDMfile)
      read(lunit,*) nact
      !> Note that the full 2-TDM has been written to file
      !> with the permutations: pqrs == qpsr
      allocate(tmp(nact,nact,nact,nact,4)); tmp = 0
      do n = 1, 4
!       !> import alpha-alpha-alpha-alpha (n == 1)
!                 alpha-beta-beta-alpha   (n == 2)
!                 beta-alpha-alpha-beta   (n == 3)
!                 beta-beta-beta-beta     (n == 4)
!          blocks separately
        do i = 1, nact**4
          read(lunit,*,iostat=io) i1,j1,k1,l1,value
          if (is_iostat_end(io)) exit
          if (io > 0) stop 'problem reading 2-TDM from file'
          tmp(i1+1,j1+1,k1+1,l1+1,n) = value * 2.0d0
        end do
      end do

      close(lunit)

!#define BLUBB
!> blubb defines an alternative way to distribute the spin-summed 2-TDM directly into the lower triangle
!> activating this strategy requires to define BLUBB also in $MOLCAS/src/rassi/mktdm2.f

#ifndef BLUBB

      !> store TDM
#ifdef _DMRG_DEBUG_
      write(6,*) ' store sorted 2-TDM for lhs|rhs',iroot,jroot
#endif
      sgnjl = 1.0d0
      sgnik = 1.0d0
      do j=1,nact
        jorba = 2*j - 1
        jorbb = 2*j
        do i=1,nact
          iorba = 2*i - 1
          iorbb = 2*i
          do l=1,nact
            lorba = 2*l - 1
            lorbb = 2*l
            do k=1,nact
              korba = 2*k - 1
              korbb = 2*k

              !> skip upper triangular part
              if((i+(j-1)*nact) < (k+(l-1)*nact)) cycle
             !write(6,*) 'start distributing ',i,j,k,l

              if(j > l)then
                SGNJL=1.0D0
                jala=((jorba-1)*(jorba-2))/2+lorba
                jalb=((jorba-1)*(jorba-2))/2+lorbb
                jbla=((jorbb-1)*(jorbb-2))/2+lorba
                jblb=((jorbb-1)*(jorbb-2))/2+lorbb
              else if(j == l)then
                jbja=((jorbb-1)*(jorbb-2))/2+jorba
              else
                SGNJL=-1.0D0
                jala=((lorba-1)*(lorba-2))/2+jorba
                jalb=((lorbb-1)*(lorbb-2))/2+jorba
                jbla=((lorba-1)*(lorba-2))/2+jorbb
                jblb=((lorbb-1)*(lorbb-2))/2+jorbb
              end if
              if(i > k)then
                SGNIK=1.0D0
                iaka=((iorba-1)*(iorba-2))/2+korba
                iakb=((iorba-1)*(iorba-2))/2+korbb
                ibka=((iorbb-1)*(iorbb-2))/2+korba
                ibkb=((iorbb-1)*(iorbb-2))/2+korbb
              else if(i == k) then
                ibia=((iorbb-1)*(iorbb-2))/2+iorba
              else
                SGNIK=-1.0D0
                iaka=((korba-1)*(korba-2))/2+iorba
                iakb=((korbb-1)*(korbb-2))/2+iorba
                ibka=((korba-1)*(korba-2))/2+iorbb
                ibkb=((korbb-1)*(korbb-2))/2+iorbb
              end if
              if(i /= k)then
                if(j /= l)then
                  iaaaa     = iaka+nrdm2*(jala-1)
                  iabba     = iakb+nrdm2*(jalb-1)
                  ibaab     = ibka+nrdm2*(jbla-1)
                  ibbbb     = ibkb+nrdm2*(jblb-1)
                  pv(iaaaa) = sgnik*sgnjl*tmp(i,k,l,j,1)
                  pv(iabba) = sgnik*sgnjl*tmp(i,k,l,j,2)
                  pv(ibaab) = sgnik*sgnjl*tmp(i,k,l,j,3)
                  pv(ibbbb) = sgnik*sgnjl*tmp(i,k,l,j,4)
                else
                  iabab     = iakb+nrdm2*(jbja-1)
                  ibaab     = ibka+nrdm2*(jbja-1)
                  pv(iabab) = sgnik*(-tmp(i,k,l,j,2))
                  pv(ibaab) = sgnik*  tmp(i,k,l,j,3)
                end if
              else
                if(j /= l)then
                  ibaba     = ibia+nrdm2*(jalb-1)
                  ibaab     = ibia+nrdm2*(jbla-1)
                  pv(ibaba) = sgnjl*(-tmp(i,k,l,j,2))
                  pv(ibaab) = sgnjl*  tmp(i,k,l,j,3)
                else
                  ibaab     = ibia+nrdm2*(jbja-1)
                  pv(ibaab) = tmp(i,k,l,j,3)
                end if
              end if
            end do
          end do
        end do
      end do

      deallocate(tmp)

#else

      write(6,*) '2-tdm full'
      allocate(tmp2(nact,nact,nact,nact)); tmp2 = 0
      do i = 1, nact
      do j = 1, nact
      do k = 1, nact
      do l = 1, nact
        tmp2(i,j,k,l)  = (tmp(i,j,k,l,1) + tmp(i,j,k,l,2) &
                       +  tmp(i,j,k,l,3) + tmp(i,j,k,l,4) )

        write(6,*) ' i, j, k, l, value == ',i,j,k,l,tmp2(i,j,k,l)

      end do
      end do
      end do
      end do

      deallocate(tmp)

! #ifdef _DMRG_DEBUG_
!       !> set # e- by hand (here == 3)
!       call bro1(tmp2,nact,3,1)
! #endif
      do i = 1, nact
        do j = 1,nact
          do k = 1, nact
            do l = 1, nact
              if((i+(j-1)*nact) < (k+(l-1)*nact)) cycle
              ijkl    = ((i+(j-1)*nact)*((i+(j-1)*nact)-1))/2 + k + (l-1)*nact
              !> contraction over kl
              pv(ijkl)= tmp2(i,k,l,j)
            end do
          end do
        end do
      end do


      deallocate(tmp2)
#endif

      end subroutine dmrg_task_import_rdmY2
!**********************************************************************

      subroutine dmrg_task_import_rdmY3(tv,nrdm3,iroot,jroot)

!     -----------------------------------------------------------------
!     purpose: import the 3p-reduced density matrix generated by DMRG
!              for post-DMRG-SCF PT2 methods
!
!     -----------------------------------------------------------------
      integer, intent(in)    :: nrdm3
      integer, intent(in)    :: iroot
      integer, intent(in)    :: jroot
      real*8,  intent(inout) :: tv(nrdm3,nrdm3,nrdm3,nrdm3,nrdm3,nrdm3)
!     -----------------------------------------------------------------
      real*8                 :: value
      character(len=2300)    :: thrRDMfile
      integer                :: irootm1, jrootm1, jj, jjj, io
      integer                :: i,j,k,l,m,n,i1,j1,k1,l1,m1,n1,nact,len_suffix
      character(len=1)       :: cJ,cJJ,cJJJ
      character(len=4)       :: suffix
      integer                :: lunit
      character(len=1997)    :: dmrg_file_prefix_save
!     -----------------------------------------------------------------

      lunit   =   140
      irootm1 = iroot - 1
      jrootm1 = jroot - 1
      thrRDMfile=""

      if(jrootm1 < 10)then
        J=jrootm1
        cJ=CHAR(J+48)
        len_suffix = 2
        suffix(1:2) = "."//cJ
      else if(jrootm1 < 100)then
        J=mod(jrootm1,10)
        JJ=jrootm1/10
        cJ=CHAR(J+48)
        cJJ=CHAR(JJ+48)
        len_suffix = 3
        suffix(1:3) = "."//cJJ//cJ
      else if(jrootm1 < 1000)then
           J=mod(jrootm1,10)
          JJ=mod(jrootm1,100)/10
         JJJ=jrootm1/100
          cJ=CHAR(J+48)
         cJJ=CHAR(JJ+48)
        cJJJ=CHAR(JJJ+48)
        len_suffix = 4
        suffix(1:4) = "."//cJJJ//cJJ//cJ
      end if

      if(irootm1 < 1000 .and. jrootm1 < 1000)then
        dmrg_file_prefix_save = dmrg_file%prefix
        dmrg_file%prefix      = ''
        call file_name_generator(irootm1,"threeparticle.tdm.",18,suffix,len_suffix,thrRDMfile)
        dmrg_file%prefix      = dmrg_file_prefix_save
      else
        print *, "There are too many states (>999) in DMRG calculation"
        stop
      end if

#ifdef _DMRG_DEBUG_
      print *, 'reading 3-TDM from file: ', trim(thrRDMfile)
#endif
      nact = nrdm3
      open(unit=lunit,file=trim(thrRDMfile))
        do
          read(lunit,*,IOSTAT=io) i1, j1, k1, l1, m1, n1, value
          if (is_iostat_end(io)) exit
          if (io>0) stop 'problem reading 3-TDM from file'
          tv(i1,j1,k1,l1,m1,n1) = value
        end do
      close(lunit)

#ifdef _DMRG_DEBUG_
      print *, " 3-TDM was imported successfully"
#endif

      end subroutine dmrg_task_import_rdmY3
!**********************************************************************

      subroutine dmrg_task_import_rdm(dv,pv,tv,fv,nrdm,mrdm,ordm,prdm,iroot,rdm1,rdm2,rdm3,rdm4)

      real*8 , optional, dimension(*), intent(inout) :: dv
      real*8 , optional, dimension(*), intent(inout) :: pv
      real*8 , optional, dimension(*), intent(inout) :: tv
      real*8 , optional, dimension(*), intent(inout) :: fv
      integer, optional,               intent(in)    :: nrdm
      integer, optional,               intent(in)    :: mrdm
      integer, optional,               intent(in)    :: ordm
      integer, optional,               intent(in)    :: prdm
      integer, optional,               intent(in)    :: iroot
      logical, optional,               intent(in)    :: rdm1
      logical, optional,               intent(in)    :: rdm2
      logical, optional,               intent(in)    :: rdm3
      logical, optional,               intent(in)    :: rdm4

      if(present(rdm1) .and. rdm1) call dmrg_task_import_rdm1(dv=dv,nrdm1=nrdm,iroot=iroot)
      if(present(rdm2) .and. rdm2) call dmrg_task_import_rdm2(pv=pv,nrdm2=mrdm,iroot=iroot)
      if(present(rdm3) .and. rdm3) call dmrg_task_import_rdm3(tv=tv,nrdm3=ordm,iroot=iroot)
      !> note: 4-RDM is expected to be precomputed at this point - FIXME: check and possibly compute first at this stage
      !> stknecht - july 2015
      if(present(rdm4) .and. rdm4)then
        call dmrg_task_import_rdm4(fv=fv,nrdm4=prdm,iroot=iroot)
      end if

      end subroutine dmrg_task_import_rdm

      subroutine dmrg_task_import_rdm1(dv,nrdm1,iroot)

!     -----------------------------------------------------------------
!     purpose: import the 1p-reduced density matrix generated by DMRG
!              into the Molcas super-CI process, in order to do the orbital
!              optimization.
!     -----------------------------------------------------------------
      integer, intent(in)    :: nrdm1
      integer, intent(in)    :: iroot
      real*8,  intent(inout) :: dv(nrdm1)
!     -----------------------------------------------------------------
      real*8, allocatable    :: tmp(:,:)
      character(len=2300)    :: oneRDMfile
      !Switcher
      integer                :: lunit
      integer                :: i,j,ij,i1,j1,nact, irootm1
!     -----------------------------------------------------------------

#ifdef _DMRG_DEBUG_
      write(*,*)"importing 1-RDMs starting, iroot is ==> ",iroot
      call flush(6)
#endif

      lunit   =  140
      irootm1 = iroot - 1

! Since DMRG start from state - 0
      oneRDMfile=""
      if(irootm1.lt.1000)then
        call file_name_generator(irootm1,"oneparticle.rdm.",16,"",0,oneRDMfile)
      else
        write(6,*)"There are too many states (>999) in DMRG calculation"
        Stop
      end if

#ifdef _DMRG_DEBUG_
      write(6,*) ' file name for import is ==> ', oneRDMfile
      write(6,*) ' target state is         ==> ', irootm1
      call flush(6)
#endif

      open(unit=lunit,file=trim(oneRDMfile))
        read(lunit,*)nact
        if(.not.allocated(tmp)) allocate(tmp(nact,nact))
        tmp=0.0d0
        do i=1,nact
          do j=1,nact
            read(lunit,*)i1,j1,tmp(i1+1,j1+1)
          end do
        end do
      close(lunit)

      ij=0
      do i=1,nact
        do j=1,i
          ij=ij+1
          dv(ij)=tmp(i,j)
        end do
      end do

#ifdef _DMRG_DEBUG_
      write(6,*)"1-RDM imported successfully in DMRG-SCF for root ",iroot
#endif

      deallocate(tmp)

      end subroutine dmrg_task_import_rdm1
!**********************************************************************

      subroutine dmrg_task_import_rdm2(pv,nrdm2,iroot)

!     -----------------------------------------------------------------
!     purpose: import the 2p-reduced density matrix generated by DMRG
!              into the Molcas super-CI process, in order to do the orbital
!              optimization.
!     -----------------------------------------------------------------
      integer, intent(in)    :: nrdm2
      integer, intent(in)    :: iroot
      real*8,  intent(inout) :: pv(nrdm2)
!     -----------------------------------------------------------------
      real*8, allocatable    :: rdm2(:,:,:,:)
      real*8, allocatable    :: rdm2T(:,:,:,:)
      real*8, allocatable    :: rdm2P(:,:,:,:)
      real*8                 :: value
      character(len=2300)    :: twoRDMfile
      integer                :: i,j,k,l,ij,jk,kl,li,ijkl,nact,ioEOF,io
!     used for dalton
!     -----------------------------------------------------------------
      integer                :: kl_max, kl_min, kl_index
      real*8, allocatable    :: pv_dalton(:,:)
!     -----------------------------------------------------------------
      !> DALTON debug
!#include "../../../DALTON/include/priunit.h"
      integer icounter, irootm1
      integer lunit
!     -----------------------------------------------------------------

#ifdef _DMRG_DEBUG_
      write(*,*)"importing 2-RDMs start"
      call flush(6)
#endif

      lunit   =   140
      irootm1 = iroot - 1

! Since DMRG start from state - 0
      twoRDMfile=""
      if(irootm1.lt.1000)then
        call file_name_generator(irootm1,"twoparticle.rdm.",16,"",0,twoRDMfile)
      else
        write(6,*)"too many states (>999) in DMRG calculation"
        Stop
      end if

      open(unit=lunit,file=trim(twoRDMfile))
      read(lunit,*)nact
! Maquis saves only permutation-symmetry unique values - recover the full 2-RDM
      allocate(rdm2(nact,nact,nact,nact)); rdm2=0
      do
        read(lunit,*,iostat=io) ij,jk,kl,li,value
        if (is_iostat_end(io)) exit
        if (io > 0) stop 'problem reading 2-RDM from file'
        !> permutations: pqrs == rspq == qpsr == srqp
        rdm2(ij+1,jk+1,kl+1,li+1) = value
        rdm2(kl+1,li+1,ij+1,jk+1) = value
        rdm2(jk+1,ij+1,li+1,kl+1) = value
        rdm2(li+1,kl+1,jk+1,ij+1) = value
      end do
      close(lunit)

#ifdef _DMRG_DEBUG_
      print *, 'DEBUG RDM2...'
      do i=1,nact
         do j=1,nact
            do k=1,nact
               do l=1,nact
                  write(6,'(4i2,3x,G20.14)') &
                    i,j,k,l,rdm2(i,j,k,l)
               end do
            end do
         end do
      end do
      print *, 'END DEBUG RDM2...'
#endif
      !> rdm2 is in Mulliken notation. DALTON expects Dirac notation in 2-rdm
      !> j<=>k such that rmd2:molcas (ijkl) => rdm2:dalton (ikjl)
      if (dmrg_host_program_settings%dmrg_host_program(1:7) == 'dalton ')then
         allocate(rdm2P(nact,nact,nact,nact))
         rdm2p = 0
         do i=1,nact
            do j=1,nact
!              do k=min(i,j),nact
               do k=1,nact
                  do l=1,nact
                       rdm2P(i,j,k,l) = rdm2(i,j,k,l)
                  end do
               end do
            end do
         end do
         rdm2 = 0
         call dscal(nact**4,2.0d0,rdm2P,1)
         do i=1,nact
            do j=1,nact
               do k=1,nact
                  do l=1,nact
                     if(i==j .or. k==l)then
                       if(i==j .and. j==k .and. k==l)then !  nothing to do
                         rdm2(i,j,k,l) = rdm2P(i,j,k,l)
                       else                       !          j <--> k
                         rdm2(i,k,j,l) = rdm2P(i,j,k,l)
                       end if
                     else if(i==k .or. j==l)then  !          check for TDMs FIXME!!!
                       if(i==k .and. j==l)then
                         rdm2(i,j,l,k) = rdm2P(i,j,k,l)!     k <--> l
                       else
                         rdm2(i,j,l,k) = rdm2P(i,j,k,l)!     k <--> l
                       end if
                     else if(i==l .or. j==k)then  !          i <--> k
                       rdm2(k,j,i,l) = rdm2P(i,j,k,l)
                     else ! all four indices are different / i <--> k; k <--> l
                       rdm2(k,j,l,i) = rdm2P(i,j,k,l)
                     end if
                  end do
               end do
            end do
         end do
         ! DALTON debug
!        do i=1,nact
!           do j=1,nact
!              do k=1,nact
!                 do l=1,nact
!                    write(lupri,'(a,i2,a1,i2,a1,i2,a1,i2,a,f16.10)') &
!                      ' pv(',i,',',j,',',k,',',l,') = ',             &
!                        rdm2(i,j,k,l)
!                 end do
!              end do
!           end do
!        end do
         deallocate(rdm2P)
         allocate(pv_dalton((nact*(nact+1))/2,(nact*(nact+1))/2))
         pv_dalton = 0
         pv        = 0
         do l = 1, nact
            do k = 1, l
               kl_max = max(k,l)
               kl_min = min(k,l)
               kl_index = (kl_max*(kl_max-1))/2 + kl_min
!              distribute elements on lower triangle
               call dgetsp_util(nact,rdm2(1,1,k,l),pv_dalton(1,kl_index))
!              write(lupri,*)'kl_ST_index, pv(kl_ST_index)', kl_ST_index,pv_dalton(1,kl_ST_index)
            end do
         end do
         call dcopy(((nact*(nact+1))/2)**2,pv_dalton(1,1), 1, pv(1), 1)
         deallocate(pv_dalton)
         deallocate(rdm2)

      else ! host program is dalton

! The j and k packing, inorder to match Molcas - I
       allocate(rdm2T(nact,nact,nact,nact))
       rdm2T=0
        do i=1,nact
          do j=1,nact
            do k=1,j
              do l=1,nact
                if(j.eq.k)then
                  rdm2T(i,j,k,l)=rdm2(i,j,k,l)
                else
                  rdm2T(i,j,k,l)=rdm2(i,j,k,l)+rdm2(i,k,j,l)
                end if
              end do
            end do
          end do
        end do

        deallocate(rdm2)

!======================================================================
          ij=0
          ijkl=1
          do i=1,nact
            do j=1,i
              ij=ij+1
              kl=0
              do k=1,nact
                do l=1,k
                  kl=kl+1
                  if(ij.ge.kl)then
                    pv(ijkl)=rdm2T(i,k,l,j)
                    ijkl=ijkl+1
                  end if
                end do
              end do
            end do
          end do
          deallocate(rdm2T)
        end if ! host program is dalton

      end subroutine dmrg_task_import_rdm2
!**********************************************************************

      subroutine dmrg_task_import_rdm3(tv,nrdm3,iroot)

!     -----------------------------------------------------------------
!     purpose: import the 3p-reduced density matrix generated by DMRG
!              for post-DMRG-SCF PT2 methods
!
!     -----------------------------------------------------------------
      integer, intent(in)    :: nrdm3
      integer, intent(in)    :: iroot
      real*8,  intent(inout) :: tv(nrdm3,nrdm3,nrdm3,nrdm3,nrdm3,nrdm3)
!     -----------------------------------------------------------------
      character(len=2300)    :: thrRDMfile
      integer                :: i,j,k,l,m,n,i1,j1,k1,l1,m1,n1, irootm1
      integer                :: lunit
!     -----------------------------------------------------------------
      integer                :: read_error

      irootm1    = iroot - 1
      thrRDMfile =   ""
      lunit      = 140

!     temporary used, in order to aim at a correct root
! Since DMRG start from state - 0
      if(irootm1.lt.1000)then
        call file_name_generator(irootm1,"threeparticle.rdm.",18,"",0,thrRDMfile)
      else
        write(6,*)"There are too many states (>999) in DMRG calculation"
        stop
      end if

#ifdef _DMRG_DEBUG_
      write(6,*)thrRDMfile
#endif
      open(unit=lunit,file=trim(thrRDMfile))

      do while (.not.is_iostat_end(read_error))
        read (lunit,*,IOSTAT=read_error) i1,j1,k1,l1,m1,n1,tv(i1,j1,k1,l1,m1,n1)
        if (read_error.gt.0) then
          stop "problem reading 3-RDM from file"
        end if
        ! if read_error < 0 -> EOF == normal termination
      end do

      close(lunit)

#ifdef _DMRG_DEBUG_
      write(6,*)"3-RDM was imported successfully"
#endif


      end subroutine dmrg_task_import_rdm3
!**********************************************************************

      subroutine dmrg_task_import_rdm4(fv,nrdm4,iroot)

!     -----------------------------------------------------------------
!     purpose: import the 4p-reduced density matrix generated by DMRG
!              for post-DMRG-SCF PT2 methods
!
!     -----------------------------------------------------------------
      integer, intent(in)    :: nrdm4
      integer, intent(in)    :: iroot
      real*8,  intent(inout) :: fv(nrdm4*(nrdm4+1)*(nrdm4+2)*(nrdm4+3)/24,nrdm4,nrdm4,nrdm4,nrdm4)
!     -----------------------------------------------------------------
      real*8                 :: value
      character(len=2300)    :: fouRDMfile
      integer                :: i,j,k,l,m,n,new_norm,old_norm,i1,j1,k1,l1,m1,nact,irootm1
      integer                :: io,lunit
      integer                :: dim1
      integer                :: four_norm
!     -----------------------------------------------------------------

      lunit=140
      lunit=141
      irootm1    = iroot - 1
      fouRDMfile = ""
      nact       = nrdm4
      dim1       = nact*(nact+1)*(nact+2)*(nact+3)/24

!     temporary used, in order to aim at a correct root
! Since DMRG start from state - 0
      if(irootm1.lt.1000)then
        call file_name_generator(irootm1,"fourparticle.rdm.",17,"",0,fouRDMfile)
      else
        print *, "There are too many states (>999) in the DMRG calculation"
        stop
      end if

#ifdef _DMRG_DEBUG_
      print *, 'reading 4-RDM from file: ', trim(fouRDMfile)
#endif

      open(unit=lunit,file=trim(fouRDMfile),status='old',iostat=io)
      if(io /= 0)then
        print *, " QCMaquis> No 4-RDM file present for state ",iroot
        stop
      end if

      do
        read(lunit,*,iostat=io) old_norm,i1,j1,k1,l1,value
        if (is_iostat_end(io)) exit
        if (io > 0) stop 'problem reading 4-TDM from file'
        fv(old_norm,i1,j1,k1,l1) = value
#ifdef _DMRG_DEBUG_
        print '(i4,4i2,3x,1e20.14)', old_norm,i1,j1,k1,l1,fv(old_norm,i1,j1,k1,l1)
        call flush(6)
#endif
      end do
      close(lunit)

      end subroutine dmrg_task_import_rdm4
!**********************************************************************
      subroutine dmrg_task_prepare_hirdm_template(task,iroot,jroot,run_qcmaquis,compress_Mmax)
!     -------------
!     writes a template file for 4-RDM or trans-3RDM evaluation
!     the template file is likely to be adapted further to allow for partial evaluation
!     but it's not part of this routine
!     -------------
        character(len=4),intent(in) :: task ! 4rdm or 3rdm
        integer,intent(in) :: iroot ! root # starting with zero
        integer,intent(in),optional :: jroot ! ket root # starting with zero for tdm
        logical,intent(in),optional :: run_qcmaquis ! Run QCMaquis directly to evaluate the RDM if true, otherwise only write the template for later evaluation
        integer,intent(in),optional :: compress_Mmax ! M for MPS compression
        !---
        character(len=*),parameter :: runDMRG = '$MOLCAS/pytools/runDMRG.py'
        character(len=*),parameter :: mpstransform = 'mps_transform_pg'
        character(len=*),parameter :: extract4RDM = '$MOLCAS/qcmaquis/lib/python/pyeval/fourptrdm_perm.py'
        character(len=*),parameter :: extract3RDM = '$MOLCAS/qcmaquis/lib/python/pyeval/transition_threeptrdm.py'
        character(len=5000) :: runstring
        character(len=2300) :: checkpointname
        character(len=2300) :: bra_checkpointname
        character(len=2300) :: su2checkpointname
        character(len=2300) :: resultname
        character(len=255) :: templatename,inputname,currdir,prj

        character(len=10)  :: ms2str,irootstr,jrootstr
        integer :: nup, ndown

        ! command exit status
        integer :: status

        status = 0

        ! let's be consistent with other parts and ask for the host program
        if (dmrg_host_program_settings%dmrg_host_program(1:7) == 'molcas ')then
          ! calculate number of up and down electrons
          ! 2nup = nactel + ms; 2ndown = nactel - ms
          nup = (dmrg_state%nactel + dmrg_state%ms2)/2
          ndown = (dmrg_state%nactel - dmrg_state%ms2)/2
          !convert integers to strings
          write(ms2str,'(i10)') dmrg_state%ms2
          ms2str = adjustl(ms2str)
          write(irootstr,'(i10)') iroot
          irootstr = adjustl(irootstr)

          ! generate the 2u1 checkpoint/result name. checkpoint with highest ms=s will be used
          ! we don't use file_name_generator b/c it uses the full path which we want to avoid

          call getenv('Project',prj)
          checkpointname = trim(prj)//".checkpoint_state."//trim(irootstr)//"."//trim(ms2str)//"."//trim(ms2str)//".h5"
          su2checkpointname = trim(prj)//".checkpoint_state."//trim(irootstr)//".h5"
          resultname = trim(prj)//".results_state."//trim(irootstr)//"."//trim(ms2str)//"."//trim(ms2str)//".h5"

          if (task .eq. '4rdm') then
            ! generate the input name and set the template name
            inputname = 'meas-4rdm.'//trim(irootstr)//'.in'
            templatename = '$MOLCAS/template-files/template-dmrg-4rdm.maquis'

            ! the code below assumes that numbers will never exceed 4 digits

                write (runstring,'(a,i4,a,i4,a,i4,a,i4,a)') ' --replace="orbital_number=',dmrg_external%norb,   &
                            '" --replace="isup=',nup,'" --replace="isdown=',ndown,  &
                            '" --replace="isirrep=',0,'" --replace="saved_checkpoint='//trim(checkpointname)//  &
             '" --replace="saved_result='//trim(resultname)

             if (present(run_qcmaquis).and.run_qcmaquis) then
                 runstring=trim(runstring)//'" --executable=dmrg_meas '//trim(adjustl(templatename))
             else
                 runstring=trim(runstring)//'" --input-only='//inputname//' '//trim(adjustl(templatename))
             end if
            ! copy checkpoint to scratch and call mps_transform_pg to produce the 2u1 checkpoint file
            call getenv('CurrDir',currdir)
            call system('cp -r '//trim(currdir)//'/'//trim(su2checkpointname)//' $PWD/'//trim(su2checkpointname))

            !> make sure this works also on Mac OS X
            !call system(trim(mpstransform)//' '//trim(su2checkpointname),status)
            call system(runDMRG//' --transformSPIN '//'--executable='//trim(mpstransform)//&
                        ' --rhs='//trim(su2checkpointname),status)

            if (status.ne.0) then
              write (6,*) "mps_transform_pg failed with errorcode", status
              stop
            end if

            ! MPS compression:
            if (present(compress_Mmax).and.(compress_Mmax.gt.0)) then
              write (6,'(a,i6)') "Enabled MPS compression, new m=", compress_Mmax
              call dmrg_task_mps_compress(iroot,compress_Mmax)
            end if

            ! call runDMRG.py to produce the template file (or to evaluate the RDM directly)
            call system(runDMRG//trim(runstring),status)

            if (status.ne.0) then
              write (6,*) "runDMRG.py call failed with errorcode", status
              write (6,*) trim(runstring)
              stop
            end if

            ! if we ran, extract the rdm from the h5 file
            if (present(run_qcmaquis).and.run_qcmaquis) then
              call system(extract4RDM//" "//trim(resultname)//" "//irootstr,status)
              if (status.ne.0) then
                write (6,*) "fourptrdm.py call failed with errorcode", status
                stop
              end if
            end if

          else
            if (task .eq. '3tdm') then
              if (present(jroot)) then
                write(jrootstr,'(i10)') jroot
                jrootstr = adjustl(jrootstr)
              else
                stop '3-TDM QCMaquis template file requested but second root not provided!'
              end if
              inputname = 'meas-3tdm.'//trim(irootstr)//'.'//trim(jrootstr)//'.in'
              templatename = '$MOLCAS/template-files/template-dmrg-trans3rdm.maquis'
              bra_checkpointname = trim(prj)//".checkpoint_state."//trim(jrootstr)//"."//trim(ms2str)//"."//trim(ms2str)//".h5"
              write (runstring,'(a,i4,a,i4,a,i4,a,i4,a)') ' --replace="orbital_number=',dmrg_external%norb,   &
                              '" --replace="isup=',nup,'" --replace="isdown=',ndown,  &
                              '" --replace="isirrep=',0,'" --replace="saved_checkpoint='//trim(checkpointname)//  &
              '" --replace="saved_result='//trim(resultname)//'" --replace="bra_checkpoint='//trim(bra_checkpointname)

              if (present(run_qcmaquis).and.run_qcmaquis) then
                 runstring=trim(runstring)//'" --executable=dmrg_meas '//trim(adjustl(templatename))
              else
                 runstring=trim(runstring)//'" --input-only='//inputname//' '//trim(adjustl(templatename))
              end if

              ! warning! we don't call mps_transform_pg here, since we assume that this routine will be called together with the 4rdm routine, where it will be called for sure!
              ! same for the MPS compression!

              ! call runDMRG.py to produce the template file
              call system(runDMRG//trim(runstring),status)

              if (status.ne.0) then
                write (6,*) "runDMRG.py call failed with errorcode", status
                stop
              end if

            ! if we ran, extract the rdm from the h5 file
              if (present(run_qcmaquis).and.run_qcmaquis) then
                call system(extract3RDM//" "//trim(resultname)//" "//irootstr//" "//jrootstr,status)
                if (status.ne.0) then
                  write (6,*) "fourptrdm.py call failed with errorcode", status
                  stop
                end if
              end if
            else
              stop 'other rdms not supported yet'
            end if
          end if
        else
          ! not supported, do nothing
          write(6,*) "host programs other than MOLCAS not supported!"
        endif
      end subroutine dmrg_task_prepare_hirdm_template

!**********************************************************************
      subroutine dmrg_task_mps_compress(iroot,Mmax)
!     -------------------------
!     Performs an MPS compression
!     ------------------------
      integer, intent(in) :: iroot
      integer, intent(in) :: Mmax
      character(len=10) :: ms2str,irootstr,Mmaxstr
      character(len=2500) :: checkpointname
      character(len=5000) :: runstring
      character(len=255) :: prj
      ! command exit status
      integer :: status

      character(len=*),parameter :: runDMRG = '$MOLCAS/pytools/runDMRG.py'
      character(len=*),parameter :: mpscompress = 'mps_compress_2u1pg'

      status = 0

      !convert integers to strings
      write(ms2str,'(i10)') dmrg_state%ms2
      ms2str = adjustl(ms2str)
      write(irootstr,'(i10)') iroot
      irootstr = adjustl(irootstr)
      write(Mmaxstr,'(i10)') Mmax
      Mmaxstr = adjustl(Mmaxstr)

      call getenv('Project',prj)
      checkpointname = trim(prj)//".checkpoint_state."//trim(irootstr)//"."//trim(ms2str)//"."//trim(ms2str)//".h5"

      runstring = " --executable="//mpscompress//" --compress="//Mmaxstr//" --rhs="//checkpointname
      call system(runDMRG//trim(runstring),status)

      if (status.ne.0) then
        write (6,*) "runDMRG.py call failed with errorcode", status
        write (6,*) trim(runstring)
        stop
      end if

      end subroutine dmrg_task_mps_compress

!**********************************************************************

      subroutine dmrg_task_run_dmrg(              &
                                    key_DMRGonly, &
                                    iter          &
                                   )

!     -----------------------------------------------------------------
!     purpose: a. run a DMRG calculation
!              b. calculate 1-rdms and 2-rdms
!              Yingjin, Stefan - ETH Zurich
!     -----------------------------------------------------------------

      !> input
      logical                ::  Key_DMRGonly
      integer                ::  iter

      !> local scratch
      character(len=300)     :: input
      character(len=300)     :: pydriver
      integer                :: i,istart,iend,norb, iloop

!     -----------------------------------------------------------------

      !> set runtime environment
      call set_dmrg_runtime_environment(dmrg_setup%nproc)

      !> set lattice size
      norb = dmrg_external%norb

      !> prepare the input string for the python runscript
      !> -------------------------------------------------
      !> a. Molcas or Dalton
      if (dmrg_host_program_settings%dmrg_host_program(1:7) == 'dalton ')then
        input    = ' $DALTON/template-files/template-dmrg-su2.maquis'
        pydriver = ' $DALTON/pytools/runDMRG.py  '
      else if (dmrg_host_program_settings%dmrg_host_program(1:7) == 'molcas ')then
        input    = ' $MOLCAS/template-files/template-dmrg-su2.maquis'
        pydriver = ' $MOLCAS/pytools/runDMRG.py  '
      end if

      istart = 1; iend = 1
      if(dmrg_warmup%dofiedler .or. dmrg_warmup%docideas) iend = 2

      iloop = 1

      do i = istart, iend
        call run_dmrg_driver(                                   &
                             norb     = norb,                   &
                             pydriver = pydriver,               &
                             input    = input,                  &
                             dmrgci   = Key_DMRGonly,           &
                             iter     = iter,                   &
                             iloop    = iloop,                  &
                             fiedler  = dmrg_warmup%dofiedler,  &
                             cideas   = dmrg_warmup%docideas    &
                            )
        !> Fiedler ordering will only be determined once
        if(i == 1 .and. dmrg_warmup%dofiedler)then
          dmrg_warmup%dofiedler = .false.
          iloop = 2
        end if
        if(i == 1 .and. dmrg_warmup%docideas)then
          dmrg_warmup%docideas  = .false.
          iloop = -2
        end if
      end do


      end subroutine dmrg_task_run_dmrg

      subroutine run_dmrg_driver(              &
                                 norb,         &
                                 pydriver,     &
                                 input,        &
                                 dmrgci,       &
                                 iter,         &
                                 iloop,        &
                                 fiedler,      &
                                 cideas        &
                                 )

      !> input
      integer,            intent(in)  :: norb, iter, iloop
      logical,            intent(in)  :: dmrgci,fiedler,cideas
      character(len=300), intent(in)  :: pydriver,input

      !> local scratch
      character(len=2)                :: local_guess
      character(len=5)                :: state_tag
      character(len=8)                :: cFmt
      character(len=100)              :: maquis_norb
      character(len=100)              :: maquis_irrep
      character(len=100)              :: maquis_dmrg_model_symmetry
      character(len=100)              :: maquis_ele_total
      character(len=100)              :: maquis_spin
      character(len=100)              :: maquis_2pdm
      character(len=100)              :: maquis_1pdm
      character(len=100)              :: A100
      character(len=300)              :: maquis_orbital_occupation
      character(len=300)              :: project
      character(len=300)              :: input_local
      character(len=300)              :: put
      character(len=300)              :: maquis_num_ortho
      character(len=500)              :: keyword
      character(len=600)              :: ssinput, sslog, ssmeas, ssorder
      character(len=2300)             :: spd_name_sp,rdm_name_sp,rdm_name_dp
      character(len=2500)             :: maquis_chkp
      character(len=2500)             :: maquis_result
      character(len=540000)           :: maquis_name_ortho_states
      character(len=*),parameter      :: det2mps = 'det2mps_su2u1pg'

      integer                         :: nele_array(200), max_bond_dim_tmp
      integer                         :: i,j,nroot,maxroot, istatus
      integer                         :: ij,iroot,jroot,naimroot,ntmp,neletol
      integer                         :: nele_alpha,nele_beta,nele_mod
      integer                         :: itmp
      integer                         :: norbt_tmp, stringdim
      real*8                          :: dtmp, ttmp
      logical                         :: do_restart,dmrgerr_exists

      integer,allocatable             :: aimroot(:)
      character(len=2300),allocatable :: maquis_name_states(:)
      character(len=2300),allocatable :: maquis_name_results(:)
!     -----------------------------------------------------------------

      !> initialize some variables
      cFmt              = "(Ixxxxx)"
      nele_array(1:200) = 0
      do_restart        = .false.
      istatus           = 0
      dmrgerr_exists    = .false.
      maxroot           = dmrg_state%maxroot
      if(doSRDFTDMRG) do_restart = .true.

      !> set initial guess
      !> -----------------
      local_guess = "de"
      if(.not.fiedler)then
        !> a. HF guess
        if(sum(dmrg_orbital_space%initial_occ) > 0)then
          local_guess = "hf"
        end if
      end if

      !> a. prepare a local input
      input_local = ' ./template-dmrg-su2.maquis '
      call system('cp '//trim(input)//' '//trim(input_local))

#ifdef _DMRG_DEBUG_
      call system('cat '//trim(input_local))
      print *, 'size of qcmaquis input: ',size(dmrg_input%qcmaquis_input)
      print *, 'qcmaquis input: ',trim(dmrg_input%qcmaquis_input(1))
#endif
      open(898,file='./template-dmrg-su2.maquis',status='old',form='formatted',    &
           action='readwrite',access='append',position='append')
      stringdim = size(dmrg_input%qcmaquis_input)
      call prepare_local_input(                                &
                               898,                            &
                               dmrg_input%qcmaquis_input,      &
                               stringdim,                      &
                               local_guess,                    &
                               do_restart,                     &
                               iter,                           &
                               E_threshold,                    &
                               norb,                           &
                               fiedler.or.cideas               &
                              )
      close(898, status='keep')

#ifdef _DMRG_DEBUG_
      print *, 'qcmaquis input written to file! '
      call system('cat '//trim(input_local))
#endif

      !> c. number of orbitals/sites
      write(maquis_norb,'(a,i3,a)')    ' --replace="orbital_number='     ,norb,'"'
      !> target symmetry of state(s)
      write(maquis_irrep,'(a,i3,a)')   ' --replace="irrep_in_pointgroup=', dmrg_state%irefsm-1,'"'

      !> d. total number of e-
      write(maquis_ele_total,'(a,i3,a)')   ' --replace="electron_number_total= ', dmrg_state%nactel,'"'
      !> spin
      write(maquis_spin,'(a,i3,a)')       ' --replace="ms2=',dmrg_state%ms2,'"'

      !> e. initial guess
      if(local_guess /= 'hf')then
        write(maquis_orbital_occupation,'(a,a)')'   --replace=hf_occ_string=','delete-line'
      end if

      !> f. use point group symmetry in DMRG?
      write(maquis_dmrg_model_symmetry,'(a,a,a)')   ' --replace="su2u1pg=','su2u1pg','"' ! Add "pg" to ensure the mps2ci work

      !> g. check for CIonly (or fielder)
      if(dmrgci .or. fiedler.or.cideas)then
        write(maquis_2pdm,'(a,a)') '  --replace="active_2rdm=','delete-line"'
      else
        write(maquis_2pdm,'(a,a,a)') '  --replace="active_2rdm=','1','"'
      end if

      if(fiedler)then
        write(maquis_1pdm,'(a,a)') '  --replace="active_1rdm=','delete-line"'
      else
        write(maquis_1pdm,'(a,a,a)') '  --replace="active_1rdm=','1','"'
      end if

      !> h. integrals
      write(put,'(a)') ' --put="$PWD/FCIDUMP"'

      !> allocate some scratch
      allocate(aimroot(maxroot),maquis_name_states(maxroot),maquis_name_results(maxroot))
      aimroot = 0; maquis_name_states = ""; maquis_name_results = ""

      ! initialize state-averaged DMRG-SCF energy
      dmrg_energy%dmrg = 0.0d0
      !> set counter to zero
      naimroot = 0
      iroot    = 0

      !> long loop to run a DMRG calculation for each state
      compute_state: do

        iroot  = iroot + 1

        if(.not.fiedler .and. dmrg_orbital_ordering%fiedler_order(iroot) /= '')then
          !> state-specific orbital ordering (if available)
          open(898,file='./template-dmrg-su2.maquis',status='old',form='formatted',    &
               action='readwrite',access='append',position='append')
          write(898,'(a)') trim(dmrg_orbital_ordering%fiedler_order(iroot))
          close(898,status='keep')
        end if

        if(.not.fiedler .or. cideas)then
          !> use the user-defined occupation string for each state of interest rather than anything automatic
          norbt_tmp = 0
          do i=1,norb
            norbt_tmp = dmrg_orbital_space%initial_occ(i,naimroot+1)
          end do


          if((local_guess /= "hf" .and. norbt_tmp > 0).and..not.cideas)then
             print *, 'init_state = "hf" is a mandatory keyword in QCMaquis if SOCCupy is requested, I quit...'
             stop     'mandatory QCMaquis input keyword missing'
          end if

          !> prepare HF initial guess
          if(norbt_tmp > 0)then

            ij = 0

            do i=1,dmrg_symmetry%nirrep
              do j=1,dmrg_orbital_space%nash(i)
                ij = ij + 1
                nele_array(ij) = dmrg_orbital_space%initial_occ(ij,naimroot+1)
              end do
            end do

!           !> possibly rearrange the electron array
             maquis_orbital_occupation=""
             do i=1,norb
               write( cFmt(3:7) , '(i5)' ) 1
               A100=""
               write(A100,cFmt) nele_array(i)
               if(i.eq.norb)then
                 maquis_orbital_occupation=trim(maquis_orbital_occupation)//trim(A100)//" "
               else
                 maquis_orbital_occupation=trim(maquis_orbital_occupation)//trim(A100)//","
               end if
             end do

            nele_alpha = 0
            nele_beta  = 0
            do i=1,norb
              select case(nele_array(i))
                case(1)
                   nele_alpha = nele_alpha
                   nele_beta  = nele_beta
                case(2)
                   nele_alpha = nele_alpha
                   nele_beta  = nele_beta  + 1
                case(3)
                   nele_alpha = nele_alpha + 1
                   nele_beta  = nele_beta
                case(4)
                   nele_alpha = nele_alpha + 1
                   nele_beta  = nele_beta  + 1
                case default
                  stop 'unknown occupation: allowed: 1 (empty), 2 (spin down), 3 (spin up), 4 (docc)'
              end select
            end do

            write(maquis_orbital_occupation,'(a,a,a)')   &
            ' --replace="hf_occ_string='//trim(maquis_orbital_occupation)//'"'
          end if

        end if !> fiedler

        call file_name_generator(iroot-1,"oneparticle.spd.",16,"",0,spd_name_sp)
        call file_name_generator(iroot-1,"oneparticle.rdm.",16,"",0,rdm_name_sp)
        call file_name_generator(iroot-1,"twoparticle.rdm.",16,"",0,rdm_name_dp)

        call file_name_generator(iroot-1,"checkpoint_state.",17,".h5",3, maquis_name_states(iroot))
        call file_name_generator(iroot-1,"results_state.",   14,".h5",3,maquis_name_results(iroot))

#ifdef _DMRG_DEBUG_
        write(6,*)"spd_name_sp : ",trim(spd_name_sp)
        write(6,*)"rdm_name_sp : ",trim(rdm_name_sp)
        write(6,*)"rdm_name_dp : ",trim(rdm_name_dp)
        write(6,*)" maquis_name_states(iroot) : ", trim(maquis_name_states(iroot))
        write(6,*)"maquis_name_results(iroot) : ",trim(maquis_name_results(iroot))
#endif
        !> do_restart == .true. for iter > 1; should also be true if cideas was done in loop 1
        !> in order to restart from this MPS in loop 2 --> coded as -2; iloop == 2 for fiedler ordering
        do_restart = (do_restart .and. (iloop == 1)) .or. (iloop == -2)

        if(do_restart)then
          call system("rm -rf "//trim(spd_name_sp)//" "//trim(rdm_name_sp)//" "//trim(rdm_name_dp)//" ")
        else
          call system("rm -rf "//trim(maquis_name_states(iroot))//" "//&
                      trim(maquis_name_results(iroot))//" "//          &
                      trim(spd_name_sp)//" "//trim(rdm_name_sp)//" "// &
                      trim(rdm_name_dp)//" ")
        end if

        !> set string for lower states to orthogonalize to
        maquis_name_ortho_states=""
        do i=1,iroot-1
          if(i.eq.1)then
            maquis_name_ortho_states=trim(maquis_name_states(i))
          else
            maquis_name_ortho_states=trim(maquis_name_ortho_states)//" "//trim(maquis_name_states(i))
          end if
        end do

        write(maquis_num_ortho,'(a,i3,a)')   ' --replace="num_ortho_states=',iroot-1,'"'
        write(maquis_name_ortho_states,'(a,a,a)')   &
        ' --replace="states.ortho_series='//trim(maquis_name_ortho_states)//'"'
        write(maquis_chkp,'(a,a,a)')   ' --replace="saved_checkpoint='//trim(maquis_name_states(iroot))//'"'
        write(maquis_result,'(a,a,a)')   ' --replace="saved_result='//trim(maquis_name_results(iroot))//'"'

#ifdef _DMRG_DEBUG_
        print *, ' my pydriver command is:'
        print *, trim(pydriver)//" "//' --tmpfull=$PWD/tmp'//         &
                 trim(maquis_irrep)//trim(maquis_norb)//trim(maquis_ele_total)//trim(maquis_spin)//       &
                 trim(maquis_chkp)//trim(maquis_result)//                                                 &
                 trim(maquis_orbital_occupation)//trim(maquis_num_ortho)//trim(maquis_name_ortho_states)//&
                 trim(maquis_dmrg_model_symmetry)//trim(maquis_2pdm)//trim(put)//                         &
                 trim(maquis_1pdm)//                                                                      &
                 trim(input_local)
#endif

        !> delete old dmrg.err if it exists
        inquire(file="dmrg.err",exist=dmrgerr_exists)
        if (dmrgerr_exists) then
          open(unit=986,file="dmrg.err",status='old')
          close(unit=986,status='delete')
        end if

        !> run DMRG
        call system(trim(pydriver)//" "//' --tmpfull=$PWD/tmp'//      &
                    trim(maquis_irrep)//trim(maquis_norb)//trim(maquis_ele_total)//trim(maquis_spin)//       &
                    trim(maquis_chkp)//trim(maquis_result)//                                                 &
                    trim(maquis_orbital_occupation)//trim(maquis_num_ortho)//trim(maquis_name_ortho_states)//&
                    trim(maquis_dmrg_model_symmetry)//trim(maquis_2pdm)//trim(put)//                         &
                    trim(maquis_1pdm)//                                                                      &
                    trim(input_local),                                                                       &
                   istatus)
        if (istatus .ne. 0 ) then
          print *, "Error: DMRG run failed with exit code ", istatus
          call system("cat dmrg.err")
          stop ! TODO: there should be a possibility to pass the exit code back to MOLCAS and exiting from there
        end if

        !> calculate obsevables
        call compute_observables(pydriver)

        if(.not.(fiedler.or.cideas))then
          !> extract density matrices (1-RDM, 1-spin-RDM and possibly 2-RDM)
          call extract_rdm(                                     &
                               trim(maquis_name_results(iroot)),&
                           len_trim(maquis_name_results(iroot)),&
                                              trim(rdm_name_sp),&
                                          len_trim(rdm_name_sp),&
                                              trim(rdm_name_dp),&
                                          len_trim(rdm_name_dp),&
                                              trim(spd_name_sp),&
                                          len_trim(spd_name_sp),&
                                                 dmrgci )

          !> extract state-specific energy
          call system("$MOLCAS/pytools/energy_last.py "//trim(maquis_name_results(iroot))//&
                      " > dmrg_energy.out")
          open(unit=987,file="dmrg_energy.out",status='old')
          read(987,*)
          read(987,*) dtmp ! read energy
          read(987,*)
          read(987,*) itmp ! read # of sweeps
          read(987,*)
          read(987,*) ttmp ! read truncated weight
          close(unit=987,status='delete')

          !> ... and save state-specific and state-average energies
          dmrg_energy%dmrg_state_specific(iroot) = dtmp
          dmrg_energy%dmrg                       =                             &
          dmrg_energy%dmrg                       + dtmp*dmrg_state%weight(iroot)

          !> Special case: 1st DMRGSCF iteration -- dmrg is executed twice so the # of sweeps will be the same
          !> in this case, take the value of itmp as is
          if (dmrg_energy%num_sweeps_old(iroot).eq.itmp) then
             dmrg_energy%num_sweeps(iroot) = itmp
             !> maximum truncated weight
             dmrg_energy%max_truncW(iroot) = ttmp
          else
            !> otherwise, itmp accumulates # of sweeps, so subtract all the remaining sweeps
            dmrg_energy%num_sweeps(iroot) = itmp - dmrg_energy%num_sweeps_old(iroot)
            dmrg_energy%max_truncW(iroot) = ttmp - dmrg_energy%max_truncW_old(iroot)
          end if
          dmrg_energy%num_sweeps_old(iroot) = itmp
          dmrg_energy%max_truncW_old(iroot)  = ttmp

          !> extract orbital ordering
          call system("$MOLCAS/pytools/ordering.py "//trim(maquis_name_results(iroot))//" > ordering.out")

          !> rename some input/log files to keep in scratch
          if(dmrg_host_program_settings%dmrg_host_program(1:7) == 'molcas ')then

            state_tag(1:5) = ' '
            call getenv("Project",Project)
            call get_state_tag(iroot,state_tag,dmrg_file%offset)

            ssinput = Project(1:index(Project,' ')-1)//'.QCMaquis.state.'//trim(state_tag)//'.inp'
            dmrg_file%qcmaquis_parameter_file(iroot) = trim(ssinput)
            call system("cp dmrg-input "//trim(ssinput))
            sslog   = Project(1:index(Project,' ')-1)//'.QCMaquis-optimization.state.'//trim(state_tag)//'.log'
            call system("cp dmrg.out "//trim(sslog))
            ssmeas  = Project(1:index(Project,' ')-1)//'.QCMaquis-measurement.state.'//trim(state_tag)//'.log'
            call system("cp dmrg.RDM_OUT "//trim(ssmeas))
            ssorder = 'internal-orbital-ordering.state.'//trim(state_tag)
            call system("cp ordering.out "//trim(ssorder))
            !> save file name for later dump in external hdf5 file, e.g., rasscf.h5
            call file_name_generator_noprefix(iroot-1,"checkpoint_state.",17,".h5",3, &
                 dmrg_file%qcmaquis_checkpoint_file(iroot))
          end if

        else

          if(fiedler)then
            call system("$MOLCAS/pytools/fiedler.py "//trim(maquis_name_results(iroot))//" > fiedler.out")
            open(unit=987,file="fiedler.out",status='old')
            read(987,'(a,a,a)') dmrg_orbital_ordering%fiedler_order(iroot)
            close(987,status='delete')

            if(cideas)then
              open(987,file='dmrg-input',status='old',form='formatted',action='write',position='append')
              write(987,'(a)') trim(dmrg_orbital_ordering%fiedler_order(iroot))
              close(987,status='keep')
            end if

          end if

          if(cideas)then

             !> add true MAX_BOND_DIMENSION and CI-level for CI-DEAS to input file
             open(unit=987,file="dmrg-input",status='old',action='write',position='append')

             j = 0; keyword = 'MAX_BOND_DIMENSION'
             call find_qcmaquis_keyword(dmrg_input%qcmaquis_input,size(dmrg_input%qcmaquis_input),keyword,j)
             read(dmrg_input%qcmaquis_input(j+1),'(i20)') max_bond_dim_tmp
             write(987,'(a,a,i20)') 'max_bond_dimension',' = ',max_bond_dim_tmp

             j = 0; keyword = 'CI_LEVEL'
             call find_qcmaquis_keyword(dmrg_input%qcmaquis_input,size(dmrg_input%qcmaquis_input),keyword,j)
             if(j < 0)then
               write(987,'(a,a,a)') 'ci_level',' = "','1,2,3,4,5,6"'
             else
               write(987,'(a,a,a,a)') trim(dmrg_input%qcmaquis_input(j)),' = "',trim(dmrg_input%qcmaquis_input(j+1)),'"'
             end if

             close(987,status='keep')
             !> setup MPS
             call system(pydriver//" --executable="//det2mps//" --cideas")
          end if
        end if !> fiedler

        !> clean up
        call system("rm dmrg.out dmrg.RDM_OUT")
        call system("rm -rf tmp")

        !> increase counter for # computed roots and check whether we have computed all target states
        naimroot          = naimroot + 1
        aimroot(naimroot) = iroot

        if(naimroot == maxroot) exit compute_state

      end do compute_state

      deallocate(aimroot,maquis_name_states,maquis_name_results)

      end subroutine run_dmrg_driver

! ===========================================================================
! This subroutine generate the file name base on the prototype_name -Yingjin
!  Input  : iroot
!         : prototype_name
!         : len_prototype_name
!         : suffix
!         : len_suffix
!  Output : generated_name
! ===========================================================================

      subroutine file_name_generator(                                  &
                                     iroot,                            &
                                     prototype_name,                   &
                                     len_prototype_name,               &
                                     suffix,                           &
                                     len_suffix,                       &
                                     generated_name                    &
                                    )

     integer, intent(in)                           :: iroot
     integer, intent(in)                           :: len_prototype_name
     integer, intent(in)                           :: len_suffix
     character(len=len_prototype_name), intent(in) :: prototype_name
     character(len=len_suffix), intent(in)         :: suffix
     character(len=2300), intent(inout)            :: generated_name

     character cI,cII,cIII
     integer ln,lx,lp, i, ii, iii, iroot_local

     generated_name = ''

     ln=len_prototype_name
     lx=len_suffix
     lp=len_trim(dmrg_file%prefix)

     iroot_local = iroot + dmrg_file%offset

#ifdef _DMRG_DEBUG_
     print *, ' iroot is          ... ', iroot
     print *, ' iroot_local is    ... ', iroot_local
     print *, ' prototype_name is ... ', trim(prototype_name)
     print *, ' prefix is         ... ', trim(dmrg_file%prefix)
     print *, ' offset is         ... ', dmrg_file%offset
#endif

     if(iroot_local.lt.10)then
                 I=iroot_local
                cI=CHAR(I+48)
       if(len_suffix.eq.0)then
            generated_name(1:ln+3) = trim(prototype_name)//cI//"."//cI
       else
         generated_name(1:lp+ln+lx+1) = trim(dmrg_file%prefix)//trim(prototype_name)//cI//suffix
       end if
     else if(iroot_local.lt.100)then
                 I=mod(iroot_local,10)
                II=iroot_local/10
                cI=CHAR(I+48)
               cII=CHAR(II+48)
       if(len_suffix.eq.0)then
            generated_name(1:ln+5) = trim(prototype_name)//cII//cI//"."//cII//cI
       else
         generated_name(1:lp+ln+lx+2) = trim(dmrg_file%prefix)//trim(prototype_name)//cII//cI//suffix
       end if
     else if(iroot_local.lt.1000)then
                 I=mod(iroot_local,10)
                II=mod(iroot_local,100)/10
               III=iroot_local/100
                cI=CHAR(I+48)
               cII=CHAR(II+48)
              cIII=CHAR(III+48)
       if(len_suffix.eq.0)then
            generated_name(1:ln+7) = trim(prototype_name)//cIII//cII//cI//"."//cIII//cII//cI
       else
         generated_name(1:lp+ln+lx+3) = trim(dmrg_file%prefix)//trim(prototype_name)//cIII//cII//cI//suffix
       end if
     else
       print *, ' error: iroot_local is too large (>= 1000): ',iroot_local
       stop 88
     end if

  end subroutine file_name_generator


! 01-12-2016 Leon: File name generator without prefix
! implemented in an extra routine for keeping the backward compatibility
! this is possibly NOT threadsafe, however should there be more than one
! thread of QCMaquis interface running at once?
  subroutine file_name_generator_noprefix(                        &
                                iroot,                            &
                                prototype_name,                   &
                                len_prototype_name,               &
                                suffix,                           &
                                len_suffix,                       &
                                generated_name                    &
                                )

     integer, intent(in)                           :: iroot
     integer, intent(in)                           :: len_prototype_name
     integer, intent(in)                           :: len_suffix
     character(len=len_prototype_name), intent(in) :: prototype_name
     character(len=len_suffix), intent(in)         :: suffix
     character(len=256), intent(inout)            :: generated_name

     character(len=2300)              :: tmp_genname
     character(len=1997)              :: tmp_prefix

     tmp_prefix = dmrg_file%prefix
     call getenv("Project", dmrg_file%prefix)
     dmrg_file%prefix = trim(dmrg_file%prefix)//"."
     call file_name_generator(                                    &
                                iroot,                            &
                                prototype_name,                   &
                                len_prototype_name,               &
                                suffix,                           &
                                len_suffix,                       &
                                tmp_genname                       &
                                )
     dmrg_file%prefix = tmp_prefix
     generated_name = trim(tmp_genname)
  end subroutine file_name_generator_noprefix
! 19-09-2016 Leon: commented bro1 b/c of version mismatch with one in NEVPT2 code
!       subroutine bro1(taa,nact,nele,metat)
!       implicit none
!       real*8, intent(in)  :: taa(metat,nact,nact,nact,nact)
!       integer, intent(in) :: nact, nele, metat
! !     scratch
!       integer             :: i, j, ll, istate
!       real*8, allocatable :: dal(:,:,:)
!
!       allocate(dal(metat,nact,nact)); dal = 0
!       do i=1,nact
!          do j=1,nact
!             do ll=1,nact
!                do istate=1,metat
!                   dal(istate,i,j)=dal(istate,i,j)+1.0d0*taa(istate,i,ll,ll,j)
!                enddo
!             enddo
!             do istate=1,metat
!                dal(istate,i,j)=dal(istate,i,j)/(nele-1)
!             enddo
!          enddo
!       enddo
!
!       print *, 'blubb DEBUG t-rho1'
!       do istate=1,metat
!       do i=1,nact
!         do j=1,nact
!           if(abs(dal(istate,i,j)) > 1.0d-16)&
!           print '(i1,1i2,3x,d19.12)',&
!           i,j,dal(istate,i,j)
!         enddo
!       enddo
!
!       enddo ! istate
!       deallocate(dal)
!
!       end subroutine bro1
! !------------------------------------------------------------------------


end module qcmaquis_interface_main


  integer function four_norm(ia,ic,ie,ig)

    implicit none

    integer, intent(in)   :: ia, ic, ie, ig

    integer, dimension(4) :: nv
    integer               :: i,j,jmax,imax,ndum

    nv(1)=ia
    nv(2)=ic
    nv(3)=ie
    nv(4)=ig

    do i = 1,4

     imax = nv(i)
     jmax = i

      do j=i+1,4
        if(nv(j) > imax)then
          imax = nv(j)
          jmax = j
        end if
      end do

      if(jmax /= i) call swapint(nv(i),nv(jmax))

    end do

    four_norm=(nv(1)-1)*nv(1)*(nv(1)+1)*(nv(1)+2)/24+ &
              (nv(2)-1)*nv(2)*(nv(2)+1)/6+            &
               nv(3)*(nv(3)-1)/2+nv(4)

    end function four_norm

! ----------------------------------------------------------------------
    subroutine swapint(a, b)
      implicit none
      INTEGER, INTENT(IN OUT) :: a, b
      INTEGER :: temp
      temp = a ; a = b ; b = temp
    end subroutine swapint

! ----------------------------------------------------------------------

