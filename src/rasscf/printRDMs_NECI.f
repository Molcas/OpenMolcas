      module print_RDMs_NECI_format
      use rasscf_data, only : NAC, NACPAR
      use index_symmetry, only : two_el_idx, one_el_idx
      private
      public:: printRDMs_NECI
      contains

      subroutine printRDMs_NECI(DMAT,NAC,PMAT,PA,NACPAR)
      implicit none
      real, parameter :: thrsh = 1.0d-12
      integer :: NAC, NACPAR, i, idx1(2)
      real*8  :: DMAT(NAC), PMAT(NACPAR), PA(NACPAR)
      integer :: p, q, r, s, idx(4)

        Write(6,*) ' In printRDMs_NECI:'

        do i = 1, (NACPAR*(NACPAR+1)/2)
          idx = two_el_idx(i, p, q, r, s)
          if(r /= s) then
             if(abs(PMAT(i)+PA(i)).gt.thrsh)
     &           write(6,'(1X,4I5,F20.12)')
     &           p, q, r, s, (PMAT(i)+PA(i))
             if (abs(PMAT(i)-PA(i)).gt.thrsh)
     &           write(6,'(1X,4I5,F20.12)')
     &           p, q, s, r, (PMAT(i)-PA(i))
          else
             if (abs(PMAT(i)*2.0d0).gt.thrsh)
     &          write(6,'(1X,4I5,F20.12)') idx, PMAT(i)*2.0d0
          end if
        end do

        do i = 1, (NAC*(NAC+1)/2)
          idx1 = one_el_idx(i, p, q)
          if(abs(DMAT(i)).gt.thrsh)
     &       write(6,'(I6,F20.12)') i, dmat(i)
!     &       write(6,'(1X,4I5,F20.12)')idx1,0,0,DMAT(i)
        end do

        CALL TRIPRT('Averaged one-body density matrix, D, in RASSCF',
     &              ' ',DMAT,NAC)
        CALL TRIPRT('Averaged two-body density matrix, P',
     &              ' ',PMAT,NACPAR)
        CALL TRIPRT('Averaged antisym 2-body density matrix PA RASSCF',
     &              ' ',PA , NACPAR)

      end subroutine printRDMs_NECI

      end module print_RDMs_NECI_format
