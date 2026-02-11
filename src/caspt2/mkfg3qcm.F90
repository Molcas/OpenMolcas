!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Stefano Battaglia                                *
!***********************************************************************
#ifdef _DMRG_

subroutine mkfg3qcm(IFF, G1, F1, G2, F2, G3, F3, idxG3)

   use Symmetry_Info, only: Mul
   use stdalloc, only: mma_allocate, mma_deallocate
   use qcmaquis_interface
   use definitions, only: wp, iwp, i1, u6
   use printLevel, only: verbose
   use gugx, only: SGS
   use caspt2_global, only: iPrGlb
   use caspt2_module, only: nAshT, jState, EPSA
      use pt2_guga, only: nG3

   implicit none


   Integer(kind=iwp), intent(in)  :: IFF
   Real(kind=wp), intent(out) :: G1(nasht, nasht), G2(nasht, nasht, nasht, nasht)
   Real(kind=wp), intent(out) :: F1(nasht, nasht), F2(nasht, nasht, nasht, nasht)
   Real(kind=wp), intent(out) :: G3(*), F3(*)
   Integer(kind=i1), intent(in) :: idxG3(6, *)
   Real(kind=wp), allocatable :: G3tmp(:, :, :, :, :, :), TG3tmp(:, :, :, :, :, :)
   ! Real(kind=wp), allocatable :: G4(:,:,:,:,:)
   Integer(kind=iwp) :: t, u, v, w, x, y, z
   Integer(kind=iwp) :: i

   ! This might be memory hungry
   call mma_allocate(G3tmp,nasht,nasht,nasht,nasht,nasht,nasht)
   call mma_allocate(TG3tmp,nasht,nasht,nasht,nasht,nasht,nasht)
   ! call mma_allocate(G4,n4,nasht,nasht,nasht,nasht,Label='G4')

   if (iPrGlb >= verbose) then
      write(u6,*) ">QCMaquis: Computing RDMs for STATE: ", JSTATE-1
   end if

   ! Remeasure RDM from rotated MPS
   call qcmaquis_interface_get_1rdm_full(G1)
   call qcmaquis_interface_get_2rdm_full(G2)
   call qcmaquis_interface_get_3rdm_full(G3tmp)

   ! Fetch rotated RDM from unrotated MPS
   ! call qcmaquis_interface_read_rdm_full(int(JSTATE-1, c_int), &
   !   int(JSTATE-1, c_int), G1, int(1, c_int))
   ! call qcmaquis_interface_read_rdm_full(int(JSTATE-1, c_int), &
   !   int(JSTATE-1, c_int), G2, int(2, c_int))
   ! call qcmaquis_interface_read_rdm_full(int(JSTATE-1, c_int), &
   !   int(JSTATE-1, c_int), G3tmp, int(3, c_int))

   ! TODO: This should be removed in the future
   TG3tmp = 0.0_wp  ! Initialize TG3tmp to zero
   call qcmaquis_interface_get_fock_contracted_4rdm_full(EPSA, 0)
   call qcmaquis_interface_read_fock_contracted_4rdm(TG3tmp, logical(.false., c_bool))

   ! call qcmaquis_interface_get_fock_contracted_4rdm_full(TG3tmp, epsa, CompressMPS)

   ! number of elements in the contracted 4-index of the 4-RDM
   ! n4 = (nasht*(nasht+1)*(nasht+2)*(nasht+3)/24)
   ! call mma_allocate(G4,n4,nasht,nasht,nasht,nasht,Label='G4')
   ! call qcmaquis_interface_get_4rdm_full(G4)

   if (iff > 0) then
      do t = 1, nasht
         do u = 1, t
            if (Mul(SGS%ism(t), SGS%ism(u)) == 1) then
               do w = 1, nasht
                  F1(t, u) = F1(t, u) + G2(t, u, w, w)*epsa(w)
               end do
               F1(u, t) = F1(t, u)
            end if
         end do
      end do
   end if

   if (iff > 0) then
      do t = 1, nasht
         do u = 1, nasht
            do v = 1, nasht
               do x = 1, nasht
                  if (Mul(SGS%ism(x), Mul(SGS%ism(v), Mul(SGS%ism(u), SGS%ism(t)))) == 1) then
                     do w = 1, nasht
                        F2(t, u, v, x) = F2(t, u, v, x) + G3tmp(t, v, w, u, x, w)*epsa(w)
                     end do
                  end if
               end do
            end do
         end do
      end do
   end if

   do i = 1, ng3
      t = idxG3(1, i)
      u = idxG3(2, i)
      v = idxG3(3, i)
      x = idxG3(4, i)
      y = idxG3(5, i)
      z = idxG3(6, i)

      G3(i) = G3tmp(t, v, y, u, x, z)
      ! F_tvyuxz = <| e_{tv,yu,xz} E_ww f_ww |> - <| e_{tv,yu,xz} |> (f_uu + f_xx + f_zz)
      ! Compute first term using transition RDM between <| and |'> = E_ww f_ww |>: <| e_{tv,yu,xz} |'>
      F3(i) = TG3tmp(t, v, y, u, x, z) - G3(i)*(epsa(u) + epsa(x) + epsa(z))

      ! do w = 1,nasht
      !   val = get_p4_element(G4,w,t,v,y,w,u,x,z,nasht)
      !   F3(i) = F3(i) + val * epsa(w)
      ! end do
   end do

   call mma_deallocate(G3tmp)
   call mma_deallocate(TG3tmp)
   ! call mma_deallocate(G4)
end subroutine mkfg3qcm

#elif defined (NAGFOR)
! Some compilers do not like empty files
subroutine empty_mkfg3qcm()
end subroutine empty_mkfg3qcm
#endif

