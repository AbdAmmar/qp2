
! ---

program compute_no_gpu

  implicit none

  print *, ' j2e_type = ', j2e_type
  print *, ' j1e_type = ', j1e_type
  print *, ' env_type = ', env_type

  my_grid_becke  = .True.
  PROVIDE tc_grid1_a tc_grid1_r 
  my_n_pt_r_grid = tc_grid1_r   
  my_n_pt_a_grid = tc_grid1_a
  touch my_grid_becke my_n_pt_r_grid my_n_pt_a_grid

  my_extra_grid_becke  = .True. 
  PROVIDE tc_grid2_a tc_grid2_r 
  my_n_pt_r_extra_grid = tc_grid2_r
  my_n_pt_a_extra_grid = tc_grid2_a
  touch my_extra_grid_becke my_n_pt_r_extra_grid my_n_pt_a_extra_grid

  call write_int(6, my_n_pt_r_grid, 'radial  external grid over')
  call write_int(6, my_n_pt_a_grid, 'angular external grid over')

  call write_int(6, my_n_pt_r_extra_grid, 'radial  internal grid over')
  call write_int(6, my_n_pt_a_extra_grid, 'angular internal grid over')

  call main()
    
end

! ---

subroutine main()

  use cutc_module

  implicit none

  integer                       :: i, j, k, l, ipoint
  double precision              :: time0, time1
  double precision              :: tt0, tt1
  double precision              :: acc_thr, err_tot, nrm_tot, err_loc
  double precision              :: noL_0e_gpu(1)
  double precision, allocatable :: int2_grad1_u12_ao(:,:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)
  double precision, allocatable :: int2_grad1_u12_bimo_t(:,:,:,:)
  double precision, allocatable :: noL_1e_gpu(:,:)
  double precision, allocatable :: noL_2e_gpu(:,:,:,:)

  PROVIDE mo_l_coef mo_r_coef
  PROVIDE mos_l_in_r_array_transp mos_r_in_r_array_transp


  call wall_time(time0)
  print*, ' start compute_no_gpu'



  allocate(int2_grad1_u12_ao(ao_num,ao_num,n_points_final_grid,3))
  print*, ' Reading int2_grad1_u12_ao from ', trim(ezfio_filename) // '/work/int2_grad1_u12_ao'
  call wall_time(tt0)
  open(unit=11, form="unformatted", file=trim(ezfio_filename)//'/work/int2_grad1_u12_ao', action="read")
    read(11) int2_grad1_u12_ao
  close(11)
  call wall_time(tt1)
  write(*,"(A,2X,F15.7)") ' wall time for reading (sec) = ', (tt1 - tt0)

  allocate(tmp(mo_num,mo_num,n_points_final_grid,3))
  allocate(int2_grad1_u12_bimo_t(n_points_final_grid,3,mo_num,mo_num))

  call wall_time(tt0)
  !$OMP PARALLEL         &
  !$OMP DEFAULT (NONE)   &
  !$OMP PRIVATE (ipoint) & 
  !$OMP SHARED (ao_num, mo_num, n_points_final_grid, int2_grad1_u12_ao, tmp)
  !$OMP DO SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
    call ao_to_mo_bi_ortho(int2_grad1_u12_ao(1,1,ipoint,1), ao_num, tmp(1,1,ipoint,1), mo_num)
    call ao_to_mo_bi_ortho(int2_grad1_u12_ao(1,1,ipoint,2), ao_num, tmp(1,1,ipoint,2), mo_num)
    call ao_to_mo_bi_ortho(int2_grad1_u12_ao(1,1,ipoint,3), ao_num, tmp(1,1,ipoint,3), mo_num)
  enddo
  !$OMP END DO
  !$OMP END PARALLEL

  deallocate(int2_grad1_u12_ao)

  !$OMP PARALLEL               &
  !$OMP DEFAULT (NONE)         &
  !$OMP PRIVATE (i, j, ipoint) & 
  !$OMP SHARED (mo_num, n_points_final_grid, tmp, int2_grad1_u12_bimo_t)
  !$OMP DO COLLAPSE(2) SCHEDULE (dynamic)
  do ipoint = 1, n_points_final_grid
    do i = 1, mo_num
      do j = 1, mo_num
        int2_grad1_u12_bimo_t(ipoint,1,j,i) = tmp(j,i,ipoint,1)
        int2_grad1_u12_bimo_t(ipoint,2,j,i) = tmp(j,i,ipoint,2)
        int2_grad1_u12_bimo_t(ipoint,3,j,i) = tmp(j,i,ipoint,3)
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  call wall_time(tt1)
  write(*,"(A,2X,F15.7)") ' wall time for 3e-tensor (sec) = ', (tt1 - tt0)

  deallocate(tmp)

  allocate(noL_2e_gpu(mo_num,mo_num,mo_num,mo_num))
  allocate(noL_1e_gpu(mo_num,mo_num))

  call cutc_no(n_points_final_grid, mo_num, elec_alpha_num, elec_beta_num, &
               final_weight_at_r_vector(1),                                &
               mos_l_in_r_array_transp(1,1), mos_r_in_r_array_transp(1,1), &
               int2_grad1_u12_bimo_t(1,1,1,1), noL_2e_gpu(1,1,1,1), noL_1e_gpu(1,1), noL_0e_gpu(1))

  deallocate(int2_grad1_u12_bimo_t)
  deallocate(noL_2e_gpu)
  deallocate(noL_1e_gpu)

  call wall_time(time1)
  write(*,"(A,2X,F15.7)") ' wall time for compute_no_gpu (sec) = ', (time1 - time0)

  return

end

! ---

