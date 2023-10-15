
! ---

BEGIN_PROVIDER [double precision, ao_Int_elecE, (ao_num, ao_num)]

  BEGIN_DOC
  !
  ! Electronic interaction with the external Field E in the |AO| basis.
  !
  END_DOC

  implicit none
  integer          :: i, j, n, l, dim1, power_A(3), power_B(3)
  double precision :: overlap_x, overlap_y, overlap_z, overlap1, overlap2
  double precision :: alpha, beta, c
  double precision :: A_center(3), B_center(3)

  dim1 = 100

  ! -- Dummy call to provide everything
  A_center(:) = 0.d0
  B_center(:) = 1.d0
  alpha = 1.d0
  beta  = .1d0
  power_A = 1
  power_B = 0
  call overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap1, dim1)


  !$OMP PARALLEL DO SCHEDULE(GUIDED)                                           &
  !$OMP DEFAULT(NONE)                                                          &
  !$OMP PRIVATE(i, j, n, l, alpha, beta, A_center, B_center, power_A, power_B, &
  !$OMP         c, overlap_x, overlap_y, overlap_z, overlap1, overlap2)        &
  !$OMP SHARED(nucl_coord, ao_power, ao_prim_num, ao_num, ao_nucl, dim1,       &
  !$OMP        ao_coef_normalized_ordered_transp, ao_expo_ordered_transp,      &
  !$OMP        ao_Int_elecE)
  do j = 1, ao_num

    A_center(1) = nucl_coord(ao_nucl(j),1)
    A_center(2) = nucl_coord(ao_nucl(j),2)
    A_center(3) = nucl_coord(ao_nucl(j),3)

    power_A(1)  = ao_power(j,1)
    power_A(2)  = ao_power(j,2)
    power_A(3)  = ao_power(j,3)

    do i = 1, ao_num

      B_center(1) = nucl_coord(ao_nucl(i),1)
      B_center(2) = nucl_coord(ao_nucl(i),2)
      B_center(3) = nucl_coord(ao_nucl(i),3)

      power_B(1)  = ao_power(i,1)
      power_B(2)  = ao_power(i,2)
      power_B(3)  = ao_power(i,3)

      ao_Int_elecE(i,j) = 0.d0

      do n = 1, ao_prim_num(j)
        alpha = ao_expo_ordered_transp(n,j)

        do l = 1, ao_prim_num(i)
          beta = ao_expo_ordered_transp(l,i)

          c = ao_coef_normalized_ordered_transp(n,j) * ao_coef_normalized_ordered_transp(l,i)

          call overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap1, dim1)

          power_A(3) = power_A(3) + 1
          call overlap_gaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, overlap_x, overlap_y, overlap_z, overlap2, dim1)
          power_A(3) = power_A(3) - 1

          ao_Int_elecE(i,j) = ao_Int_elecE(i,j) + c * (A_center(3) * overlap1 + overlap2)
        enddo
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO

  ao_Int_elecE = ao_Int_elecE * elec_field

END_PROVIDER

! ---


