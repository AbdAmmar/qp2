! ---

BEGIN_PROVIDER [double precision, ao_overlap_cgtos, (ao_num, ao_num)]

  implicit none

  integer          :: i, j, m, n, l, ii, jj, dim1, power_A(3), power_B(3)
  double precision :: c, overlap
  double precision :: phiA, phiB
  complex*16       :: alpha, alpha_inv, A_center(3)
  complex*16       :: beta, beta_inv, B_center(3)
  complex*16       :: C1, C2
  complex*16       :: overlap1, overlap_x1, overlap_y1, overlap_z1
  complex*16       :: overlap2, overlap_x2, overlap_y2, overlap_z2

  ao_overlap_cgtos   = 0.d0

  dim1 = 100

 !$OMP PARALLEL DO SCHEDULE(GUIDED)                                     &
 !$OMP DEFAULT(NONE)                                                    &
 !$OMP PRIVATE(i, j, m, n, l, ii, jj, c, C1, C2,                        &
 !$OMP         alpha, alpha_inv, A_center, power_A, phiA,               &
 !$OMP         beta, beta_inv, B_center, power_B, phiB, overlap,        &
 !$OMP         overlap_x1, overlap_y1, overlap_z1, overlap1,            &
 !$OMP         overlap_x2, overlap_y2, overlap_z2, overlap2)            &
 !$OMP SHARED(nucl_coord, ao_power, ao_prim_num, ao_num, ao_nucl, dim1, &
 !$OMP        ao_coef_cgtos_norm_ord_transp, ao_expo_cgtos_ord_transp,  &
 !$OMP        ao_expo_pw_ord_transp, ao_expo_phase_ord_transp,          &
 !$OMP        ao_overlap_cgtos)

  do j = 1, ao_num

    jj = ao_nucl(j)
    power_A(1) = ao_power(j,1)
    power_A(2) = ao_power(j,2)
    power_A(3) = ao_power(j,3)

    do i = 1, ao_num

      ii = ao_nucl(i)
      power_B(1) = ao_power(i,1)
      power_B(2) = ao_power(i,2)
      power_B(3) = ao_power(i,3)

      do n = 1, ao_prim_num(j)

        alpha = ao_expo_cgtos_ord_transp(n,j)
        alpha_inv = (1.d0, 0.d0) / alpha
        phiA = ao_expo_phase_ord_transp(4,n,j)
        do m = 1, 3
          A_center(m) = dcmplx(nucl_coord(jj,m))
        enddo

        do l = 1, ao_prim_num(i)

          beta = ao_expo_cgtos_ord_transp(l,i)
          beta_inv = (1.d0, 0.d0) / beta
          phiB = ao_expo_phase_ord_transp(4,l,i)
          do m = 1, 3
            B_center(m) = dcmplx(nucl_coord(ii,m))
          enddo

          c = ao_coef_cgtos_norm_ord_transp(n,j) * ao_coef_cgtos_norm_ord_transp(l,i)

          C1 = zexp((0.d0, 1.d0) * (-phiA - phiB))
          C2 = zexp((0.d0, 1.d0) * ( phiA - phiB))

          call overlap_cgaussian_xyz(A_center, B_center, alpha, beta, power_A, power_B, &
                                     A_center, B_center, overlap_x1, overlap_y1, overlap_z1, overlap1, dim1)

          call overlap_cgaussian_xyz(conjg(A_center), B_center, conjg(alpha), beta, power_A, power_B, &
                                     conjg(A_center), B_center, overlap_x2, overlap_y2, overlap_z2, overlap2, dim1)

          overlap = 0.5d0 * real(C1 * overlap1 + C2 * overlap2)

          ao_overlap_cgtos(i,j) = ao_overlap_cgtos(i,j) + c * overlap

          if(isnan(ao_overlap_cgtos(i,j))) then
            print*,'i, j', i, j
            print*,'l, n', l, n
            print*,'c, overlap', c, overlap
            stop
          endif
        enddo
      enddo
    enddo
  enddo
 !$OMP END PARALLEL DO

END_PROVIDER

! ---


