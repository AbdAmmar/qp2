
BEGIN_PROVIDER [logical, use_cgtos]

  implicit none

  BEGIN_DOC
  ! If true, use cgtos for AO integrals
  END_DOC

  logical                        :: has
  PROVIDE ezfio_filename
  use_cgtos = .False.
  if (mpi_master) then
    call ezfio_has_ao_basis_use_cgtos(has)
    if (has) then
!      write(6,'(A)') '.. >>>>> [ IO READ: use_cgtos ] <<<<< ..'
      call ezfio_get_ao_basis_use_cgtos(use_cgtos)
    else
      call ezfio_set_ao_basis_use_cgtos(use_cgtos)
    endif
  endif
  IRP_IF MPI_DEBUG
    print *,  irp_here, mpi_rank
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
  IRP_ENDIF
  IRP_IF MPI
    include 'mpif.h'
    integer :: ierr
    call MPI_BCAST( use_cgtos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
      stop 'Unable to read use_cgtos with MPI'
    endif
  IRP_ENDIF

!  call write_time(6)

END_PROVIDER

! ---

 BEGIN_PROVIDER [complex*16, ao_expo_cgtos_ord_transp, (ao_prim_num_max, ao_num)]
&BEGIN_PROVIDER [double precision, ao_expo_phase_ord_transp, (4, ao_prim_num_max, ao_num)]

  implicit none

  integer :: i, j, m

  do j = 1, ao_num
    do i = 1, ao_prim_num_max

      ao_expo_cgtos_ord_transp(i,j) = ao_expo_cgtos_ord(j,i)

      do m = 1, 4
        ao_expo_phase_ord_transp(m,i,j) = ao_expo_phase_ord(m,j,i)
      enddo
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_coef_norm_cgtos_ord, (ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [complex*16      , ao_expo_cgtos_ord, (ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [double precision, ao_expo_phase_ord, (4, ao_num, ao_prim_num_max)]

  implicit none

  integer          :: i, j, m
  integer          :: iorder(ao_prim_num_max)
  double precision :: d(ao_prim_num_max,7)

  d = 0.d0

  do i = 1, ao_num

    do j = 1, ao_prim_num(i)
      iorder(j) = j

      d(j,1) = ao_expo(i,j)
      d(j,2) = ao_coef_norm_cgtos(i,j)
      d(j,3) = ao_expo_im(i,j)

      do m = 1, 3
        d(j,3+m) = ao_expo_phase(m,i,j)
      enddo
      d(j,7) = d(j,4) + d(j,5) + d(j,6)
    enddo

    call dsort(d(1,1), iorder, ao_prim_num(i))
    do j = 2, 7
      call dset_order(d(1,j), iorder, ao_prim_num(i))
    enddo

    do j = 1, ao_prim_num(i)

      ao_expo_cgtos_ord     (i,j) = d(j,1) + (0.d0, 1.d0) * d(j,3)
      ao_coef_norm_cgtos_ord(i,j) = d(j,2)

      do m = 1, 4
        ao_expo_phase_ord(m,i,j) = d(j,3+m)
      enddo
    enddo
  enddo

END_PROVIDER

! ---

BEGIN_PROVIDER [double precision, ao_coef_cgtos_norm_ord_transp, (ao_prim_num_max, ao_num)]

  implicit none

  integer :: i, j

  do j = 1, ao_num
    do i = 1, ao_prim_num_max
      ao_coef_cgtos_norm_ord_transp(i,j) = ao_coef_norm_cgtos_ord(j,i)
    enddo
  enddo

END_PROVIDER

! ---

 BEGIN_PROVIDER [double precision, ao_coef_norm_cgtos, (ao_num, ao_prim_num_max)]
&BEGIN_PROVIDER [double precision, ao_coef_norm_factor_cgtos, (ao_num)]

  implicit none

  integer          :: i, j, k, ii, m, powA(3), nz
  double precision :: norm, c
  double precision :: phiA1, phiA2
  complex*16       :: C_A(3), expo1, expo2
  complex*16       :: overlap_x, overlap_y, overlap_z
  complex*16       :: integ1, integ2, C1, C2

  double precision, allocatable :: self_overlap(:)


  nz = 100

  C_A(1) = (0.d0, 0.d0)
  C_A(2) = (0.d0, 0.d0)
  C_A(3) = (0.d0, 0.d0)

  ao_coef_norm_cgtos = 0.d0

  ! GAMESS convention for primitive factors
  if(primitives_normalized) then

    do i = 1, ao_num

      ii = ao_nucl(i)
      powA(1) = ao_power(i,1)
      powA(2) = ao_power(i,2)
      powA(3) = ao_power(i,3)

      do j = 1, ao_prim_num(i)

        expo1 = ao_expo(i,j) + (0.d0, 1.d0) * ao_expo_im(i,j)
        phiA1 = ao_expo_phase(1,i,j) + ao_expo_phase(2,i,j) + ao_expo_phase(3,i,j)

        C1 = zexp(-(0.d0, 2.d0) * phiA1)
        C2 = (1.d0, 0.d0)

        call overlap_cgaussian_xyz(C_A, C_A, expo1, expo1, powA, powA, &
                                   C_A, C_A, overlap_x, overlap_y, overlap_z, integ1, nz)

        call overlap_cgaussian_xyz(conjg(C_A), C_A, conjg(expo1), expo1, powA, powA, &
                                   conjg(C_A), C_A, overlap_x, overlap_y, overlap_z, integ2, nz)

        norm = 0.5d0 * real(C1 * integ1 + C2 * integ2)

        ao_coef_norm_cgtos(i,j) = ao_coef(i,j) / dsqrt(norm)
      enddo
    enddo

  else

    do i = 1, ao_num
      do j = 1, ao_prim_num(i)
        ao_coef_norm_cgtos(i,j) = ao_coef(i,j)
      enddo
    enddo

  endif ! primitives_normalized



  if (ao_normalized) then

    allocate(self_overlap(ao_num))

    do i = 1, ao_num
      ii = ao_nucl(i)
      powA(1) = ao_power(i,1)
      powA(2) = ao_power(i,2)
      powA(3) = ao_power(i,3)
  
      self_overlap(i) = 0.d0
      do j = 1, ao_prim_num(i)

        expo1 = ao_expo(i,j) + (0.d0, 1.d0) * ao_expo_im(i,j)
        phiA1 = ao_expo_phase(1,i,j) + ao_expo_phase(2,i,j) + ao_expo_phase(3,i,j)

        do k = 1, j-1

          expo2 = ao_expo(i,k) + (0.d0, 1.d0) * ao_expo_im(i,k)
          phiA2 = ao_expo_phase(1,i,k) + ao_expo_phase(2,i,k) + ao_expo_phase(3,i,k)

          C1 = zexp((0.d0, 1.d0) * (-phiA1 - phiA2))
          C2 = zexp((0.d0, 1.d0) * ( phiA1 - phiA2))

          call overlap_cgaussian_xyz(C_A, C_A, expo1, expo2, powA, powA, &
                                     C_A, C_A, overlap_x, overlap_y, overlap_z, integ1, nz)
          call overlap_cgaussian_xyz(conjg(C_A), C_A, conjg(expo1), expo2, powA, powA, &
                                     conjg(C_A), C_A, overlap_x, overlap_y, overlap_z, integ2, nz)
          c = 0.5d0 * real(C1 * integ1 + C2 * integ2)

          self_overlap(i) = self_overlap(i) + 2.d0 * c * ao_coef_norm_cgtos(i,j) * ao_coef_norm_cgtos(i,k)
        enddo

        C1 = zexp(-(0.d0, 2.d0) * phiA1)
        C2 = (1.d0, 0.d0)

        call overlap_cgaussian_xyz(C_A, C_A, expo1, expo1, powA, powA, &
                                   C_A, C_A, overlap_x, overlap_y, overlap_z, integ1, nz)
        call overlap_cgaussian_xyz(conjg(C_A), C_A, conjg(expo1), expo1, powA, powA, &
                                   conjg(C_A), C_A, overlap_x, overlap_y, overlap_z, integ2, nz)
        c = 0.5d0 * real(C1 * integ1 + C2 * integ2)

        self_overlap(i) = self_overlap(i) + c * ao_coef_norm_cgtos(i,j) * ao_coef_norm_cgtos(i,k)
      enddo
    enddo

    do i = 1, ao_num
      ao_coef_norm_factor_cgtos(i) = 1.d0 / dsqrt(self_overlap(i))
    enddo

    deallocate(self_overlap)

  else

    do i = 1, ao_num
      ao_coef_norm_factor_cgtos(i) = 1.d0
    enddo
  endif

  do i = 1, ao_num
    do j = 1, ao_prim_num(i)
      ao_coef_norm_cgtos(i,j) = ao_coef_norm_cgtos(i,j) * ao_coef_norm_factor_cgtos(i)
    enddo
  enddo

END_PROVIDER

! ---


