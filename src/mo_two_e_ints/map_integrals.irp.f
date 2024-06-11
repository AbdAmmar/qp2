use map_module

!! MO Map
!! ======

BEGIN_PROVIDER [ type(map_type), mo_integrals_map ]
  implicit none
  BEGIN_DOC
  ! MO integrals
  END_DOC
  integer(key_kind)              :: key_max
  integer(map_size_kind)         :: sze
  call two_e_integrals_index(mo_num,mo_num,mo_num,mo_num,key_max)
  sze = key_max
  call map_init(mo_integrals_map,sze)
  print*, 'MO map initialized: ', sze
END_PROVIDER

subroutine insert_into_mo_integrals_map(n_integrals,                 &
      buffer_i, buffer_values, thr)
  use map_module
  implicit none

  BEGIN_DOC
  ! Create new entry into MO map, or accumulate in an existing entry
  END_DOC

  integer, intent(in)                :: n_integrals
  integer(key_kind), intent(inout)   :: buffer_i(n_integrals)
  real(integral_kind), intent(inout) :: buffer_values(n_integrals)
  real(integral_kind), intent(in)    :: thr
  call map_update(mo_integrals_map, buffer_i, buffer_values, n_integrals, thr)
end

 BEGIN_PROVIDER [ integer*4, mo_integrals_cache_min ]
&BEGIN_PROVIDER [ integer*4, mo_integrals_cache_max ]
&BEGIN_PROVIDER [ integer*4, mo_integrals_cache_shift]
&BEGIN_PROVIDER [ integer*4, mo_integrals_cache_size ]
 implicit none
 BEGIN_DOC
 ! Min and max values of the MOs for which the integrals are in the cache
 END_DOC
 mo_integrals_cache_shift = 7  ! 7 = log(128). Max 7
 mo_integrals_cache_size  = 2**mo_integrals_cache_shift

 mo_integrals_cache_min = max(1,elec_alpha_num - (mo_integrals_cache_size/2 - 1) )
 mo_integrals_cache_max = min(mo_num, mo_integrals_cache_min + mo_integrals_cache_size - 1)

END_PROVIDER

BEGIN_PROVIDER [ double precision, mo_integrals_cache, (0:mo_integrals_cache_size**4) ]
 implicit none
 BEGIN_DOC
 ! Cache of MO integrals for fast access
 END_DOC
 PROVIDE mo_two_e_integrals_in_map
 integer                        :: i,j,k,l
 integer                        :: ii
 integer(key_kind)              :: idx
 real(integral_kind)            :: integral
 FREE ao_integrals_cache
 if (do_mo_cholesky) then

   call set_multiple_levels_omp(.False.)
   !$OMP PARALLEL DO PRIVATE (k,l,ii)
   do l=mo_integrals_cache_min,mo_integrals_cache_max
     do k=mo_integrals_cache_min,mo_integrals_cache_max
         ii = l-mo_integrals_cache_min
         ii = ior( shiftl(ii,mo_integrals_cache_shift), k-mo_integrals_cache_min)
         ii = shiftl(ii,mo_integrals_cache_shift)
         ii = shiftl(ii,mo_integrals_cache_shift)
         call dgemm('T','N', mo_integrals_cache_max-mo_integrals_cache_min+1, &
                             mo_integrals_cache_max-mo_integrals_cache_min+1, &
           cholesky_mo_num, 1.d0, &
           cholesky_mo_transp(1,mo_integrals_cache_min,k), cholesky_mo_num, &
           cholesky_mo_transp(1,mo_integrals_cache_min,l), cholesky_mo_num, 0.d0, &
           mo_integrals_cache(ii), mo_integrals_cache_size)
     enddo
   enddo
   !$OMP END PARALLEL DO

 else
   !$OMP PARALLEL DO PRIVATE (i,j,k,l,idx,ii,integral)
   do l=mo_integrals_cache_min,mo_integrals_cache_max
     do k=mo_integrals_cache_min,mo_integrals_cache_max
       do j=mo_integrals_cache_min,mo_integrals_cache_max
         do i=mo_integrals_cache_min,mo_integrals_cache_max
           !DIR$ FORCEINLINE
           call two_e_integrals_index(i,j,k,l,idx)
           !DIR$ FORCEINLINE
           call map_get(mo_integrals_map,idx,integral)
           ii = l-mo_integrals_cache_min
           ii = ior( shiftl(ii,mo_integrals_cache_shift), k-mo_integrals_cache_min)
           ii = ior( shiftl(ii,mo_integrals_cache_shift), j-mo_integrals_cache_min)
           ii = ior( shiftl(ii,mo_integrals_cache_shift), i-mo_integrals_cache_min)
           mo_integrals_cache(ii) = integral
         enddo
       enddo
     enddo
   enddo
   !$OMP END PARALLEL DO
 endif

END_PROVIDER


double precision function get_two_e_integral(i,j,k,l,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  integer(key_kind)              :: idx
  integer                        :: ii
  type(map_type), intent(inout)  :: map
  real(integral_kind)            :: tmp
  integer                        :: kk

  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache do_mo_cholesky

  if (use_banned_excitation) then
    if (banned_excitation(i,k)) then
      get_two_e_integral = 0.d0
      return
    endif
    if (banned_excitation(j,l)) then
      get_two_e_integral = 0.d0
      return
    endif
  endif


  ii = l-mo_integrals_cache_min
  ii = ior(ii, k-mo_integrals_cache_min)
  ii = ior(ii, j-mo_integrals_cache_min)
  ii = ior(ii, i-mo_integrals_cache_min)

  if (iand(ii, -mo_integrals_cache_size) /= 0) then
    ! Integral is not in the cache

    if  (do_mo_cholesky) then

      double precision, external :: ddot
      get_two_e_integral = ddot(cholesky_mo_num, cholesky_mo_transp(1,i,k), 1, cholesky_mo_transp(1,j,l), 1)

    else
      ! Integrals is in the map

      !DIR$ FORCEINLINE
      call two_e_integrals_index(i,j,k,l,idx)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
      get_two_e_integral = dble(tmp)
    endif

  else
    ! Integrals is in the cache

    ii = l-mo_integrals_cache_min
    ii = ior( shiftl(ii,mo_integrals_cache_shift), k-mo_integrals_cache_min)
    ii = ior( shiftl(ii,mo_integrals_cache_shift), j-mo_integrals_cache_min)
    ii = ior( shiftl(ii,mo_integrals_cache_shift), i-mo_integrals_cache_min)
    get_two_e_integral = mo_integrals_cache(ii)

  endif
end


double precision function mo_two_e_integral(i,j,k,l)
  implicit none
  BEGIN_DOC
  ! Returns one integral <ij|kl> in the MO basis
  END_DOC
  integer, intent(in)            :: i,j,k,l
  double precision               :: get_two_e_integral
  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache
  !DIR$ FORCEINLINE
  mo_two_e_integral = get_two_e_integral(i,j,k,l,mo_integrals_map)
  return
end

subroutine get_mo_two_e_integrals(j,k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i for j,k,l fixed.
  END_DOC
  integer, intent(in)            :: j,k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  double precision, external :: get_two_e_integral

  integer                        :: ii
  real(integral_kind)            :: tmp
  integer(key_kind)              :: i1, idx
  integer(key_kind)              :: p,q,r,s,i2
  PROVIDE mo_two_e_integrals_in_map mo_integrals_cache


  out_val(1:sze) = 0.d0
  if (banned_excitation(j,l)) then
    return
  endif

  ii0 = l-mo_integrals_cache_min
  ii0 = ior(ii0, k-mo_integrals_cache_min)
  ii0 = ior(ii0, j-mo_integrals_cache_min)

  integer :: ii0, ii0_8, ii_8
  ii0_8 = l-mo_integrals_cache_min
  ii0_8 = ior( shiftl(ii0_8,mo_integrals_cache_shift), k-mo_integrals_cache_min)
  ii0_8 = ior( shiftl(ii0_8,mo_integrals_cache_shift), j-mo_integrals_cache_min)

  q = min(j,l)
  s = max(j,l)
  q = q+shiftr(s*s-s,1)

  do i=1,sze
    if (banned_excitation(i,k)) cycle
    ii = ior(ii0, i-mo_integrals_cache_min)
    if (iand(ii, -mo_integrals_cache_size) == 0) then
      ii_8 = ior( shiftl(ii0_8,mo_integrals_cache_shift), i-mo_integrals_cache_min)
      out_val(i) = mo_integrals_cache(ii_8)
    else
      p = min(i,k)
      r = max(i,k)
      p = p+shiftr(r*r-r,1)
      i1 = min(p,q)
      i2 = max(p,q)
      idx = i1+shiftr(i2*i2-i2,1)
      !DIR$ FORCEINLINE
      call map_get(map,idx,tmp)
      out_val(i) = dble(tmp)
    endif
  enddo

!  if (banned_excitation(j,l)) then
!      out_val(1:sze) = 0.d0
!      return
!  endif
!
!  if (mo_integrals_cache_min > 1) then
!
!    if (do_mo_cholesky) then
!
!      call dgemv('T', cholesky_mo_num, mo_integrals_cache_min-1, 1.d0, &
!         cholesky_mo_transp(1,1,k), cholesky_mo_num, &
!         cholesky_mo_transp(1,j,l), 1, 0.d0, &
!         out_val, 1)
!
!    else
!
!      q = min(j,l)
!      s = max(j,l)
!      q = q+shiftr(s*s-s,1)
!
!      do i=1,mo_integrals_cache_min-1
!        if (banned_excitation(i,k)) then
!          out_val(i) = 0.d0
!          cycle
!        endif
!        p = min(i,k)
!        r = max(i,k)
!        p = p+shiftr(r*r-r,1)
!        i1 = min(p,q)
!        i2 = max(p,q)
!        idx = i1+shiftr(i2*i2-i2,1)
!        !DIR$ FORCEINLINE
!        call map_get(map,idx,tmp)
!        out_val(i) = dble(tmp)
!      enddo
!
!    endif
!
!  endif
!
!
!  ii = l-mo_integrals_cache_min
!  ii = ior( shiftl(ii, mo_integrals_cache_shift), k-mo_integrals_cache_min)
!  ii = ior( shiftl(ii, mo_integrals_cache_shift), j-mo_integrals_cache_min)
!  ii = shiftl(ii, mo_integrals_cache_shift)
!  do i=mo_integrals_cache_min, mo_integrals_cache_max
!    ii = ii+1
!    out_val(i) = mo_integrals_cache(ii)
!  enddo
!
!
!  if (mo_integrals_cache_max < mo_num) then
!
!    if (do_mo_cholesky) then
!
!      call dgemv('T', cholesky_mo_num, mo_num-mo_integrals_cache_max, 1.d0, &
!         cholesky_mo_transp(1,mo_integrals_cache_max+1,k), cholesky_mo_num, &
!         cholesky_mo_transp(1,j,l), 1, 0.d0, &
!         out_val(mo_integrals_cache_max+1), 1)
!
!    else
!
!      q = min(j,l)
!      s = max(j,l)
!      q = q+shiftr(s*s-s,1)
!
!      do i=mo_integrals_cache_max+1,mo_num
!        if (banned_excitation(i,k)) then
!          out_val(i) = 0.d0
!          cycle
!        endif
!        p = min(i,k)
!        r = max(i,k)
!        p = p+shiftr(r*r-r,1)
!        i1 = min(p,q)
!        i2 = max(p,q)
!        idx = i1+shiftr(i2*i2-i2,1)
!        !DIR$ FORCEINLINE
!        call map_get(map,idx,tmp)
!        out_val(i) = dble(tmp)
!      enddo
!
!    endif
!
!  endif

end

subroutine get_mo_two_e_integrals_ij(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ij|kl> in the MO basis, all
  ! i(1)j(2) 1/r12 k(1)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: j
  real(integral_kind), allocatable :: tmp_val(:)

  if (do_mo_cholesky) then
    call dgemm('T', 'N', mo_num, mo_num, cholesky_mo_num, 1.d0, &
       cholesky_mo_transp(1,1,k), cholesky_mo_num, &
       cholesky_mo_transp(1,1,l), cholesky_mo_num, 0.d0, &
       out_array, sze)
  else
    do j=1,sze
      call get_mo_two_e_integrals(j,k,l,sze,out_array(1,j),map)
    enddo
  endif
end

subroutine get_mo_two_e_integrals_i1j1(k,l,sze,out_array,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ik|jl> in the MO basis, all
  ! i(1)j(1) 1/r12 k(2)l(2)
  ! i, j for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_array(sze,sze)
  type(map_type), intent(inout)  :: map
  integer                        :: j
  PROVIDE mo_two_e_integrals_in_map

  if (do_mo_cholesky) then

    call dgemv('T', cholesky_mo_num, mo_num*mo_num, 1.d0, &
       cholesky_mo_transp(1,1,1), cholesky_mo_num, &
       cholesky_mo_transp(1,k,l), 1, 0.d0, &
       out_array, 1)

  else

    do j=1,sze
      call get_mo_two_e_integrals(k,j,l,sze,out_array(1,j),map)
    enddo

  endif

end


subroutine get_mo_two_e_integrals_coulomb_ii(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|li>
  ! k(1)i(2) 1/r12 l(1)i(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  double precision, external     :: get_two_e_integral
  PROVIDE mo_two_e_integrals_in_map

  if (do_mo_cholesky) then

    call dgemv('T', cholesky_mo_num, mo_num, 1.d0, &
       cholesky_mo_transp(1,1,1), cholesky_mo_num*(mo_num+1), &
       cholesky_mo_transp(1,k,l), 1, 0.d0, &
       out_val, 1)

  else

    do i=1,sze
      out_val(i) = get_two_e_integral(k,i,l,i,map)
    enddo

  endif

end

subroutine get_mo_two_e_integrals_exch_ii(k,l,sze,out_val,map)
  use map_module
  implicit none
  BEGIN_DOC
  ! Returns multiple integrals <ki|il>
  ! k(1)i(2) 1/r12 i(1)l(2) :: out_val(i1)
  ! for k,l fixed.
  END_DOC
  integer, intent(in)            :: k,l, sze
  double precision, intent(out)  :: out_val(sze)
  type(map_type), intent(inout)  :: map
  integer                        :: i
  double precision, external     :: get_two_e_integral
  PROVIDE mo_two_e_integrals_in_map

  do i=1,sze
    out_val(i) = get_two_e_integral(k,i,i,l,map)
  enddo

end

 BEGIN_PROVIDER [ logical, banned_excitation, (mo_num,mo_num) ]
&BEGIN_PROVIDER [ logical, use_banned_excitation  ]
 implicit none
 use map_module
 BEGIN_DOC
 ! If true, the excitation is banned in the selection. Useful with local MOs.
 END_DOC
 banned_excitation = .False.
 use_banned_excitation = .False.

 integer :: i,j, icount
 integer(key_kind)              :: idx
 double precision :: tmp

!icount = 1 ! Avoid division by zero
!do j=1,mo_num
!  do i=1,j-1
!   call two_e_integrals_index(i,j,j,i,idx)
!   !DIR$ FORCEINLINE
!   call map_get(mo_integrals_map,idx,tmp)
!   banned_excitation(i,j) = dabs(tmp) < 1.d-14
!   banned_excitation(j,i) = banned_excitation(i,j)
!   if (banned_excitation(i,j)) icount = icount+2
! enddo
!enddo
!use_banned_excitation =  (mo_num*mo_num) / icount <= 100  !1%
!if (use_banned_excitation) then
!  print *, 'Using sparsity of exchange integrals'
!endif

END_PROVIDER



integer*8 function get_mo_map_size()
  implicit none
  BEGIN_DOC
  ! Return the number of elements in the MO map
  END_DOC
  get_mo_map_size = mo_integrals_map % n_elements
end


subroutine dump_mo_integrals(filename)
  use map_module
  implicit none
  BEGIN_DOC
  ! Save to disk the |MO| integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer*8                      :: i,j, n
  if (.not.mpi_master) then
    return
  endif
  call ezfio_set_work_empty(.False.)
  open(unit=66,file=filename,FORM='unformatted')
  write(66) integral_kind, key_kind
  write(66) mo_integrals_map%sorted, mo_integrals_map%map_size,    &
      mo_integrals_map%n_elements
  do i=0_8,mo_integrals_map%map_size
    write(66) mo_integrals_map%map(i)%sorted, mo_integrals_map%map(i)%map_size,&
        mo_integrals_map%map(i)%n_elements
  enddo
  do i=0_8,mo_integrals_map%map_size
    key => mo_integrals_map%map(i)%key
    val => mo_integrals_map%map(i)%value
    n = mo_integrals_map%map(i)%n_elements
    write(66) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  close(66)

end


integer function load_mo_integrals(filename)
  implicit none
  BEGIN_DOC
  ! Read from disk the |MO| integrals
  END_DOC
  character*(*), intent(in)      :: filename
  integer*8                      :: i
  integer(cache_key_kind), pointer :: key(:)
  real(integral_kind), pointer   :: val(:)
  integer                        :: iknd, kknd
  integer*8                      :: n, j
  load_mo_integrals = 1
  open(unit=66,file=filename,FORM='unformatted',STATUS='UNKNOWN')
  call lock_io()
  read(66,err=98,end=98) iknd, kknd
  if (iknd /= integral_kind) then
    print *,  'Wrong integrals kind in file :', iknd
    stop 1
  endif
  if (kknd /= key_kind) then
    print *,  'Wrong key kind in file :', kknd
    stop 1
  endif
  read(66,err=98,end=98) mo_integrals_map%sorted, mo_integrals_map%map_size,&
      mo_integrals_map%n_elements
  do i=0_8, mo_integrals_map%map_size
    read(66,err=99,end=99) mo_integrals_map%map(i)%sorted,          &
        mo_integrals_map%map(i)%map_size, mo_integrals_map%map(i)%n_elements
    call cache_map_reallocate(mo_integrals_map%map(i),mo_integrals_map%map(i)%map_size)
  enddo
  do i=0_8, mo_integrals_map%map_size
    key => mo_integrals_map%map(i)%key
    val => mo_integrals_map%map(i)%value
    n = mo_integrals_map%map(i)%n_elements
    read(66,err=99,end=99) (key(j), j=1,n), (val(j), j=1,n)
  enddo
  call unlock_io()
  call map_sort(mo_integrals_map)
  load_mo_integrals = 0
  return
  99 continue
  call map_deinit(mo_integrals_map)
  98 continue
  stop 'Problem reading mo_integrals_map file in work/'

end

