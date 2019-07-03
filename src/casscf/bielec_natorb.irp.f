 BEGIN_PROVIDER [real*8, bielec_PQxx_no, (mo_num, mo_num,n_core_inact_act_orb,n_core_inact_act_orb)]
  BEGIN_DOC
  ! integral (pq|xx) in the basis of natural MOs
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q,pp
  double precision, allocatable  :: f(:,:,:), d(:,:,:)


  bielec_PQxx_no(:,:,:,:) = bielec_PQxx(:,:,:,:)

  allocate (f(n_act_orb,mo_num,n_core_inact_act_orb), &
      d(n_act_orb,mo_num,n_core_inact_act_orb))

  do l=1,n_core_inact_act_orb

    do k=1,n_core_inact_act_orb
      do j=1,mo_num
        do p=1,n_act_orb
          f(p,j,k)=bielec_PQxx_no(list_act(p),j,k,l)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,mo_num*n_core_inact_act_orb,n_act_orb,1.d0,  &
          natorbsCI, size(natorbsCI,1),                              &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do k=1,n_core_inact_act_orb
      do j=1,mo_num
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          bielec_PQxx_no(list_act(p),j,k,l)=d(pp,j,k)
        end do
      end do

      do j=1,mo_num
        do p=1,n_act_orb
          f(p,j,k)=bielec_PQxx_no(j,list_act(p),k,l)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,mo_num*n_core_inact_act_orb,n_act_orb,1.d0,  &
          natorbsCI, n_act_orb,                                      &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do k=1,n_core_inact_act_orb
      do p=1,n_act_orb
        pp=n_act_orb-p+1
        do j=1,mo_num
          bielec_PQxx_no(j,list_act(p),k,l)=d(pp,j,k)
        end do
      end do
    end do
  end do

  deallocate (f,d)

  allocate (f(mo_num,mo_num,n_act_orb),d(mo_num,mo_num,n_act_orb))

  do l=1,n_core_inact_act_orb

    do p=1,n_act_orb
      do k=1,mo_num
        do j=1,mo_num
          f(j,k,p) = bielec_PQxx_no(j,k,n_core_inact_orb+p,l)
        end do
      end do
    end do
    call dgemm('N','N',mo_num*mo_num,n_act_orb,n_act_orb,1.d0,       &
          f, mo_num*mo_num,                                          &
          natorbsCI, n_act_orb,                                      &
          0.d0,                                                      &
          d, mo_num*mo_num)
    do p=1,n_act_orb
      pp=n_act_orb-p+1
      do k=1,mo_num
        do j=1,mo_num
          bielec_PQxx_no(j,k,n_core_inact_orb+p,l)=d(j,k,pp)
        end do
      end do
    end do
  end do

  do l=1,n_core_inact_act_orb
    do p=1,n_act_orb
      do k=1,mo_num
        do j=1,mo_num
          f(j,k,p) = bielec_PQxx_no(j,k,l,n_core_inact_orb+p)
        end do
      end do
    end do
    call dgemm('N','N',mo_num*mo_num,n_act_orb,n_act_orb,1.d0,       &
          f, mo_num*mo_num,                                          &
          natorbsCI, n_act_orb,                                      &
          0.d0,                                                      &
          d, mo_num*mo_num)
    do p=1,n_act_orb
      pp=n_act_orb-p+1
      do k=1,mo_num
        do j=1,mo_num
          bielec_PQxx_no(j,k,l,n_core_inact_orb+p)=d(j,k,pp)
        end do
      end do
    end do
  end do

  deallocate (f,d)

END_PROVIDER



BEGIN_PROVIDER [real*8, bielec_PxxQ_no, (mo_num,n_core_inact_act_orb,n_core_inact_act_orb, mo_num)]
  BEGIN_DOC
  ! integral (px|xq) in the basis of natural MOs
  ! indices are unshifted orbital numbers
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q,pp
  double precision, allocatable  :: f(:,:,:), d(:,:,:)

  bielec_PxxQ_no(:,:,:,:) = bielec_PxxQ(:,:,:,:)

  allocate (f(n_act_orb,n_core_inact_act_orb,n_core_inact_act_orb), &
      d(n_act_orb,n_core_inact_act_orb,n_core_inact_act_orb))

  do j=1,mo_num
    do l=1,n_core_inact_act_orb
      do k=1,n_core_inact_act_orb
        do p=1,n_act_orb
            f(p,k,l) = bielec_PxxQ_no(list_act(p),k,l,j)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,n_core_inact_act_orb**2,n_act_orb,1.d0,  &
          natorbsCI, size(natorbsCI,1),                              &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do l=1,n_core_inact_act_orb
      do k=1,n_core_inact_act_orb
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          bielec_PxxQ_no(list_act(p),k,l,j)=d(pp,k,l)
        end do
      end do
    end do
  end do

  deallocate (f,d)

  allocate (f(n_act_orb,mo_num,n_core_inact_act_orb), &
    d(n_act_orb,mo_num,n_core_inact_act_orb))

  do k=1,mo_num
    do l=1,n_core_inact_act_orb
      do j=1,mo_num
        do p=1,n_act_orb
          f(p,j,l) = bielec_PxxQ_no(j,n_core_inact_orb+p,l,k)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,mo_num*n_core_inact_act_orb,n_act_orb,1.d0,  &
          natorbsCI, size(natorbsCI,1),                              &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do l=1,n_core_inact_act_orb
      do j=1,mo_num
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          bielec_PxxQ_no(j,n_core_inact_orb+p,l,k)=d(pp,j,l)
        end do
      end do
    end do
  end do

  deallocate(f,d)

  allocate(f(mo_num,n_core_inact_act_orb,n_act_orb), &
    d(mo_num,n_core_inact_act_orb,n_act_orb) )

  do k=1,mo_num
    do p=1,n_act_orb
      do l=1,n_core_inact_act_orb
        do j=1,mo_num
          f(j,l,p) = bielec_PxxQ_no(j,l,n_core_inact_orb+p,k)
        end do
      end do
    end do
    call dgemm('N','N',mo_num*n_core_inact_act_orb,n_act_orb,n_act_orb,1.d0,  &
          f, mo_num*n_core_inact_act_orb,                                       &
          natorbsCI, size(natorbsCI,1),                              &
          0.d0,                                                      &
          d, mo_num*n_core_inact_act_orb)
    do p=1,n_act_orb
      pp=n_act_orb-p+1
      do l=1,n_core_inact_act_orb
        do j=1,mo_num
          bielec_PxxQ_no(j,l,n_core_inact_orb+p,k)=d(j,l,pp)
        end do
      end do
    end do
  end do


  do l=1,n_core_inact_act_orb
    do p=1,n_act_orb
      do k=1,n_core_inact_act_orb
        do j=1,mo_num
          f(j,k,p) = bielec_PxxQ_no(j,k,l,n_core_inact_orb+p)
        end do
      end do
    end do
    call dgemm('N','N',mo_num*n_core_inact_act_orb,n_act_orb,n_act_orb,1.d0,  &
          f, mo_num*n_core_inact_act_orb,                                       &
          natorbsCI, size(natorbsCI,1),                              &
          0.d0,                                                      &
          d, mo_num*n_core_inact_act_orb)
    do p=1,n_act_orb
      pp=n_act_orb-p+1
      do k=1,n_core_inact_act_orb
        do j=1,mo_num
          bielec_PxxQ_no(j,k,l,n_core_inact_orb+p)=d(j,k,pp)
        end do
      end do
    end do
  end do

END_PROVIDER


BEGIN_PROVIDER [real*8, bielecCI_no, (n_act_orb,n_act_orb,n_act_orb, mo_num)]
  BEGIN_DOC
  ! integrals (tu|vp) in the basis of natural MOs
  ! index p runs over the whole basis, t,u,v only over the active orbitals
  END_DOC
  implicit none
  integer                        :: i,j,k,l,t,u,p,q,pp
  double precision, allocatable  :: f(:,:,:), d(:,:,:)
  
  bielecCI_no(:,:,:,:) = bielecCI(:,:,:,:)

  allocate (f(n_act_orb,n_act_orb,mo_num), &
      d(n_act_orb,n_act_orb,mo_num))

  do l=1,mo_num
    do k=1,n_act_orb
      do j=1,n_act_orb
        do p=1,n_act_orb
          f(p,j,k)=bielecCI_no(p,j,k,l)
        end do
      end do
    end do
    call dgemm('T','N',n_act_orb,n_act_orb*n_act_orb,n_act_orb,1.d0,  &
          natorbsCI, size(natorbsCI,1),                              &
          f, n_act_orb,                                              &
          0.d0,                                                      &
          d, n_act_orb)
    do k=1,n_act_orb
      do j=1,n_act_orb
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          bielecCI_no(p,j,k,l)=d(pp,j,k)
        end do
      end do
    end do
  end do

  do j=1,n_act_orb
    do k=1,n_act_orb
      do l=1,mo_num
        do p=1,n_act_orb
          d(p,1,1)=0.D0
        end do
        do p=1,n_act_orb
          do q=1,n_act_orb
            d(p,1,1)+=bielecCI_no(j,q,k,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          bielecCI_no(j,p,k,l)=d(pp,1,1)
        end do
      end do
    end do
  end do
  ! 3rd quarter
  do j=1,n_act_orb
    do k=1,n_act_orb
      do l=1,mo_num
        do p=1,n_act_orb
          d(p,1,1)=0.D0
        end do
        do p=1,n_act_orb
          do q=1,n_act_orb
            d(p,1,1)+=bielecCI_no(j,k,q,l)*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          bielecCI_no(j,k,p,l)=d(pp,1,1)
        end do
      end do 
    end do  
  end do  

!  do l=1,mo_num
!    do k=1,n_act_orb
!      do p=1,n_act_orb
!        do j=1,n_act_orb
!          f(j,p,k)=bielecCI_no(j,p,k,l)
!        end do
!      end do
!    end do
!    call dgemm('T','N',n_act_orb,n_act_orb*n_act_orb,n_act_orb,1.d0,  &
!          natorbsCI, n_act_orb,                                      &
!          f, n_act_orb,                                              &
!          0.d0,                                                      &
!          d, n_act_orb)
!    do k=1,n_act_orb
!      do p=1,n_act_orb
!        pp=n_act_orb-p+1
!        do j=1,n_act_orb
!          bielecCI_no(j,p,k,l)=d(j,pp,k)
!        end do
!      end do
!    end do
!
!    do p=1,n_act_orb
!      do k=1,n_act_orb
!        do j=1,n_act_orb
!          f(j,k,p)=bielecCI_no(j,k,p,l)
!        end do
!      end do
!    end do
!    call dgemm('N','N',n_act_orb*n_act_orb,n_act_orb,n_act_orb,1.d0,       &
!          f, n_act_orb*n_act_orb,                                          &
!          natorbsCI, n_act_orb,                                      &
!          0.d0,                                                      &
!          d, n_act_orb*n_act_orb)
!
!    do p=1,n_act_orb
!      pp=n_act_orb-p+1
!      do k=1,n_act_orb
!        do j=1,n_act_orb
!          bielecCI_no(j,k,p,l)=d(j,k,pp)
!        end do
!      end do 
!    end do  
!  end do  
!
  do l=1,n_act_orb
    do k=1,n_act_orb
      do j=1,n_act_orb
        do p=1,n_act_orb
          d(p,1,1)=0.D0
        end do
        do p=1,n_act_orb
          do q=1,n_act_orb
            d(p,1,1)+=bielecCI_no(j,k,l,list_act(q))*natorbsCI(q,p)
          end do
        end do
        do p=1,n_act_orb
          pp=n_act_orb-p+1
          bielecCI_no(j,k,l,list_act(p))=d(pp,1,1)
        end do
      end do 
    end do  
  end do  

END_PROVIDER

