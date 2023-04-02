module regcoil_init_Fourier_modes_mod

  implicit none

contains

  subroutine regcoil_init_Fourier_modes(mpol, ntor, mnmax, xm, xn, include_00,helicity_ratio)

    implicit none

    integer :: mpol, ntor, mnmax, helicity_ratio
    integer, dimension(:), allocatable :: xm, xn
    logical, intent(in) :: include_00
    
    integer :: jn, jm, index, iflag
    integer, dimension(:), allocatable :: xm_temp, xn_temp
    
    ! xm is nonnegative.
    ! xn can be negative, zero, or positive.
    ! When xm is 0, xn must be 0 or positive.
    mnmax = mpol*(ntor*2+1) + ntor! make this just the number of basis modes I have...
    if (helicity_ratio .ne. 0) mnmax = ntor ! we only will have ntor modes, only ntor matters

    if (include_00) mnmax = mnmax + 1! ignore for now

   ! helicity_ratio is s.t. m = - (helicity_ratio) * n

    if (allocated(xm)) deallocate(xm)
    allocate(xm(mnmax),stat=iflag)
    if (iflag .ne. 0) stop 'init_Fourier Allocation error!'
 
    if (allocated(xn)) deallocate(xn)
    allocate(xn(mnmax),stat=iflag)
    if (iflag .ne. 0) stop 'init_Fourier Allocation error!'

! how this will work is we will have 
! xn go from 1:ntor, then xm = -c * xn
! this way we ensure that every fourier basis mode
! is a fxn of the form cos(xn*(-c*u - v)) i.e. they all have the same helicity xm/xn=-c

    xm = 0
    xn = 0

    ! Handle the xm=0 modes:
    index = 0
    if (include_00) index = 1
    if (helicity_ratio .eq. 0) then
       do jn = 1, ntor
          index = index + 1
          xn(index)=jn
       end do
    end if
    ! Handle the xm>0 modes:
    
    if (helicity_ratio .eq. 0) then
       do jm = 1,mpol
          do jn = -ntor, ntor
             index = index + 1
             xn(index) = jn
             xm(index) = jm
          end do
       end do
    end if
    !do jm = 1,mpol
! if helicity ratio is not zero
    if (helicity_ratio .ne. 0) then
       do jm = 1, mpol
          index = index + 1
          xm(index) = jm
          xn(index) = -helicity_ratio * jm
       end do
    end if
    !end do
    
    if (index .ne. mnmax) then
       print *,"Error!  index=",index," but mnmax=",mnmax
       stop
    end if

  end subroutine regcoil_init_Fourier_modes

subroutine regcoil_init_Fourier_modes_sensitivity &
      (mpol,ntor,mnmax,nomega,xm,xn,omega,sensitivity_symmetry_option)

    implicit none

    integer :: mpol, ntor, mnmax, iomega, i, nomega
    integer, dimension(:), allocatable :: xm, xn, omega
    integer :: minSymmetry, maxSymmetry, sensitivity_symmetry_option

    integer :: jn, jm, iflag

    ! xm is nonnegative.
    ! xn can be negative, zero, or positive.
    ! When xm is 0, xn must be non-negative
    mnmax = mpol*(ntor*2+1) + ntor+1

    ! nomega is the length of the number of fourier coefficients
    select case (sensitivity_symmetry_option)
      case (1) ! stellarator symmetric
        nomega = mnmax*2
        minSymmetry = 1
        maxSymmetry = 2
      case (2) ! even in theta and zeta
        nomega = mnmax*2
        minSymmetry = 3
        maxSymmetry = 4
      case (3) ! no symmetry
        nomega = mnmax*4
        minSymmetry = 1
        maxSymmetry = 4
    end select

    allocate(xm(nomega),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(xn(nomega),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'
    allocate(omega(nomega),stat=iflag)
    if (iflag .ne. 0) stop 'Allocation error!'

    ! Handle the xm=0 modes:
    xm=0
    iomega = 0
    do jn=1,ntor+1
      do i=minSymmetry,maxSymmetry
        iomega = iomega + 1
        omega(iomega) = i
        xn(iomega)=jn-1
      enddo
    end do

    ! Handle the xm>0 modes:
    do jm = 1,mpol
      do jn = -ntor, ntor
        do i=minSymmetry,maxSymmetry
          iomega = iomega + 1
          xn(iomega) = jn
          xm(iomega) = jm
          omega(iomega) = i
        enddo
      end do
    end do

    if (iomega .ne. nomega) then
      print *,"Error!  iomega=",iomega," but nomega=",nomega
      stop
    end if

  end subroutine regcoil_init_Fourier_modes_sensitivity
  
end module regcoil_init_Fourier_modes_mod
