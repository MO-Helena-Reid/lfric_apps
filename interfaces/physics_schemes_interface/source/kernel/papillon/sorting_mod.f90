! Created by  on 17/11/2023.

module sorting_mod
use ml_constants_mod, only: r_ml
implicit none
contains
subroutine insertionsort(array)
    implicit none
    real(kind=r_ml), dimension(:), intent(inout) :: array
    integer :: i,j
    real(kind=r_ml) :: temp

    do i=2,size(array)
        temp=array(i)
        do j=i-1,1,-1
          if (array(j)<=temp) exit
          array(j+1)=array(j)
        enddo
        array(j+1)=temp
    enddo
end subroutine insertionsort
subroutine insertionsort_withindices(array, indices)
    implicit none
    real(kind=r_ml), dimension(:), intent(inout) :: array
    integer, dimension(:), intent(out) :: indices
    integer :: i,j,itemp
    real(kind=r_ml) :: temp
    if (size(array)/=size(indices)) then
        print*, 'ERROR: size of array does not match size of indices:', size(array), size(indices)
    end if
!   initialise indices
    do i=1, size(array)
        indices(i) = i
    end do

    do i=2,size(array)
        temp=array(i)
        itemp=indices(i)
        do j=i-1,1,-1
            if (array(j)<=temp) exit
            array(j+1)=array(j)
            indices(j+1)=indices(j)
        enddo
        array(j+1)=temp
        indices(j+1)=itemp
    enddo
end subroutine insertionsort_withindices
end module sorting_mod