

! Скалярное произведение векторов
real function scalar_vector_multipilication(A, B, n)
integer, intent(in) :: n
real, intent(in) :: A(n), B(n)
integer :: i
scalar_vector_multipilication = 0
do i=1,n
    scalar_vector_multipilication = scalar_vector_multipilication + A(i) * B(i)
end do
end function scalar_vector_multipilication


program name
    implicit none
    real, allocatable :: A(:,:), B(:), X(:), F(:), rk(:), tmp(:), tmpResult(:), tmp1(:)
	integer :: n, i, j, iter, endVal
	real :: tau, SM, SMD, eps=0.01, tmpReal, scalar_vector_multipilication
	logical :: check
    interface 
        subroutine show_matrix(matrix, n)
            implicit none
            integer :: n,i,j
            real, intent(in) :: matrix(:,:)
        end subroutine show_matrix
    end interface
    interface 
        subroutine show_vector(vector, n)
            implicit none
            integer :: n,i
            real, intent(in) :: vector(:)
        end subroutine show_vector
    end interface
    interface
        subroutine multiplication_matrix_vector(matrix, vector, n, resultVector)
            implicit none
            integer :: n,i,j
            real, intent(in) :: vector(:), matrix(:,:)
            real, intent(in out) :: resultVector(:)
        end subroutine multiplication_matrix_vector
    end interface
    interface
        subroutine vectors_difference(A,B,n, result)
            implicit none
            integer :: n,i
            real, intent(in) :: A(:), B(:)
            real, intent(in out) :: result(:)
        end subroutine vectors_difference
    end interface
    interface
        subroutine calculate_f(A,X,B, n,result)
            implicit none
            integer :: n,i,j
            real, intent(in) :: A(:,:), B(:), X(:)
            real, intent(in out) :: result(:)
            real :: tmp(n)
        end subroutine calculate_f
    end interface
    interface
        subroutine check_result(vector,n, result)
            implicit none
            integer :: n,i,number,result
            real, intent(in) :: vector(:)
        end subroutine check_result
    end interface
    interface
        subroutine multiplication_vector_number(vector, number,n, result)
            implicit none
            integer :: n,i
            real :: number
            real, intent(in) :: vector(:)
            real, intent(in out) :: result(:)
        end subroutine multiplication_vector_number
    end interface


    write (*, *) "Начало программы"
	
	! Выделяем память
	open(1,file="matrix.f90")
	read(1,*) n

	allocate(A(n,n))
	allocate(B(n))
	allocate(X(n))
    allocate(F(n))
    allocate(rk(n))
    allocate(tmp(n))
    allocate(tmp1(n))
    allocate(tmpResult(n))
    tmpReal = 0

	! Считываем данные
	do i=1,n
		read(1,*)(A(i,j), j=1,n), B(i)
	end do


    open(1,file="firstApproach.f90")
	
    do i=1,n
        X(i) = B(i)/A(i,i)
        if (A(i,i) == 0) then
            X(i) = B(i)
        endif
        
    end do

	! Вывод матрицы
	write(*,*) "Размерность матрицы = ", n
	write(*,*) ""
	
    write(*,*) "Исходная матрица A"
	call show_matrix(A,n)

    write(*,*) "Исходный вектор B"
	call show_vector(B,n)

    write(*,*) "Первое приближение"
    call show_vector(X,n)

    
    call calculate_f(A,X,B,n,F)

    write(*,*) "F"
    endVal = 0
    iter = 0
    do while(endVal == 0)
        tmp = 0
        tmp1 = 0
        tmpReal = 0
        tmpResult = 0
        write(*,*) "Проверка условия"
        write(*,*) "Результат:"
        write(*,*) (endVal)
        write(*,*) "Рассчёт F"
        call show_vector(F,n)
        write(*,*) "rk = F*-1"
        call multiplication_vector_number(F, -1.0, n , rk)
        call show_vector(rk,n)
        write(*,*) "Перемножаем A на Rk"
        call multiplication_matrix_vector(A, rk,n, tmp)
        call show_vector(tmp,n)
        write(*,*) "Частное скалярных произведений"
        tmpReal = dot_product(rk,rk)/dot_product(rk,tmp)
        write(*,*) (tmpReal)
        write(*,*) "F * скалярное произведение"
        call multiplication_vector_number(F, tmpReal, n, tmp1)
        write(*,*) "tmp:"
        call show_vector(tmp1, n)
        call vectors_difference(X, tmp1, n, tmpResult)
        X = tmpResult
        write(*,*) "Результат шага X:"
        call show_vector(X,n)
        iter = iter + 1
        call calculate_f(A,X,B,n,F)
        call check_result(F,n,endVal)
    enddo


    write(*,*) "F"
    call show_vector(F,n)
end program name

!Объявляемые процедуры


!Вывод матрицы в терминал
subroutine show_matrix(matrix, n)
    implicit none
    integer :: n,i,j
    real, intent(in) :: matrix(:,:)

    do, i=1,n
        write(*,*)(matrix(i,j), j=1,n)
    enddo
end subroutine show_matrix

!Вывод вектора в терминал
subroutine show_vector(vector,n)
    implicit none
    integer :: n,i
    real, intent(in) :: vector(:)

    do, i=1,n
        write(*,*)(vector(i))
    enddo
end subroutine show_vector


!Перемножение матрицы на вектор
subroutine multiplication_matrix_vector(matrix, vector, n, resultVector)
    implicit none
    integer :: n,i,j
    real, intent(in) :: vector(:), matrix(:,:)
    real, intent(in out) :: resultVector(:)

    do, i=1,n
        resultVector(i) = 0
        do, j=1,n
            resultVector(i) = resultVector(i) + matrix(i,j) * vector(j)
        enddo
    end do

end subroutine multiplication_matrix_vector

!Вычисление результата Ax-b
subroutine calculate_f(A,X,B, n,result)
    implicit none
    integer :: n,i
    real, intent(in) :: A(:,:), B(:), X(:)
    real, intent(in out) :: result(:)
    real :: tmp(n)
    interface
        subroutine multiplication_matrix_vector(matrix, vector, n, resultVector)
            implicit none
            integer :: n,i,j
            real, intent(in) :: vector(:), matrix(:,:)
            real, intent(in out) :: resultVector(:)
        end subroutine multiplication_matrix_vector
    end interface
    interface
        subroutine vectors_difference(A,B,n, result)
            implicit none
            integer :: n,i
            real, intent(in) :: A(:), B(:)
            real, intent(in out) :: result(:)
        end subroutine vectors_difference
    end interface

    do i=1,n
        call multiplication_matrix_vector(A,X,n,tmp)
        call vectors_difference(tmp, B,n,result)
    enddo
end subroutine calculate_f


!Разность векторов
!A и B - не исходные данные решения СЛАУ, а просто переменные
subroutine vectors_difference(A,B,n, result)
    implicit none
    integer :: n,i
    real, intent(in) :: A(:), B(:)
    real, intent(in out) :: result(:)
    do i=1,n
        result(i) = A(i) - B(i)
    enddo
end subroutine vectors_difference

subroutine multiplication_vector_number(vector, number,n, result)
    implicit none
    integer :: n,i
    real :: number
    real, intent(in) :: vector(:)
    real, intent(in out) :: result(:)

    do, i=1,n
        result(i) = vector(i) * number
    enddo
end subroutine

subroutine check_result(vector,n, result)
    implicit none
    integer :: n,i,number,result
    real, intent(in) :: vector(:)
    result = 1
    do, i=1,n
        if (vector(i) > 0.001 .or. vector(i) < -0.001) then
            result = 0
        end if
    enddo
end subroutine check_result


