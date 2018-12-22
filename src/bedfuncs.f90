!==============================================================================================================
! functions and subroutines used for bedfiles
!==============================================================================================================
!
!	https://www.cog-genomics.org/plink2/formats
!	00 Homozygous for first allele in .bim file
!	01 Missing genotype
!	10 Heterozygous
!	11 Homozygous for second allele in .bim file
!
!==============================================================================================================

!==============================================================================================================
  module bedfuncs

    implicit none

    contains

    function matvec(a,b) result(c)
    implicit none
    external dgemm
    real*8, dimension(:,:), intent(in)  :: a
    real*8, dimension(:), intent(in)  :: b
    real*8 :: c(size(a,1))
    call dgemm('n', 'n', size(a,1), 1, size(a,2), 1.0D0, a, size(a,1), b, size(a,2), 0.0D0, c, size(a,1))
    end function matvec

    function raw2real(n,nbytes,raw) result(w)
    implicit none
    integer, intent(in) :: nbytes,n
    integer, parameter :: byte = selected_int_kind(1)
    integer(byte), intent(in) :: raw(nbytes)
    integer*4 :: i,j,k,rawbits
    real*8 :: w(n)
    real*8, dimension(4) :: rawcodes
    rawcodes = (/ 0.0D0, 3.0D0, 1.0D0, 2.0D0 /)
    ! 00 01 10 11
    w=0.0D0
    k=0
    do i=1,nbytes
      do j=0,6,2
        k = k + 1
        rawbits = ibits(raw(i), j, 2)
        w(k) = rawcodes(rawbits+1)
        if (k==n) exit
      enddo
      if (k==n) exit
    enddo
    end function raw2real

    end module bedfuncs
!==============================================================================================================

!==============================================================================================================
  subroutine ginv(a,n,tol,rank)
!==============================================================================================================
! returns generalized inverse of x(n,n). tol is working zero.
! rework of rohan fernando's f77 subroutine by i. misztal 05/05/87-05/23/00

 implicit none
 integer n, rank
 real*8 :: a(n,n),w(n),tol
 integer i,ii,j

 rank=n
 do i=1,n
    do  j=1,i-1
         a(i:n,i)=a(i:n,i)-a(i,j)*a(i:n,j)
    enddo
    if (a(i,i).lt.tol) then
         a(i:n,i)=0.0
      rank=rank-1
      else
        a(i,i)=sqrt(a(i,i))
        a(i+1:n,i)=a(i+1:n,i)/a(i,i)
    endif
 enddo

 do i=1,n
    if (a(i,i).eq.0.) then
         a(i+1:n,i)=0
      else
         a(i,i)=1.0/ a(i,i)
         w(i+1:n)=0
         do  ii=i+1,n
             w(ii:n)=w(ii:n)-a(ii:n,ii-1)*a(ii-1,i)
             if (a(ii,ii).eq.0.) then
                 a(ii,i)=0.
               else
                 a(ii,i)=w(ii)/a(ii,ii)
             endif
         enddo
     endif
 enddo

 do j=1,n
    do i=j,n
  a(i,j)=dot_product(a(i:n,j),a(i:n,i))
    enddo
 enddo

 do i=1,n
    a(i,i+1:n)=a(i+1:n,i)
 enddo
 end subroutine ginv
!==============================================================================================================

!==============================================================================================================
  subroutine bed2raw(m,cls,nbytes,append,fnBED,fnRAW)
!==============================================================================================================
  use bedfuncs
  implicit none
  integer*4 :: m,cls(m),nbytes,append
  character(len=1000) :: fnBED,fnRAW
  integer*4, parameter :: byte = selected_int_kind(1)
  integer(byte) :: raw(nbytes), magic(3)
  integer*4 :: i,offset,nchar
  integer, parameter :: k14 = selected_int_kind(14)
  integer (kind=k14) :: pos14, nbytes14, offset14,i14
  offset=3
  nchar=index(fnBED, '.bed')
  open(unit=13, file=fnBED(1:(nchar+3)), status='old', access='stream', form='unformatted', action='read')
  nchar=index(fnRAW, '.raw')
  if (append==0) open(unit=14, file=fnRAW(1:(nchar+3)), status='new', access='stream', action='write')
  if (append==1) open(unit=14, file=fnRAW(1:(nchar+3)), status='old', access='stream', action='write', position='append')
  nbytes14 = nbytes
  offset14 = offset
  read(13) magic
  do i=1,m
    if(cls(i)==1) then
    i14=i
    pos14 = 1 + offset14 + (i14-1)*nbytes14
    read(13, pos=pos14) raw
      write(14) raw
    endif
  enddo
  close(unit=13)
  close(unit=14)
  end subroutine bed2raw
!==============================================================================================================

!==============================================================================================================
  subroutine readbed(n,nr,rws,nc,cls,W,nbytes,fnRAW)
!==============================================================================================================
  use bedfuncs
  implicit none
  integer*4 :: n,nr,nc,rws(nr),cls(nc),nbytes
  real*8 :: W(nr,nc),gr(n)
  character(len=1000) :: fnRAW
  integer*4, parameter :: byte = selected_int_kind(1)
  integer(byte) :: raw(nbytes)
  integer*4 :: i, stat,nchar,offset
  offset=0
  nchar=index(fnRAW, '.bed')
  if(nchar>0) offset=3
  if(nchar==0) nchar=index(fnRAW, '.raw')
  open(unit=13, file=fnRAW(1:(nchar+3)), status='old', access='direct', form='unformatted', recl=nbytes)
  W=0.0D0
  do i=1,nc
    read(13, iostat=stat, rec=cls(i)) raw
    gr = raw2real(n,nbytes,raw)
    where(gr==3.0D0) gr=0.0D0
    W(1:nr,i) = gr(rws)
  enddo
  close(unit=13)
  end subroutine readbed
!==============================================================================================================

!==============================================================================================================
  subroutine ltgebed(n,m,rws,cls,nbytes,fnRAW,varG,varE,gebv,Effect,SNPvar,nI,tol)
!==============================================================================================================
  use bedfuncs
  implicit none
  integer :: i,n,rws(n),cls(n),nbytes,m,ncw,irank
  real*8 :: G(n,n),Caa(n,n),nI(n,n),W(n,m),tol,WGa(m,n),varG,varE,One(n),gebv(n),SNPvar(m),Effect(m),P(m)
  character(len=1000) :: fnRAW
  W = 0.0D0
  G = 0.0D0
  do i=1,m,m
    if((i+m-1)<m) ncw = size(cls(i:(i+m-1)))
    if((i+m-1)>=m) ncw = size(cls(i:m))
    call readbed(n,n,rws,ncw,cls(i:(i+ncw-1)),W(:,1:ncw),nbytes,fnRAW)
  enddo
  P = (sum(W, dim=1)/(2.0*dble(n)))
  do i=1,n
    W(i,1:m) = (W(i,1:m) - 2.0*P)/sqrt(dble(m)*2.0*P*(1.0-P))
  end do
  G = matmul(W,transpose(W))
  irank=0
  call ginv(G,n,tol,irank)
  One = 1.0D0
  Caa = 0.0D0
  Caa = nI + (G*varE/varG)
  irank=0
  call ginv(Caa,n,tol,irank)
  Caa = Caa*varE
  WGa = matmul(transpose(W),G)
  SNPvar = matvec(WGa*transpose(W),One)*varG - matvec(matmul(matmul(WGa,Caa),G)*transpose(W),One)
  Effect = matvec(WGa,gebv)
  end subroutine ltgebed
!==============================================================================================================
