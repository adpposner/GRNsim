! // mainfort.f90 - a fortran-only interface for testing the analysis method
! /*
    
!     Copyright (C) 2018, Russell Posner

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
PROGRAM MAIN
use MKL_DFTI
use htmod
IMPLICIT NONE

integer:: fnamelen,hdr
character(len=256):: fname
DOUBLE PRECISION :: messPre,microPre,protPre,messPost,microPost,protPost
double precision,dimension(:),allocatable :: ttdarray,maxesPre,meansPre,maxesPost,meansPost
integer :: nElem
double precision :: globalttd
integer :: ttdIndex(1)
integer :: nval,nrec,i,dumdum,status
integer :: ios
double complex,dimension(:,:),allocatable :: A,B
double precision,dimension(:),allocatable :: T
double precision,dimension(:,:),allocatable :: SYNTHRATES
logical :: fileexists

call GETARG(1,fname,status)
nElem = 0
!print*,trim(fname)
inquire(file=trim(fname),exist=fileexists)
if (fileexists .eq. .true.) then
open(unit=2322,file=trim(fname),action='read',iostat=ios)
call getData(2322,A,T,SYNTHRATES,ios,1)
if (ios .ne. 0) then
  nElem = 0
  close(2322)
  print*,fname
  print*,"IOSERROR"
stop
endif
close(2322)
ALLOCATE(B(size(A,2),size(A,1)))

call DO_HILB_TRANS(A,B,ttdIndex)
nElem = size(B,2)
globalttd = T(ttdIndex(1))
CALL getSynthStats(SYNTHRATES,T,ttdindex(1),messPre,microPre,protPre,messPost,microPost,protPost)
allocate(ttdarray(nElem),maxesPre(nElem),meansPre(nElem),maxesPost(nElem),meansPost(nElem))
do i=1,nElem
  dtemp = MINLOC(DIMAG(B(:,i)),1)

  ttdarray(i) = T(dtemp)
  maxesPre(i) = MAXVAL(DREAL(A(i,1:ttdindex(1))),1)
  maxesPost(i) = MAXVAL(DREAL(A(i,ttdindex(1):size(A,2))),1)

  meansPre(i) = sum(DREAL(A(i,1:ttdindex(1))))/DBLE(ttdindex(1))
  meansPost(i) = sum(DREAL(A(i,ttdindex(1):size(A,2))))/DBLE(size(A,2)-ttdindex(1))
enddo
! print*,maxes
! print*,means
! print*,ttdarray
! print*,nMax
! print*,nElem
! print*,globalttd
!print*,B
call WRITEMAT_DBLCOMPLEX(B,"ht.txt")

DEALLOCATE(A,B,T,SYNTHRATES)
DEALLOCATE(ttdarray,maxesPre,maxesPost,meansPre,meansPost)
else
	print*,fname
	stop
endif

contains

subroutine WRITEMAT_DBLCOMPLEX(data,filename)
  double complex,dimension(:,:) :: data
  character(len=3), dimension(size(data,1),size(data,2)) :: imag_unit
  character(len=45) :: rowfmt
  integer :: i,j
  character(len=*) :: filename
  integer :: numrows,numcols
  numrows = size(data,1)
  numcols = size(data,2)
  
  write(rowfmt,'(A,I4,A)') '(',numcols,'(F10.8,SP,F10.8,"j",1x))'
print*,rowfmt
  where(aimag(data)<0.) 
    imag_unit = '-i*'
elsewhere 
  imag_unit = '+i*'
  end where

  open(unit=2414,file=trim(filename),action='write')!,recl =(15*numcols)+10)
  
  do i=1,numrows
      write(2414,rowfmt) (data(i,j),j=1,numcols)
      
  enddo


  close(2414)

end subroutine WRITEMAT_DBLCOMPLEX

END PROGRAM MAIN