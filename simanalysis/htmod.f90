! // htmod.f90 - sim analysis subroutines
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
! */
MODULE HTMOD
use MKL_DFTI
  IMPLICIT NONE

  CONTAINS

  FUNCTION GET_HT_EXTREME(V)
  double complex, dimension(:) :: V
  double precision :: startval,endval,midval,endptsavg
  integer,dimension(1) :: GET_HT_EXTREME
  integer :: mididx, endidx

  endidx = size(V)
  mididx = endidx / 2
  startval = DIMAG(V(1))
  midval = DIMAG(V(mididx))
  endval = DIMAG(V(endidx))
  endptsavg = (startval + endval) / 2.0d0
  if (endptsavg > midval) then
    GET_HT_EXTREME = MINLOC(DIMAG(V))
  else
    GET_HT_EXTREME = MAXLOC(DIMAG(V))
  endif
END FUNCTION GET_HT_EXTREME

SUBROUTINE GETTD(fname,fnamelen,ttdarray,maxesPre,maxesPost,meansPre,meansPost,nMax,nElem, &
globalttd,messPre,microPre,protPre,messPost,microPost,protPost,hdr)
integer,intent(in) :: fnamelen,nMax,hdr
character(len=fnamelen),intent(in) :: fname
DOUBLE PRECISION,intent(out) :: messPre,microPre,protPre,messPost,microPost,protPost
double precision,dimension(nMax),intent(out) :: ttdarray,maxesPre,meansPre,maxesPost,meansPost
integer,intent(out) :: nElem
double precision,intent(out) :: globalttd
integer :: ttdIndex(1)
integer :: nval,nrec,i
integer,dimension(1) :: dumdum
integer :: ios
double complex,dimension(:,:),allocatable :: A,B
double precision,dimension(:),allocatable :: T
double precision,dimension(:,:),allocatable :: SYNTHRATES
logical :: fileexists

nElem = 0
!print*,fname
inquire(file=trim(fname),exist=fileexists)
if (fileexists .eq. .true.) then
open(unit=2322,file=trim(fname),action='read',iostat=ios)
call getData(2322,A,T,SYNTHRATES,ios,hdr)
if (ios .ne. 0) then
  nElem = 0
  close(2322)
  print*,fname
return
endif
close(2322)
ALLOCATE(B(size(A,2),size(A,1)))

call DO_HILB_TRANS(A,B,ttdIndex)
nElem = size(B,2)
globalttd = T(ttdIndex(1))
CALL getSynthStats(SYNTHRATES,T,ttdindex(1),messPre,microPre,protPre,messPost,microPost,protPost)

do i=1,nElem

  dumdum = GET_HT_EXTREME(B(:,i))
  !print*,dtemp
  ttdarray(i) = T(dumdum(1))
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
DEALLOCATE(A,B,T,SYNTHRATES)
else
	print*,fname
	return
endif
END SUBROUTINE GETTD


SUBROUTINE getSynthStats(SYNTHRATES,T,ttdindex,messPre,microPre,protPre,messPost,microPost,protPost)
  DOUBLE PRECISION,DIMENSION(:,:),intent(in) :: SYNTHRATES
  DOUBLE PRECISION,DIMENSION(:),intent(in) :: T
  INTEGER,intent(in) :: ttdindex
  DOUBLE PRECISION,intent(out) :: messPre,microPre,protPre,messPost,microPost,protPost
  DOUBLE PRECISION :: TspanPre,TspanPost
  INTEGER :: N

  N=size(SYNTHRATES,2)

  TspanPre = T(ttdindex)
  TspanPost =  T(N)-TspanPre

  messPre = SYNTHRATES(1,ttdindex)/TspanPre
  microPre = SYNTHRATES(2,ttdindex)/TspanPre
  protPre = SYNTHRATES(3,ttdindex)/TspanPre


  messPost = (SYNTHRATES(1,N)-SYNTHRATES(1,ttdindex))/TspanPost
  microPost = (SYNTHRATES(2,N)-SYNTHRATES(2,ttdindex))/TspanPost
  protPost = (SYNTHRATES(3,N)-SYNTHRATES(3,ttdindex))/TspanPost

END SUBROUTINE getSynthStats

SUBROUTINE writeMtx(fname,M,real_1_imag_0)
  integer :: numrows,numcols,i,j,real_1_imag_0
  double complex,dimension(:,:) :: M
  CHARACTER(LEN=30) :: rowfmt
  character(len=*) :: fname


  numrows=size(M,1)
  numcols=size(M,2)
  write(rowfmt,'(A,I4,A)') '(',numcols,'(E))'
  open(unit=231,file=fname,action='write',status='replace')!,recl=(7*numcols + 50))
  if(real_1_imag_0 .eq. 1) then
    do i=1,numrows
      write(231,FMT=rowfmt) (DREAL(M(i,j)),j=1,numcols)
      
    enddo
  else
    do i=1,numrows
      write(231,FMT=rowfmt) (DIMAG(M(i,j)),j=1,numcols)
      
    enddo
  endif


  close(231)
END SUBROUTINE writeMtx

SUBROUTINE normalizeCols(B)
  double complex,dimension(:,:) :: B
  double precision :: cn
  integer :: j
  do j=1,size(B,2)
    cn = sum(abs(b(:,j))**2)
    cn = sqrt(cn)
    b(:,j)=b(:,j) / cn
  enddo

END SUBROUTINE normalizeCols

SUBROUTINE normalizeRows(A)
  double complex,dimension(:,:) :: A
  double precision :: rownorm,rn
  integer :: i,j
  !Skip dnrm2 and do the safe way
  do i=1,size(A,1)
    rn = 0.0d0
    !rownorm = dnrm2(size(A,2),A(i,:),1)
    !print*,rownorm
    do j=1,size(A,2)
      rn = rn + abs(A(i,j))**2
    enddo
    rn = sqrt(rn)
    A(i,:)=A(i,:)/rn
  enddo

END SUBROUTINE normalizeRows

FUNCTION getNumRecords(unitno,hdr)
  integer :: unitno,getNumRecords,hdr
  rewind(unitno)
  getNumRecords = 0

  DO
    read(unitno,*,end=8243)
    getNumRecords = getNumRecords + hdr
  END DO

8243 CONTINUE
  rewind(unitno)


END FUNCTION getNumRecords


FUNCTION getNumFields(unitno,skipheader)
  integer :: unitno,getNumFields,skipheader
  CHARACTER(1) :: DUMMY_CHAR
  CHARACTER(1) :: DELIM_CHAR = CHAR(9)
  
  getNumFields = 1
  rewind(unitno)
  if(skipheader .ne. 0)   read(unitno,*)
!SKIP FIRST LINE WHICH MAY BE HEADER

  DO
    read(unitno,'(A)',advance='no',eor=18294) DUMMY_CHAR
    if (DUMMY_CHAR .eq. DELIM_CHAR) then
      getNumFields = getNumFields + 1
    endif
  END DO


18294 CONTINUE

  rewind(unitno)

END FUNCTION getNumFields

SUBROUTINE getData(unitno,A,T,SYNTHRATES,ios,hdr)
  integer :: unitno, nRecs,nVals,nHeaderEntries,nProts,nSynths,hdr
  double complex,dimension(:,:),allocatable :: A
  double precision,dimension(:,:),allocatable :: dummy
  double precision,dimension(:),allocatable :: T
  double precision,dimension(:,:),allocatable :: SYNTHRATES
  integer :: ios

  nHeaderEntries = getNumFields(unitno,0)
  nProts = nHeaderEntries-1

  nRecs = getNumRecords(unitno,hdr) -1
  nVals = getNumFields(unitno,hdr)
  nSynths = nVals - nHeaderEntries

  if (ALLOCATED(A)) then
    DEALLOCATE(A)
  endif

  ALLOCATE(A(nProts,nRecs),dummy(nVals,nRecs))
  ALLOCATE(T(nRecs))
  if (nSynths .gt. 0) ALLOCATE(SYNTHRATES(nSynths,nRecs))

  rewind(unitno)
  if (hdr .ne. 0) then
    read(unitno,*,iostat=ios)
  endif
  if (ios .ne. 0) then 
    DEALLOCATE(A,dummy,T,SYNTHRATES)
    return
  endif
  read(unitno,*,iostat=ios) dummy
  if (ios .ne. 0) then 
    DEALLOCATE(A,dummy,T,SYNTHRATES)
    return
  endif
  rewind(unitno)
!slice up dummy
  T=dummy(1,:)
  A=DCMPLX(dummy(2:nProts+1,:),0.0d0)
  SYNTHRATES=dummy(nProts+2:nProts+2+nSynths,:)
  DEALLOCATE(dummy)
END SUBROUTINE getData

SUBROUTINE FFTROWS(A,fwd)
  double complex,dimension(:,:) :: A
  double precision :: scaleFactor 
  integer :: fwd
  TYPE(DFTI_DESCRIPTOR),POINTER :: dhandl
  integer :: status

  scaleFactor = 1.0d0 / sqrt(DBLE(size(A,1)))

  status = DftiCreateDescriptor(dhandl,DFTI_DOUBLE,DFTI_COMPLEX,1,size(A,1))
  status = DftiSetValue(dhandl,DFTI_NUMBER_OF_TRANSFORMS,size(A,2))
  status = DftiSetValue(dhandl,DFTI_INPUT_DISTANCE,size(A,1))
  status = DftiSetValue(dhandl,DFTI_FORWARD_SCALE,scaleFactor)
  status = DftiSetValue(dhandl,DFTI_BACKWARD_SCALE,scaleFactor)
  status = DftiCommitDescriptor(dhandl)
  if (fwd .eq. 1 ) then
    status = DftiComputeForward(dhandl,A(:,1))
  else
    status = DftiComputeBackward(dhandl,A(:,1))
  endif 

  status = DftiFreeDescriptor(dhandl)


END SUBROUTINE FFTROWS



SUBROUTINE DO_HILB_TRANS(A,B,globalttdindex)
  double complex,dimension(:,:) :: A,B
  double complex,dimension(:),allocatable :: Y
  integer :: n2
  integer :: globalttdindex(1)

  !call normalizeRows(B)
  B = transpose(A)

  allocate(Y(size(B,1)))

  call normalizeCols(B)
  Y=sum(B,2)

  n2 = (size(B,1) / 2) + 1
  call FFTROWS(B,1)
    B(2:n2,:) = B(2:n2,:) * 2.0d0
    B(n2+1:size(B,1),:) = 0.0d0
  call FFTROWS(B,0)

  call FFTROW(Y,1)
  Y(2:n2) = Y(2:n2) * 2.0d0
  Y(n2+1:size(Y)) = 0.0d0
  call FFTROW(Y,0)

  !There may be a reflection because the discrete HT can be may be interpreted differently 
  !for different FFT constants
  globalttdindex = GET_HT_EXTREME(Y)

  DEALLOCATE(Y)

END SUBROUTINE DO_HILB_TRANS



SUBROUTINE FFTROW(B,fwd)
  double complex,dimension(:) :: B
  double precision :: scaleFactor 
  integer :: fwd
  TYPE(DFTI_DESCRIPTOR),POINTER :: dhandl
  integer :: status

  scaleFactor = 1.0d0 / sqrt(DBLE(size(B)))
  
  status = DftiCreateDescriptor(dhandl,DFTI_DOUBLE,DFTI_COMPLEX,1,size(B))
  status = DftiSetValue(dhandl,DFTI_FORWARD_SCALE,scaleFactor)
  status = DftiSetValue(dhandl,DFTI_BACKWARD_SCALE,scaleFactor)
  status = DftiCommitDescriptor(dhandl)
  if (fwd .eq. 1 ) then
    status = DftiComputeForward(dhandl,B)
  else
    status = DftiComputeBackward(dhandl,B)
  endif 

  status = DftiFreeDescriptor(dhandl)
  

END SUBROUTINE FFTROW

END MODULE HTMOD