      implicit none
      integer nconf
      integer ntot
      parameter (ntot=1000000)
      real*8 plq(ntot),plqave
      real*8 jack(0:ntot),plqt(ntot),plqterr
      integer i,j,ibin,iconf,nnconf,maxbin,m
      character*80 arg1,arg2,arg3

      call getarg(1,arg1)
c      call getarg(2,arg2)
c      call getarg(3,arg3)
c      read(arg3,*)nconf
      open(22,file=arg1,status="old")
      nconf=0
      do
        nconf=nconf+1
        read(22,*,end=10 )plq(nconf)
      enddo
10    close(22)
      nconf=nconf-1

      plqave=0.0d0
      do i=1,nconf
        plqave=plqave+plq(i)
      enddo
      plqave=plqave/nconf
c      write(*,'(ES24.15)')plqave

*******************
      do ibin=1,nconf/2
      if (mod(nconf,ibin).NE.0) goto 100

      i=1
      iconf=1
      do while (ibin*i .LE. nconf)  ! .LE. -> <=
        plqave=0.0d0
        do j=1+ibin*(i-1),ibin*i
c          write(*,*)ibin,j,iconf
          plqave=plqave+plq(j)
        enddo
        plqave=plqave/ibin
c        write(*,*)plqave
        plqt(iconf)=plqave
        i = i+1
        iconf = iconf+1
      enddo
      nnconf=iconf-1

      plqave=0.0d0
      do iconf = 1,nnconf
        jack(iconf)=plqt(iconf)
        plqave = plqave + plqt(iconf)
      enddo
      jack(0)=plqave/DBLE(nnconf)
      do iconf = 1,nnconf
        jack(iconf)=(jack(0)*dble(nnconf)-jack(iconf))/DBLE(nnconf-1)
      enddo
      plqave=0.0d0
      do iconf = 1,nnconf
        plqave = plqave + jack(iconf)
      enddo
      plqave = plqave/DBLE(nnconf)
      plqterr = 0.0d0
      do iconf = 1,nnconf
        plqterr = plqterr + (plqave-jack(iconf))**2
      enddo
      plqterr = plqterr*DBLE((nnconf-1))/DBLE(nnconf)
      plqterr = dsqrt(plqterr)
      plqave = jack(0)

c      write(*,'("# bin size=",I5," totalconf=",I5)')ibin,nnconf
c      write(*,'("# PLQ :",2F32.16)') plqave,plqterr
c      write(*,'(I11,2F32.14,"  #",I11)')ibin,plqave,plqterr,nnconf
      if (ibin == 1) then
      write(*,'(2ES24.15)') plqave,plqterr
      end if
100   continue
      enddo

      stop
      end
