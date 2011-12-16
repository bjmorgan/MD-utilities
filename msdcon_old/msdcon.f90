program msdcalc

! Last altered 25/10/12 
! Converted to Fortran 90, 26.3.10

implicit none 

integer, parameter :: nmsdcorrmax=1000
integer, parameter :: nummax=2000
integer, parameter :: nspmax=3
 
double precision :: xdisp(nummax),xdispstore(nummax,nmsdcorrmax), &
                    ydisp(nummax),ydispstore(nummax,nmsdcorrmax), &
                    zdisp(nummax),zdispstore(nummax,nmsdcorrmax)
double precision :: xmsd(nummax,0:nmsdcorrmax), &
                    ymsd(nummax,0:nmsdcorrmax), &
                    zmsd(nummax,0:nmsdcorrmax)
double precision :: xdispabs(nummax,0:nmsdcorrmax), &
                    ydispabs(nummax,0:nmsdcorrmax), &
                    zdispabs(nummax,0:nmsdcorrmax)
integer :: norm(nummax,0:nmsdcorrmax)
double precision :: xds(0:nmsdcorrmax,nspmax,nspmax), &
                    yds(0:nmsdcorrmax,nspmax,nspmax), &
                    zds(0:nmsdcorrmax,nspmax,nspmax)
integer :: normtot(0:nmsdcorrmax)
double precision :: xtot(0:nmsdcorrmax, nspmax), &
                    ytot(0:nmsdcorrmax, nspmax), &
                    ztot(0:nmsdcorrmax, nspmax)
double precision :: msd(0:nmsdcorrmax, nspmax), &
                    msdx(0:nmsdcorrmax, nspmax), &
                    msdy(0:nmsdcorrmax, nspmax), &
                    msdz(0:nmsdcorrmax, nspmax), &
                    msdxabs(0:nmsdcorrmax, nspmax), &
                    msdyabs(0:nmsdcorrmax, nspmax), &
                    msdzabs(0:nmsdcorrmax, nspmax), &
                    msdcoll(0:nmsdcorrmax, nspmax, nspmax)
double precision :: suma(nspmax), sumx(nspmax), sumy(nspmax), sumz(nspmax)
double precision :: sumxabs(nspmax), sumyabs(nspmax), sumzabs(nspmax)
integer :: nmsdlength,nmsdcalltime,mcorrtime
integer :: ntype(nummax),numspc(nspmax)
integer :: num
integer :: nstep,nrun
double precision :: dtime
integer :: ncation1,ncation2,ncation3,nanion
double precision :: z(nspmax)
integer :: nspecies
double precision :: rn, time, work, wnernst

character(len=200) :: filein
logical :: overflow
logical :: restart
logical :: endrun
logical :: readfrominpt_log

integer :: i, j, i1, n, ns, ns2, i2, ispec, nt, ipt, ip, ip1, ip2, n1, n2
integer :: ipoint
double precision :: dum1, dum2

real :: bin
integer :: nbin
integer, parameter :: DISPIO = 11

character(len=20) :: rstfile,filename
rstfile='msdrst.dat'

restart=.false.
overflow=.false.
endrun=.false.
readfrominpt_log=.false.

!==================================================
! Get parameters
!==================================================


open(10,file='msd.inpt')
read(10,'(a)') filein
open(DISPIO,file=filein,status='old',form='formatted')
read(10,*) num
read(10,*) nanion
read(10,*) ncation1
read(10,*) ncation2
read(10,*) nmsdlength
read(10,*) nspecies

if (nspecies.gt.nspmax.or.num.gt.nummax.or. &
nmsdlength.gt.nmsdcorrmax) then
    write(6,*) '***Problems with array dimensions***'
    stop
endif

if (num /= (nanion+ncation1+ncation2)) then
    write(6,*) '*** Problems with ion numbers ***'
    stop
endif
do ns=1,nspecies
    read(10,*) z(ns)
enddo

do ns=1,nspecies
    filename='msd'//char(ns+48)//'.dat'
    open (52+ns,file=filename)
    filename='msdx'//char(ns+48)//'.dat'
    open(80+ns*3,file=filename)
    filename='msdy'//char(ns+48)//'.dat'
    open(81+ns*3,file=filename) 
    filename='msdz'//char(ns+48)//'.dat'
    open(82+ns*3,file=filename) 
    do ns2=ns,nspecies
        filename='msdcollect'//char(ns+48)//char(ns2+48)//'.dat'
        open (70+(ns*nspecies)+ns2,file=filename)
    enddo
enddo
open (98,file='work.dat')
open (99,file='nernst.dat')
  
!==================================================
! Zero arrays
!==================================================

xdisp = 0.0d0
ydisp = 0.0d0
zdisp = 0.0d0

xdispstore = 0.0d0
ydispstore = 0.0d0
zdispstore = 0.0d0
       
xmsd = 0.0d0
ymsd = 0.0d0
zmsd = 0.0d0

xdispabs = 0.0d0
ydispabs = 0.0d0
zdispabs = 0.0d0

norm = 0

xds = 0.0d0
yds = 0.0d0
zds = 0.0d0
normtot =0
 
read(10,*) restart

if (restart) then

    open(12,file=rstfile,status='OLD',form='formatted')

    do i=0,nmsdlength
        do n=1,num
            read(12,*) xmsd(n,i)
            read(12,*) ymsd(n,i)
            read(12,*) zmsd(n,i)
            read(12,*) norm(n,i)
        enddo
        do i1=1,nspecies
            do i2=1,nspecies
                read(12,*) xds(i,i1,i2)
                read(12,*) yds(i,i1,i2)
                read(12,*) zds(i,i1,i2)
            enddo
        enddo
        read(12,*) normtot(i)
    enddo
    close(12)

endif

read(10,*) readfrominpt_log

if (readfrominpt_log) then
    read (10,*) nmsdcalltime
    read (10,*) dtime
    read (10,*) nrun

    do n=1,nanion
        ntype(n)=1
    enddo

    if (ncation1.ne.0) then
        do n=nanion+1,nanion+ncation1
            ntype(n)=2
        enddo
    endif

    if (ncation2.ne.0) then
        do n=nanion+ncation1+1,num
            ntype(n)=3
        enddo
    endif

    read (DISPIO,*) nbin,bin,nbin
    do i=1,num
        read(DISPIO,*) nbin
    enddo
else
    read (DISPIO,*) nmsdcalltime,dtime,nrun
    do i=1,num
        read(DISPIO,*) ntype(i)
    enddo
endif

do ispec=1,nspecies
    numspc(ispec) = 0
enddo

do n=1,num
    ispec = ntype(n)
    numspc(ispec) = numspc(ispec) + 1
enddo

close(10)

nstep=0
mcorrtime=1

endrun = .false.

do while ( endrun == .false. )

    nstep=nstep+nmsdcalltime

    if (nstep.ge.nrun) endrun=.true.

    do i=1,num
        read(DISPIO,*) xdisp(i), ydisp(i), zdisp(i)
    enddo

    print *,'entering msdcalc',nstep
    do i=1,num
        xdispstore(i,mcorrtime) = 0.0d0
        ydispstore(i,mcorrtime) = 0.0d0
        zdispstore(i,mcorrtime) = 0.0d0
    enddo

    ! accumulate displacements

    do j=1,mcorrtime
        do i=1,num
            xdispstore(i,j) = xdispstore(i,j) + xdisp(i)
            ydispstore(i,j) = ydispstore(i,j) + ydisp(i)
            zdispstore(i,j) = zdispstore(i,j) + zdisp(i)
        enddo
    enddo

    if (overflow) then

        do j=mcorrtime+1,nmsdlength
            do i=1,num
                xdispstore(i,j) = xdispstore(i,j) + xdisp(i)
                ydispstore(i,j) = ydispstore(i,j) + ydisp(i)
                zdispstore(i,j) = zdispstore(i,j) + zdisp(i)
            enddo
        enddo
    endif

    ! zero displacement arrays

    do i=1,num
        xdisp(i) = 0.0d0
        ydisp(i) = 0.0d0
        zdisp(i) = 0.0d0
    enddo
     
    do j=1,mcorrtime
        nt=mcorrtime-j
        do ipt=1,nspecies
             xtot(j,ipt) = 0.0d0
             ytot(j,ipt) = 0.0d0
             ztot(j,ipt) = 0.0d0
        enddo
        do i=1,num
            ipoint=ntype(i)
              
            xtot(j,ipoint) = xtot(j,ipoint) + xdispstore(i,j)
            ytot(j,ipoint) = ytot(j,ipoint) + ydispstore(i,j)
            ztot(j,ipoint) = ztot(j,ipoint) + zdispstore(i,j)

            norm(i,nt) = norm(i,nt) + 1

            xmsd(i,nt) = xmsd(i,nt) + xdispstore(i,j)**2
            ymsd(i,nt) = ymsd(i,nt) + ydispstore(i,j)**2
            zmsd(i,nt) = zmsd(i,nt) + zdispstore(i,j)**2
 
            xdispabs(i,nt) = xdispabs(i,nt) + sign(xdispstore(i,j)**2, xdispstore(i,j))
            ydispabs(i,nt) = ydispabs(i,nt) + sign(ydispstore(i,j)**2, ydispstore(i,j))
            zdispabs(i,nt) = zdispabs(i,nt) + sign(zdispstore(i,j)**2, zdispstore(i,j))
       enddo
        normtot(nt)=normtot(nt)+1
        do i1=1,nspecies
            do i2=i1,nspecies
                xds(nt,i1,i2) = xds(nt,i1,i2) + xtot(j,i1) * xtot(j,i2)
                yds(nt,i1,i2) = yds(nt,i1,i2) + ytot(j,i1) * ytot(j,i2)
                zds(nt,i1,i2) = zds(nt,i1,i2) + ztot(j,i1) * ztot(j,i2)
            enddo
        enddo

    enddo
     
    if (overflow) then
     
        do j=mcorrtime+1,nmsdlength

            nt=mcorrtime-j+nmsdlength

            do ipt=1,nspecies
                xtot(j,ipt)=0.0d0
                ytot(j,ipt)=0.0d0
                ztot(j,ipt)=0.0d0
            enddo

            do i=1,num
                ipoint=ntype(i)

                xtot(j,ipoint) = xtot(j,ipoint) + xdispstore(i,j)
                ytot(j,ipoint) = ytot(j,ipoint) + ydispstore(i,j)
                ztot(j,ipoint) = ztot(j,ipoint) + zdispstore(i,j)

                norm(i,nt) = norm(i,nt) + 1

                xmsd(i,nt) = xmsd(i,nt) + xdispstore(i,j)**2
                ymsd(i,nt) = ymsd(i,nt) + ydispstore(i,j)**2
                zmsd(i,nt) = zmsd(i,nt) + zdispstore(i,j)**2
 
                xdispabs(i,nt) = xdispabs(i,nt) + sign( xdispstore(i,j)**2, xdispstore(i,j) )
                ydispabs(i,nt) = ydispabs(i,nt) + sign( ydispstore(i,j)**2, ydispstore(i,j) )
                zdispabs(i,nt) = zdispabs(i,nt) + sign( zdispstore(i,j)**2, zdispstore(i,j) )
            enddo
 
            normtot(nt)=normtot(nt)+1
     
            do i1=1,nspecies
                do i2=i1,nspecies
      
                    xds(nt,i1,i2) = xds(nt,i1,i2)+ xtot(j,i1)*xtot(j,i2)
                    yds(nt,i1,i2) = yds(nt,i1,i2)+ ytot(j,i1)*ytot(j,i2)
                    zds(nt,i1,i2) = zds(nt,i1,i2)+ ztot(j,i1)*ztot(j,i2)
                enddo
            enddo

        enddo

    endif

    !========================================================================
    ! Update array counters
    !========================================================================
    if (mod(float(mcorrtime),float(nmsdlength)).eq.0) then
        overflow=.true.
    endif

    mcorrtime=int(mod(float(mcorrtime),float(nmsdlength)))
    mcorrtime=mcorrtime+1

enddo

close(DISPIO)

!=========================================================
! Write out restart file
!=========================================================

msd = 0.0d0
msdcoll = 0.0d0
  
!   Writing out values for restart file

open (14,file=rstfile,form='formatted')

do i=0,nmsdlength
    do n=1,num
        write(14,*) xmsd(n,i)
        write(14,*) ymsd(n,i)
        write(14,*) zmsd(n,i)
        write(14,*) norm(n,i)
    enddo
    do i1=1,nspecies
        do i2=1,nspecies
            write(14,*) xds(i,i1,i2),"xds"
            write(14,*) yds(i,i1,i2),"yds"
            write(14,*) zds(i,i1,i2),"zds"
        enddo
    enddo

    write(14,*) normtot(i)
enddo

close(14)

!======================================================================
! Average over components and number of molecules
!======================================================================

do i=1, num
    do j=0, nmsdlength
        if (norm(i,j) > 0) then

            xmsd(i,j) = xmsd(i,j) / float(norm(i,j))
            ymsd(i,j) = ymsd(i,j) / float(norm(i,j))
            zmsd(i,j) = zmsd(i,j) / float(norm(i,j))

            xdispabs(i,j) = xdispabs(i,j) / float(norm(i,j))
            ydispabs(i,j) = ydispabs(i,j) / float(norm(i,j))
            zdispabs(i,j) = zdispabs(i,j) / float(norm(i,j))

        endif
    enddo
enddo

do j=0,nmsdlength
    do n1=1,nspecies
        do n2=1,nspecies
            if (normtot(j) > 0) then
                rn=sqrt(float(numspc(n1)*numspc(n2)))*float(normtot(j))

                xds(j,n1,n2) = xds(j,n1,n2)/rn
                yds(j,n1,n2) = yds(j,n1,n2)/rn
                zds(j,n1,n2) = zds(j,n1,n2)/rn
                msdcoll(j,n1,n2) = xds(j,n1,n2)+yds(j,n1,n2)+zds(j,n1,n2)

            endif
        enddo
    enddo
enddo

do i=0,nmsdlength
    do ip=1,nspecies
        suma(ip)=0.0d0
        sumx(ip)=0.0d0
        sumy(ip)=0.0d0
        sumz(ip)=0.0d0
        sumxabs(ip)=0.0d0
        sumyabs(ip)=0.0d0
        sumzabs(ip)=0.0d0
    enddo
    do j=1,num
        ip=ntype(j)
        suma(ip) = suma(ip) + xmsd(j,i) + ymsd(j,i) + zmsd(j,i)
        sumx(ip) = sumx(ip) + xmsd(j,i)
        sumy(ip) = sumy(ip) + ymsd(j,i)
        sumz(ip) = sumz(ip) + zmsd(j,i)

        sumxabs(ip) = sumxabs(ip) + xdispabs(j,i)
        sumyabs(ip) = sumyabs(ip) + ydispabs(j,i)
        sumzabs(ip) = sumzabs(ip) + zdispabs(j,i)
    enddo

    do ip=1,nspecies
        if (numspc(ip) > 0) then
            msd(i,ip) = suma(ip) / float(numspc(ip))

            msdx(i,ip) = sumx(ip) / float(numspc(ip))
            msdy(i,ip) = sumy(ip) / float(numspc(ip))
            msdz(i,ip) = sumz(ip) / float(numspc(ip))

            msdxabs(i,ip) = sumxabs(ip) / float(numspc(ip))
            msdyabs(i,ip) = sumyabs(ip) / float(numspc(ip))
            msdzabs(i,ip) = sumzabs(ip) / float(numspc(ip))
        endif
    enddo
enddo


dum1=0.0
dum2=0.0

do ip=1,nspecies
    write (52+ip,*) dum1,dum2
enddo

write(98,*)dum1,dum2
write(99,*)dum1,dum2

do ip1=1,nspecies
    if (numspc(ip1) > 0) then
        do ip2 = ip1, nspecies
            write(70+(ip1*nspecies)+ip2,*) dum1,dum2
        enddo
    endif
enddo

do i=0,nmsdlength-1
    time=(dble(i)+1)*dble(nmsdcalltime)*dtime*2.418d-5
    do ip=1, nspecies
        write(80+ip*3,*) time, msdx(i, ip)
        write(81+ip*3,*) time, msdy(i, ip)
        write(82+ip*3,*) time, msdz(i, ip)
        write(52+ip,*) time, msd(i, ip)
        write(200+ip,'(4(e10.3,1x))') time, msdxabs( i, ip ), msdyabs( i, ip ), msdzabs( i, ip )
    enddo

    work=0.0d0
    wnernst=0.0d0

    do ip1=1,nspecies
        if (numspc(ip1).gt.0) then
            do ip2=ip1,nspecies
                write(70+(ip1*nspecies)+ip2,*) time, &
                               msdcoll(i,ip1,ip2)
                msdcoll(i,ip1,ip2)=msdcoll(i,ip1,ip2)* &
                                   z(ip1)*z(ip2)* &
                sqrt(float(numspc(ip1)*numspc(ip2)))/float(num)
                if (ip1.ne.ip2) then 
                    work=work+msdcoll(i,ip1,ip2)
                endif
                work=work+msdcoll(i,ip1,ip2)
            enddo
        endif

        wnernst=wnernst+msd(i,ip1)*z(ip1)*z(ip1)*float(numspc(ip1))/float(num)

    enddo
    write(98,*)time,work
    write(99,*)time,wnernst
enddo
  
write (6,*)
write (6,*) '*** Mean squared displacements written out. ***'
write (6,*)

do ns=1,nspecies
     close(102+ns)
enddo

stop
end
