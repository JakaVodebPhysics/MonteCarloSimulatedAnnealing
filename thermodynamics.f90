!include "num_to_txt.f90"

program main
!Variable declaration
implicit none
integer(kind=8)               :: ni, nit, lt, intv, intvt, intt, inttt
integer                       :: i, j, k, x, y, le, side1, side2, size, count, radij, mrad, range1&
&, range2, answer
double precision              :: u, delez, tempz, tempi, tempk, ls, ena, dva, tri, stiri, jay, c, start, finish
double precision, allocatable :: temp(:), v(:,:), enm(:,:), cvm(:,:)
integer, allocatable          :: s0(:,:), enke0(:), nbar(:), tarr(:), s0out(:,:,:)
character(len=30)             :: str, strl, intname, tempzzi, tempzzd, tempkzi, tempkzd, tempizi&
&, tempizd, lszi, lszd, jayzi, jayzd, czi, czd, lattice, pwd
!2D polaron model with nearest neighbour strain and cutoff Coulomb interaction between polarons located on a uniformly charged plate in the shape of a rhombus
!(periodic boundary conditions) influenced by an external field- Jaka Vodeb
!ENTER THE WORKING DIRECTORY:
pwd=""
!Interface
call CPU_TIME(start)
open (10,file=""//trim(pwd)//"input.txt")
!Do you want to read a premade configuration?
read(10,*) answer
!The lattice topology
read(10,*) lattice
!The initial temperature of the system
read(10,*) tempz
!The final temperature of the system
read(10,*) tempk
!The temperature interval of the system
read(10,*) tempi
!The interaction name
read(10,*) intname
!The nearest-neighbour coupling
read(10,*) jay
!The screened Coulomb interaction magnitude
read(10,*) c
!The radius of the interaction
read(10,*) radij
!The movement radius of the local update scheme
read(10,*) mrad
!Characteristic length of the screened Coulomb
read(10,*) ls
!The number of sites on each side of the paralelogram
read(10,*) side1
read(10,*) side2
!Total number of sites
size=side1*side2
!The number of fractions of polaron occupation and sites in the system to evaluate
read(10,*) le
allocate (nbar(le))
!Enter the fractions to be evaluated
do i=1,le
! i-th fraction
  read(10,*) delez
  nbar(i)=idnint(delez*dble(size))
end do
!At 0.1 temperature relaxation of system measured at 5*10^5 for 30, 1*10^6 for 60 and 5*10^6 for 120 as the number of sites on the side of the rhombic lattice.
!Enter the number of sweeps at each temperature
read(10,*) nit
!Enter the sampling interval in sweeps
read(10,*) intvt
!Enter the thermalization interval in sweeps
read(10,*) inttt
close(10)
lt=idnint(dabs(tempk-tempz)/dabs(tempi))
ena=1.
dva=2.
tri=3.
stiri=4.
allocate (s0(side1,side1))
s0=0
do i = 1, side1
   do j = 1, side1
       s0(i,j) = MINVAL((/ABS(i-j), ABS(i-j+side1), ABS(i-j-side1)/))
   end do
end do
range1 = MAXVAL(s0)
deallocate(s0)
allocate (s0(side2,side2))
do i = 1, side2
   do j = 1, side2
       s0(i,j) = MINVAL((/ABS(i-j), ABS(i-j+side2), ABS(i-j-side2)/))
   end do
end do
range2 = MAXVAL(s0)
deallocate(s0)
allocate (s0(side1,side2))
!Array allocation
allocate (temp(lt))
allocate (v(-range1:range1,-range2:range2))
allocate (enke0(1))
allocate (enm(lt,3))
allocate (cvm(lt,3))
allocate (s0out(lt,side1,side2))
do i=1,lt
  temp(i)=tempz+dble(i)*tempi
end do
do i=-range1,range1
    do j=-range2,range2
        if (lattice == 'tri') then
            u=dsqrt((dble(i)+dble(j)/dva)*(dble(i)+dble(j)/dva)+tri*dble(j)*dble(j)/stiri)
        else if (lattice == 'sq') then
            u=dsqrt(dble(i)*dble(i)+dble(j)*dble(j))
        else
            print *, 'THE LATTICE PARAMETER DOESNT MATCH ANY OF THE AVAILABLE OPTIONS!'
            go to 20
        end if
        if (u>dble(0).and.u<=dble(radij)) then
            if (u<=1.1) then
                v(i,j)=-jay
            else
                v(i,j)=c*dexp(-u/ls)/u
            end if
        else
            v(i,j)=0.
        end if
    end do
end do
!Klicanje novega zacetka izbiranja random stevil
call random_seed
!Poimenovanje indeksov, ki spremljata izvajanje programa
!print *, '          delez'
!Zacetek zanke, ki tece po vseh zeljenih delezih polaronov v sistemu
do i=1,le
    ni=nit*nbar(i)
    intv=intvt*nbar(i)
    intt=inttt*nbar(i)
    !Pisanje mest v zanki stevila polaronov v sistemu z namenom lazjega spremljanja delovanja programa
    !print *, idnint(dble(nbar(i))/dble(size))
    !Pretvorba temperature v zapis za output
    call dec2str (tempz, tempzzi, tempzzd)
    call dec2str (tempk, tempkzi, tempkzd)
    call dec2str (tempi, tempizi, tempizd)
    call dec2str (ls, lszi, lszd)
    call dec2str (jay, jayzi, jayzd)
    call dec2str (c, czi, czd)
    !Spreminjanje velikosti arrayjev polaronov in praznih mest
    deallocate(enke0)
    !deallocate(nicle0)
    allocate (enke0(2*nbar(i)))
    !allocate (nicle0(2*(size-nbar(i))))
    if (answer==1) then
        count=1
        !index=1
        enke0=0
        s0=0
        allocate(tarr(side2))
        open(10,file=""//trim(pwd)//"konf_int_barr_1_radij_24_ls_4_5_side1_&
		&91_side2_104_N_728_ni_30000000_tempz_1_000000000000000E-002_tempk_8_999999999999999E-003_tempi_-1_000000000000000E-0&
		&03.dat")
        do j=1,side1
            read(10,*) tarr
            do k=1,side2
                if (tarr(k)==1) then
                    enke0(count)=j
                    enke0(count+1)=k
                    count=count+2
                    s0(j,k)=1
                !else
                    !nicle0(index)=j
                    !nicle0(index+1)=k
                    !index=index+2
                end if
            end do
        end do
        close(10)
    else if (answer==0) then
        !Dodeljevanje mest delcem v levi zgornji kot resetke
        count=1
        x=1
        y=1
        s0=0
        do while (y<=side2)
            do while (x<=side1)
                if (count<2*nbar(i)) then
                    enke0(count)=x
                    enke0(count+1)=y
                    s0(x,y)=1
                    x=x+1
                    count=count+2
                else
                    go to 30
                end if
            end do
            x=1
            y=y+1
        end do
    else
        !print *, 'WRONG ANSWER!'
        go to 20
    end if
    !Klicanje dinamike polaronov
30  call spini(lt, temp, v, side1, side2, range1, range2, ni, intv, intt, nbar(i), enke0, s0, mrad, enm, cvm, s0out)
    !Pisanje koncne konfiguracije polaronov v sistemu v output
    open(10,file=""//trim(pwd)//"konf_latt&
    &_"//trim(lattice)//"_int_"//trim(intname)//"_J_"//trim(jayzi)//"_"//trim(jayzd)//"_C_&
    &"//trim(czi)//"_"//trim(czd)//"_radij_"//trim(str(radij))//"_mrad_"//trim(str(mrad))//"_ls&
    &_"//trim(lszi)//"_"//trim(lszd)//"_side1_"//trim(str(side1))//"_side2_"//trim(str(side2))//"_&
    &N_"//trim(str(nbar(i)))//"_ni_"//trim(strl(nit))//"_intv_"//trim(strl(intvt))//"_intt&
    &_"//trim(strl(inttt))//"_tempz_"//trim(tempzzi)//"_"//trim(tempzzd)//"_tempk_"//trim(tempkzi)//"&
    &_"//trim(tempkzd)//"_tempi_"//trim(tempizi)//"_"//trim(tempizd)//".dat")
    do j=1,side1
        write(10,*)(s0out(lt,j,k),k=1,side2)
    end do
    close(10)
    open(10,file=""//trim(pwd)//"all_konf_latt&
    &_"//trim(lattice)//"_int_"//trim(intname)//"_J_"//trim(jayzi)//"_"//trim(jayzd)//"_C_&
    &"//trim(czi)//"_"//trim(czd)//"_radij_"//trim(str(radij))//"_mrad_"//trim(str(mrad))//"_ls&
    &_"//trim(lszi)//"_"//trim(lszd)//"_side1_"//trim(str(side1))//"_side2_"//trim(str(side2))//"_&
    &N_"//trim(str(nbar(i)))//"_ni_"//trim(strl(nit))//"_intv_"//trim(strl(intvt))//"_intt&
    &_"//trim(strl(inttt))//"_tempz_"//trim(tempzzi)//"_"//trim(tempzzd)//"_tempk_"//trim(tempkzi)//"&
    &_"//trim(tempkzd)//"_tempi_"//trim(tempizi)//"_"//trim(tempizd)//".dat")
    do x=1,lt
      do j=1,side1
          write(10,*)(s0out(x,j,k),k=1,side2)
      end do
    end do
    close(10)
    open(10,file=""//trim(pwd)//"en_latt&
    &_"//trim(lattice)//"_int_"//trim(intname)//"_J_"//trim(jayzi)//"_"//trim(jayzd)//"_C_&
    &"//trim(czi)//"_"//trim(czd)//"_radij_"//trim(str(radij))//"_mrad_"//trim(str(mrad))//"_ls&
    &_"//trim(lszi)//"_"//trim(lszd)//"_side1_"//trim(str(side1))//"_side2_"//trim(str(side2))//"&
    &_N_"//trim(str(nbar(i)))//"_ni_"//trim(strl(nit))//"_intv_"//trim(strl(intvt))//"_intt&
    &_"//trim(strl(inttt))//"_tempz_"//trim(tempzzi)//"_"//trim(tempzzd)//"_tempk&
    &_"//trim(tempkzi)//"_"//trim(tempkzd)//"_tempi_"//trim(tempizi)//"_"//trim(tempizd)//".dat")
    do j=1,lt
        write (10,*) (enm(j,k),k=1,3)
    end do
    close (10)
    open(10,file=""//trim(pwd)//"cv_latt&
    &_"//trim(lattice)//"_int_"//trim(intname)//"_J_"//trim(jayzi)//"_"//trim(jayzd)//"_C&
    &_"//trim(czi)//"_"//trim(czd)//"_radij_"//trim(str(radij))//"_mrad_"//trim(str(mrad))//"_&
    &ls_"//trim(lszi)//"_"//trim(lszd)//"_side1_"//trim(str(side1))//"_side2_"//trim(str(side2))//"_&
    &N_"//trim(str(nbar(i)))//"_ni_"//trim(strl(nit))//"_intv_"//trim(strl(intvt))//"_intt&
    &_"//trim(strl(inttt))//"_tempz_"//trim(tempzzi)//"_"//trim(tempzzd)//"_tempk_"//trim(tempkzi)//"&
    &_"//trim(tempkzd)//"_tempi_"//trim(tempizi)//"_"//trim(tempizd)//".dat")
    do j=1,lt
        write (10,*) (cvm(j,k),k=1,3)
    end do
    close (10)
!Konec zanke stevila polaronov v sistemu
    end do
20 radij=0
call CPU_TIME(finish)
print *, finish-start
end program main

subroutine spini(lt, temp, v, side1, side2, range1, range2, ni, intv, intt, nbar, enke, s0, mrad, enm, cvm, s0out)
!Napoved spremenljivk
implicit none
integer(kind=8)               :: ni, i, j, k, lt, intv, intt, ns
integer                       :: o, p, x, side1, side2, xr, yr, nbar, count, range1, range2, mrad
integer                       :: enke(2*nbar), s0(side1,side2), s0out(lt,side1,side2)
double precision              :: v(2*range1+1,2*range2+1), temp(lt), enm(lt,3), cvm(lt,3)
double precision, allocatable :: ar(:)
double precision              :: u, deltaen, t, faktor, nbardi, en, cv, enerr, cverr
!Zacetek dinamike polaronov
ns=(ni-intt)/intv
allocate(ar(ns))
nbardi=1./dble(nbar)
do j=1,lt
  count=1
  t=temp(j)
  do i=1,ni
    if (i > intt .and. mod(i,intv)==0 .and. count<=ns) then
      call energija(nbar, range1, range2, side1, side2, enke, v, u)
      ar(count)=u*nbardi
      count=count+1
    end if
    !Zamenjava delca in njegove vrzeli v okviru local update scheme
    call RANDOM_NUMBER(u)
    x=1+2*floor(nbar*u)
    call RANDOM_NUMBER(u)
    xr=enke(x)-mrad+floor((2*mrad+1)*u)
    if (xr<=0) then
        xr=xr+side1
    else if (xr>side1) then
        xr=xr-side1
    else
        continue
    end if
    call RANDOM_NUMBER(u)
    yr=enke(x+1)-mrad+floor((2*mrad+1)*u)
    if (yr<=0) then
        yr=yr+side2
    else if (yr>side2) then
        yr=yr-side2
    else
        continue
    end if
    if (s0(xr,yr)==0) then
        call faktorfe(x, xr, yr, side1, side2, nbar, range1, range2, enke, v, faktor)
        deltaen=faktor
        call faktorfn(x, enke(x), enke(x+1), side1, side2, nbar, range1, range2, enke, v, faktor)
        deltaen=deltaen+faktor
        !Random odlocanje o sprejemu oziroma zavracanju zamenjave polaronov na podlagi energijske ugodnosti spremembe
        if (deltaen<=0.) then
            s0(xr,yr)=1
            s0(enke(x),enke(x+1))=0
            enke(x)=xr
            enke(x+1)=yr
        else
            call random_number(u)
            if (u<dexp(-deltaen/t)) then
                s0(xr,yr)=1
                s0(enke(x),enke(x+1))=0
                enke(x)=xr
                enke(x+1)=yr
            else
                continue
            end if
        end if
    else
        continue
    end if
  !Konec dinamike
  end do
  cv=0.
  cverr=0.
  do i=1,10000
    en=0.
    enerr=0.
    do k=1,ns
      call random_number(u)
      x=1+floor(ns*u)
      en=en+ar(x)
      enerr=enerr+ar(x)*ar(x)
    end do
    en=en/dble(ns)
    enerr=enerr/dble(ns)
    cv=cv+dble(nbar)*(enerr-en*en)/t/t
    cverr=cverr+nbar*(enerr-en*en)/t/t*nbar*(enerr-en*en)/t/t
  end do
  cv=cv/dble(10000)
  cverr=cverr/dble(10000)
  cverr=dsqrt(cverr-cv*cv)
  en=0.
  enerr=0.
  do i=1,ns
    en=en+ar(i)
    enerr=enerr+ar(i)*ar(i)
  end do
  en=en/dble(ns)
  enerr=enerr/dble(ns)
  cv=nbar*(enerr-en*en)/t/t
  enerr=dsqrt((enerr-en*en)/dble(ns-1))
  enm(j,1)=t
  enm(j,2)=en
  enm(j,3)=enerr
  cvm(j,1)=t
  cvm(j,2)=cv
  cvm(j,3)=cverr
  do o=1,side1
    do p=1,side2
      s0out(j,o,p)=s0(o,p)
    end do
  end do
  !Konec temperaturne zanke
end do
end subroutine spini
    
subroutine energija(nbar, range1, range2, side1, side2, enke, v, en)
implicit none
integer          :: i, j, nbar, range1, range2, side1, side2, di, dj
integer          :: enke(2*nbar)
double precision :: en
double precision :: v(2*range1+1,2*range2+1)
en=0.
do i=1,2*nbar,2
    do j=i+2,2*nbar,2
        di=enke(j)-enke(i)
        if (di>range1) then
            di=di-side1
        else if (di<-range1) then
            di=di+side1
        end if
        dj=enke(j+1)-enke(i+1)
        if (dj>range2) then
            dj=dj-side2
        else if (dj<-range2) then
            dj=dj+side2
        end if
        en=en+v(di+range1+1,dj+range2+1)
    end do
end do
    end subroutine
    
subroutine faktorfe(x, xr, yr, side1, side2, nbar, range1, range2, enke, v, faktor)
implicit none
integer          :: i, x, xr, yr, di, dj, side1, side2, nbar, range1, range2, enke(2*nbar)
double precision :: faktor, v(2*range1+1,2*range2+1)
faktor=0.
do i=1,x-2,2
    di=enke(i)-xr
    if (di>range1) then
        di=di-side1
    else if (di<-range1) then
        di=di+side1
    end if
    dj=enke(i+1)-yr
    if (dj>range2) then
        dj=dj-side2
    else if (dj<-range2) then
        dj=dj+side2
    end if
    faktor=faktor+v(di+range1+1,dj+range2+1)
end do
do i=x+2,2*nbar,2
    di=enke(i)-xr
    if (di>range1) then
        di=di-side1
    else if (di<-range1) then
        di=di+side1
    end if
    dj=enke(i+1)-yr
    if (dj>range2) then
        dj=dj-side2
    else if (dj<-range2) then
        dj=dj+side2
    end if
    faktor=faktor+v(di+range1+1,dj+range2+1)
end do
end subroutine
    
subroutine faktorfn(x, xk, yk, side1, side2, nbar, range1, range2, enke, v, faktor)
implicit none
integer          :: i, x, xk, yk, di, dj, side1, side2, nbar, range1, range2, enke(2*nbar)
double precision :: faktor, v(2*range1+1,2*range2+1)
faktor=0.
do i=1,x-2,2
    di=enke(i)-xk
    if (di>range1) then
        di=di-side1
    else if (di<-range1) then
        di=di+side1
    end if
    dj=enke(i+1)-yk
    if (dj>range2) then
        dj=dj-side2
    else if (dj<-range2) then
        dj=dj+side2
    end if
    faktor=faktor-v(di+range1+1,dj+range2+1)
end do
do i=x+2,2*nbar,2
    di=enke(i)-xk
    if (di>range1) then
        di=di-side1
    else if (di<-range1) then
        di=di+side1
    end if
    dj=enke(i+1)-yk
    if (dj>range2) then
        dj=dj-side2
    else if (dj<-range2) then
        dj=dj+side2
    end if
    faktor=faktor-v(di+range1+1,dj+range2+1)
end do
end subroutine

subroutine dec2str(decin, si, sd)
double precision :: decin
character(len=*) :: si, sd
character*30 :: s
character*30 :: strr
s=strr(decin)
call SplitString(s,si,sd,".")
call ReplaceText(sd,"0","",sd)
end subroutine

!Sprememba integerja v string
character(len=30) function str(k)
  integer, intent(in) :: k
  write (str,*) k
  str = adjustl(str)
end function str

!Sprememba integerja v string
character(len=30) function strl(k)
  integer(kind=8), intent(in) :: k
  write (strl,*) k
  strl = adjustl(strl)
end function strl

!Sprememba reala v string
character(len=30) function strr(k)
  integer, parameter  :: ikind=selected_real_kind(p=15)
  real (kind=ikind)   :: k
  write (strr, '(F16.8)') k
  strr = adjustl(strr)
end function strr

!Razclenitev stringa na dva dela, ki ju lucuje predpisan znak
subroutine SplitString(instring, string1, string2, delim)
implicit none
character(len=*)              :: instring,delim
character(len=*), intent(out) :: string1,string2
integer                       :: index
instring = trim(instring)
index = scan(instring,delim)
string1 = instring(1:index-1)
string2 = instring(index+1:)
end subroutine SplitString

!Zamenjava vseh zankov enakega predpisanega tipa za drug predpisan tip
subroutine ReplaceText (s,text,rep,outs)
!Napoved spremenljivk
implicit none
character(len=*) :: s,text,rep
character(len=*) :: outs
integer          :: i, nt, nr
!V primeru samih nicel v stringu zelimo imeti za outs samo 0
if (scan(outs,"123456789")==0) then
  outs="0"
else
  !Definiranje outputa in pregled dolzine brez presledkov stringa, ki moramo zamenjati in stringa, ki pride na njegovo mesto
  outs = s
  nt = len_trim(text)
  nr = len_trim(rep)
  !Zamenjava vseh stringov, ki ustrezajo vnosu
  do i=len(outs),1,-1
    if (scan(outs(i:len(outs)),"123456789")/=0) then
      continue
    else
      outs = outs(:i-1) // rep(:nr) // outs(i+nt:)
    end if
  end do
end if
end subroutine ReplaceText

