;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;get the 1D model results 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
f=addfile(pathin+"/"+filein,"r")

tavg=-1.
if (casename.eq."CPBL2".or.casename.eq."CPBL4") then
tavg=600.
t_ind=(/3.,5./)*3600.
end if
if (casename.eq."BOMEX") then
tavg=1800.
t_ind=(/2.5/)*3600.
end if

if (tavg.eq.-1.) then
  print("ERROR: Case not specified in NCL postproc!")
  exit
end if

z=(f->z)
zi=(f->z2)
time=f->time
nzm=dimsizes(z)

thl1D = new((/dimsizes(t_ind),nzm/),float)
thv1D = new((/dimsizes(t_ind),nzm/),float)
qt1D = new((/dimsizes(t_ind),nzm/),float)
thvd1D = new((/dimsizes(t_ind),nzm/),float)
thvd1D@_FillValue=-999.
thvd1D=0.0
qtd1D = new((/dimsizes(t_ind),nzm/),float)
qtd1D@_FillValue=-999.
qtd1D=0.0

wthl1D = new((/dimsizes(t_ind),nzm+1/),float)
wqt1D = new((/dimsizes(t_ind),nzm+1/),float)

if (casename.eq."BOMEX") then
qc1D = new((/1,nzm/),float)
;rhoi1D= new(nzm+1,float)
else
thgrad1D=new((/dimsizes(t_ind)/),float)
end if
plkmax1D=new((/dimsizes(t_ind)/),float)
plkmax1D=int2flt(nzm-2)

do t=0,dimsizes(t_ind)-1
    tmin_vec = ind(abs(time-t_ind(t)+tavg).eq.min(abs(time-t_ind(t)+tavg)))
    tmax_vec = ind(abs(time-t_ind(t)-tavg).eq.min(abs(time-t_ind(t)-tavg)))
    tmin=tmin_vec(0)
    tmax=tmax_vec(0)
    nt  =tmax-tmin+1
    avgind = ispan(tmin,tmax,1)
    thl1D(t,:)=(dim_avg_n(f->thl(tmin:tmax,:),0))
    thv1D(t,:)=(dim_avg_n(f->thr(tmin:tmax,:),0))
    qt1D(t,:)=(dim_avg_n(f->qt(tmin:tmax,:),0))

    upthv=f->upthv(tmin:tmax,:,:)
    upqt=f->upqt(tmin:tmax,:,:)
    fraca=f->upa(tmin:tmax,:,:)
    fractot=dim_sum_n(fraca,1)

    do k=0,nzm-1
    navg=0
    do indi=0,nt-1
      if (fractot(indi,k).ge.0.01) then
         fracsmall=0.
         fracsmallold=0.
         np2=39
         lfrac=True
         do while (np2.ge.0.and.lfrac)
          fracsmall=fracsmallold+fraca(indi,np2,k)
          if (fracsmall.gt.0.01) then
            wgt = (0.01-fracsmallold)/0.01
          else
            wgt = fraca(indi,np2,k)/0.01
          end if
          thvd1D(t,k) = thvd1D(t,k) + upthv(indi,np2,k) * wgt
          qtd1D(t,k) = qtd1D(t,k) + upqt(indi,np2,k) * wgt
          fracsmallold=fracsmall
          if (fracsmallold.gt.0.01) then
            lfrac=False
          end if
          np2=np2-1
         end do
         navg=navg+1
      end if
    end do
    if (navg.gt.0) then
       thvd1D(t,k) = thvd1D(t,k) / int2flt(navg)
       qtd1D(t,k) = qtd1D(t,k) / int2flt(navg)
    else
       thvd1D(t,k) = -999.
       qtd1D(t,k) = -999.
    end if
    end do
  

    do k=0,nzm-2
        if ((.not.ismissing(thvd1D(t,k))).and.(.not.ismissing(thvd1D(t,k+1)))) then
           tmp = 0.5*(thvd1D(t,k)+thvd1D(t,k+1))
           thvd1D(t,k)= tmp - thv1D(t,k)
           tmp = 0.5*(qtd1D(t,k)+qtd1D(t,k+1))
           qtd1D(t,k)= tmp - qt1D(t,k)
           delete(tmp)
        else
           thvd1D(t,k)= -999.
           qtd1D(t,k)= -999.
           if (plkmax1D(t).eq.int2flt(nzm-2)) then
              plkmax1D(t) = int2flt(k-1)
           end if
        end if
    end do
    thvd1D(t,nzm-1)= -999.
    qtd1D(t,nzm-1)= -999.
    if (plkmax1D(t).lt.0.) then
      print("Plume maximum level below lowest level")
      exit
    end if

  
    wthl1D(t,:)=(dim_avg_n(f->tflx(avgind,:),0) )
    wqt1D(t,:)=(dim_avg_n(f->qtflx(avgind,:),0)) 

    if (casename.eq."BOMEX") then
    qc1D(0,:)=(dim_avg_n(f->qcl(avgind,:),0))
    else
      tmp=thl1D(t,:)
      zlow=1200.+t*200.
      zhigh=zlow+400.
      zlow=700.
      tmp=where(z.gt.zhigh.or.z.lt.zlow,-999.,tmp)
      tmp@_FillValue=-999.
      thgrad1D(t)=regline(doubletofloat(z),tmp)
      delete(tmp)
    end if
    delete(avgind)
    delete(upthv)
    delete(upqt)
    delete(fraca)
    delete(fractot)
    delete(tmin_vec)
    delete(tmax_vec)
end do


system("/bin/rm -f "+ pathin +"/mminput.nc" )   ; remove any pre-existing file
ncdf = addfile(pathin+"/mminput.nc" ,"c")  ; open output netCDF file

fAtt               = True
fAtt@title         = "Input to metamodel"
fAtt@source_file   =  "single column model"
fAtt@Conventions   = "None"
fAtt@creation_date = systemfunc ("date")
fileattdef( ncdf, fAtt )            ; copy file attributes

filedimdef(ncdf,"time",dimsizes(t_ind),False)
filedimdef(ncdf,"z",nzm,False)
filedimdef(ncdf,"zi",nzm+1,False)



thl1D!0 = "time"
thl1D!1 = "z"
thl1D&time = t_ind
thl1D&z  = z
thl1D@long_name = "liquid potential temperature"           ; assign attributes
thl1D@_FillValue = -999.           ; assign attributes
ncdf->THL  = thl1D

qt1D!0 = "time"
qt1D!1 = "z"
qt1D&time = t_ind
qt1D&z  = z
qt1D@long_name = "total water mass fraction"           ; assign attributes
qt1D@_FillValue = -999.           ; assign attributes
ncdf->QT  = qt1D

qtd1D!0 = "time"
qtd1D!1 = "z"
qtd1D&time = t_ind
qtd1D&z  = z
qtd1D@long_name = "perturbation in total water mass fraction"           ; assign attributes
qtd1D@_FillValue = -999.           ; assign attributes
ncdf->QTD  = qtd1D

thvd1D!0 = "time"
thvd1D!1 = "z"
thvd1D&time = t_ind
thvd1D&z  = z
thvd1D@long_name = "perturbation in virtual potential temperature"           ; assign attributes
thvd1D@_FillValue = -999.           ; assign attributes
ncdf->THVD  = thvd1D

wthl1D!0 = "time"
wthl1D!1 = "zi"
wthl1D&time = t_ind
wthl1D&zi  = zi
wthl1D@long_name = "flux of liquid water potential temperature"           ; assign attributes
wthl1D@_FillValue = -999.           ; assign attributes
ncdf->WTHL  = wthl1D

wqt1D!0 = "time"
wqt1D!1 = "zi"
wqt1D&time = t_ind
wqt1D&zi  = zi
wqt1D@long_name = "flux of water"           ; assign attributes
wqt1D@_FillValue = -999.           ; assign attributes
ncdf->WQT  = wqt1D

plkmax1D!0 = "time"
plkmax1D&time = t_ind
plkmax1D@long_name = "index of highest level with plumes; 0=1st mass level"           ; assign attributes
plkmax1D@_FillValue = -999.           ; assign attributes
ncdf->PLKMAX  = plkmax1D


if (casename.eq."BOMEX") then
;add ac
qc1D!0 = "time"
qc1D!1 = "z"
qc1D&time = t_ind
qc1D&z  = z
qc1D@long_name = "liquid water mass fraction"           ; assign attributes
qc1D@_FillValue = -999.           ; assign attributes
ncdf->QC  = qc1D

else
;add gradient of thl

thgrad1D!0 = "time"
thgrad1D&time = t_ind
thgrad1D@long_name = "potential temperature gradient"           ; assign attributes
thgrad1D@_FillValue = -999.           ; assign attributes
ncdf->THGRAD  = thgrad1D


end if

