pmin=(/0.,0.,0.,0.,100.,0./)
pmax=(/4.,4.,2000.,2.,1000.,2./)


colors=(/"red","blue","blue","red"/)

best01=0.473029

is=0
ie=71
ofs=asciiread("./opti_new/of_srbf.txt",(/300,20/),"float")
ofd=asciiread("./opti_new/of_dy.txt",(/300,20/),"float")

maske=new((/2,300,20/),float)
maske=-999.
maske(0,is:ie,:)=where(ofs(is:ie,:).lt.best01,1.,-999.)
maske(1,is:ie,:)=where(ofd(is:ie,:).lt.best01,1.,-999.)
maske@_FillValue=-999.


nums=0
numd=0
do i=0,19
  if (sum(maske(0,is:ie,:)).ge.1) then
     nums=nums+1
  end if
  if (sum(maske(1,is:ie,:)).ge.1) then
     numd=numd+1
  end if
end do

print(nums)
print(numd)


filebase=(/"x_srbf","x_dycors"/)

npara=6
data=new((/dimsizes(filebase),npara,300,20/),float)
data@_FillValue=-999.
do j=0,dimsizes(filebase)-1
do i=0,npara-1
  filename="./opti_new/"+filebase(j) +"_"+i+".txt"
  tmp= asciiread(filename,(/300,20/),"float")
  data(j,i,:,:) = pmin(i) + maske(j,:,:) * tmp * (pmax(i)-pmin(i))
end do
end do
print(num(ismissing(data(:,0,is:ie,:))))
print(num(.not.ismissing(data(0,0,is:ie,:))))
print(num(.not.ismissing(data(1,0,is:ie,:))))
print((ie-is+1)*2*20)

rmsedir="/global/scratch/langhans/metamodel_large2/"
f=addfile(rmsedir+"/RMSE_mmthvdqtdrmse.nc","r")


p1=f->alphathv
p2=f->alphaqt
p3=f->ctketau
p4=f->pwmin
p5=f->tauneggers
p6=f->Wc

pnames=(/"a!BF8!q!F!v!N!","a!B!qt!N!","!F8!t!FB!l!NF!","w!B!min!N!","!F8!t!F!!BF8!e!NF!","b!BF8!e!NF!"/)

p=new((/6,dimsizes(p1)/),float)
do i=0,5
  p(i,:) = fspan(pmin(i),pmax(i),dimsizes(p1))
end do


pdf = new((/dimsizes(filebase),npara,(dimsizes(p1)-1)/),float)
pdfcenters = new((/6,(dimsizes(p1)-1)/),float)
optpdf=True

do j=0,dimsizes(filebase)-1
do i=0,npara-1
  optpdf@bin_min=pmin(i)
  optpdf@bin_max=pmax(i)
  zpdf = pdfx(data(j,i,:,:),(dimsizes(p1)-1),optpdf)
  pdfcenters(i,:)=doubletofloat(zpdf@bin_center)
  zpdf = zpdf / (pdfcenters(i,1)-pdfcenters(i,0))
  pdf(j,i,:) = doubletofloat(zpdf)
  delete(zpdf)
end do
end do

pdf=pdf/100.


plotpdf=new(6,graphic)

  wks = gsn_open_wks("eps","pdf_stochrbf")               ; send graphics to PNG file
;**********************************************
  res=True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tiYAxisString=" "
  res@xyLineThicknessF=3.
;  ;res@gsnYRefLine=(/best1,best5/)
  res@xyMarkLineMode="MarkLines"
  res@xyLineColors= colors
  res@xyMonoDashPattern = False
  res@xyDashPatterns=(/0,0,2,2/)
  res@xyMarkerColors=colors
  res@xyMarker = 1
  res@xyMarkerSizeF=0.02
  res@tiXAxisFontHeightF=0.03
  res@trYMinF = 0.
  res@gsnXRefLineColor="grey40"
  res@gsnXRefLineThicknessF=4.

  opt1=(/0.2,0.6,0.0,0.21,449.,0.55/)
  opt2=(/0.6,0.88,1187.,0.72,407.,0.48/)

  ymaxvec=(/2.,1.2,0.002,2.,0.0085,4.0/)
  do i=0,5
;  res@gsnXRefLine=(/opt1(i),opt2(i)/)
  res@tiXAxisString=pnames(i)
  res@tiYAxisString="f("+pnames(i)+")"+"/!F8!D!F!"+pnames(i)
  res@tmXBLabels =  p(i,:) ; labels for each box
  res@trXMinF =p(i,0)
  res@trXMaxF = p(i,dimsizes(p1)-1)
  res@trYMaxF = ymaxvec(i)  ;1./(2.*pdfcenters(i,0))

  plotpdf(i)  = gsn_csm_xy (wks,pdfcenters(i,:),pdf(0,i,:),res)
  end do

  resp=True
  resp@gsnPanelXWhiteSpacePercent=5.

  gsn_panel(wks,plotpdf,(/2,3/),resp)









