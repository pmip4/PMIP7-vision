;This script cycles through all the files in a "data" directory and calculates the standard devations of the timeseries and their correlations/regressions

;History:
; 1. It was originally used to create the supplementary table for Brierley & Wainer (2018) - although that was then tidied up in Excel to make it presentable
; 2. It was then modified to output amplitude of EOF-based modes, as well as index-based ones. The output of this script was then plotted as Fig 3 in Rehfeld et al (2020,Earth Sys. Dynam.)
; 3. Now it's being overwritten to provide some tables as inputs for my PMIP4_past2future_analyzer repo

; run with:
; conda activate ncl-nco; ncl -n compile_tidy_numbers.ncl

load "cvdp_data.functions.ncl"

data_dir="cvdp_output_by_experiment/"
;only the main per-simulation file, whose year range follows "cvdp_data." directly.
;CVDP v6.1 also writes auxiliary per-variable files (.pdv., .siconc., .amv.tas. etc)
;that hold no timeseries we want and break the filename parsing below.
files=systemfunc("ls "+data_dir+"{1pctCO2,abrupt4xCO2,lig127k,lgm,midHolocene,midPliocene-eoi400,piControl}/*cvdp_data.[0-9]*-[0-9]*.nc")
files = str_sub_str(files, data_dir, "")
;files=(/"INM-CM4-8_historical.cvdp_data.1850-2014.nc"/) ;select a single file for testing...

;define the filename for the outputting
fname="tidy_numbers_cvdp6.1.csv"
system("rm -f "+fname)


; loop over the various CVDP files
do fi = 0, dimsizes(files)-1
  print("Working on "+files(fi))
  ;set up what to look over
  variables=new(50, string)
  headers=new(50, string)
  headers(0)="model"
  headers(1)="experiment"
  headers(2)="years"
  file_parts = str_split(files(fi),"/")
  expt=file_parts(0)
  filename=file_parts(1)
  name_parts = str_split(filename,".")
  name_strs0 = str_split(name_parts(0),"_")
  name_strs2 = str_split(name_parts(2),"-")
  variables(0)=name_strs0(0) ;model
  if dimsizes(name_strs0).eq.1 then
    variables(1)="observations" ;experiment
  else
    variables(1)=name_strs0(1) ;experiment
  end if
  variables(2)=toint(name_strs2(1))-toint(name_strs2(0))+1 ;years
  delete([/file_parts,expt,filename,name_parts,name_strs0,name_strs2/])

  ts_opt=True
  ts_opt@NO_EXIT_ON_ERROR=True
  if any(variables(1).eq.(/"1pctCO2","abrupt4xCO2","rcp26","rcp85"/)) then
    ts_opt@subset=True
    ts_opt@subset_syear=-29
    ts_opt@subset_eyear=0
    ts_opt@DETREND=True
  end if

  ;Global Mean TAS
  headers(3)="GMTemp"
  fil=addfile(data_dir+files(fi),"r")
  if .not.isfilevar(fil,"tas_spatialmean_ann") then
    variables(3)=variables@_FillValue
  else
    variables(3)=extract_latlon_areastat(data_dir+files(fi),"tas_spatialmean_ann",(/-90.,90./),(/0.,360./),"mean")
  end if
  delete(fil)

  ;IOD mean gradient
  ts_opt_abs=ts_opt
  ts_opt_abs@make_absolute=True
  ts_opt_abs@DETREND=False
  headers(4)="ATL3 mean"
  variables(4)=stat_ts_var(data_dir+files(fi),"atl3_timeseries_mon","mean",ts_opt_abs)

  ;IOD Amplitude
  headers(5)="IOD stddev"
  variables(5)=stat_ts_var(data_dir+files(fi),"indian_ocean_dipole","stddev",ts_opt)

  ;Nino3.4 Amplitude
  headers(6)="Nino34 stddev"
  variables(6)=stat_ts_var(data_dir+files(fi),"nino34","stddev",ts_opt)

  ;atlantic_nino mean. Read from the *_timeseries_mon copy of the index: it is
  ;bit-identical but, unlike the bare index, still carries the AnnCycle attribute
  ;that make_absolute/AnnCycAmp need under CVDP v6.1.
  headers(7)="AtlNino mean"
  variables(7)=stat_ts_var(data_dir+files(fi),"atl3_timeseries_mon","mean",ts_opt_abs)

  ;atlantic_nino Mode
  headers(8)="AtlNino stddev"
  variables(8)=stat_ts_var(data_dir+files(fi),"atlantic_nino","stddev",ts_opt)

  ;atlantic_nino Annual Cycle
  headers(9)="AtlNino AnnCycAmp"
  variables(9)=stat_ts_var(data_dir+files(fi),"atl3_timeseries_mon","AnnCycAmp",ts_opt)

  ;Atlantic Meridional Mode
  headers(10)="AMM mean"
  variables(10)=stat_ts_var(data_dir+files(fi),"amm_timeseries_mon","mean",ts_opt_abs)

  ;Atlantic Meridional Mode
  headers(11)="AMM stddev"
  variables(11)=stat_ts_var(data_dir+files(fi),"atlantic_meridional_mode","stddev",ts_opt)

  ;AMM Annual Cycle
  headers(12)="AMM AnnCycAmp"
  variables(12)=stat_ts_var(data_dir+files(fi),"amm_timeseries_mon","AnnCycAmp",ts_opt)

  patterns=(/"tas_spatialmean_ann","pr_spatialmean_ann","nao_pattern_djf","pdv_pattern_mon"/)
  patterns_regions=(/(/(/-90.,90./),(/0.,359./)/),\ ;tas_spatialmean_ann - global mean 
                   (/(/-90.,90./),(/0.,359./)/),\ ;pr_spatialmean_ann - global mean
                   (/(/20.,80./),(/-90.,40./)/),\ ;nao_pattern_djf - NAO summer land NH 
                   (/(/20.,70./),(/110.,260./)/)/) ;pdo_pattern_djf - PDO
  patterns_stat=(/"mean","mean","stddev","stddev"/)
  patterns_human_names=(/"global_mean_absolute_temperatures","global_mean_precipitation","NAO DJF","PDO mon"/)
 
  do patt_i=0,dimsizes(patterns)-1
    f=addfile(data_dir+files(fi),"r")
    headers(13+patt_i)=patterns_human_names(patt_i)
    if isfilevar(f,patterns(patt_i)) then
      variables(13+patt_i)=extract_latlon_areastat(data_dir+files(fi),patterns(patt_i),patterns_regions(patt_i,0,:),patterns_regions(patt_i,1,:),patterns_stat(patt_i))
    else
      variables(13+patt_i)=variables@_FillValue
    end if
    delete(f)
  end do  
 
  ;write out to file
  exists=ind(.not.ismissing(headers))
  if fi.eq.0 then
    hlist=[/str_join(headers(exists),",")/]
    write_table(fname, "w", hlist, "%s")
  end if 
  vlist=[/str_join(variables(exists),",")/]
  write_table(fname, "a", vlist, "%s")

  delete([/headers,variables,ts_opt,ts_opt_abs/])
end do
