% program reanal_amp_vs_time
%
% makes probability density function plots 
%
% Steven Cavallo
% February 2, 2009

clear all; close all; beep off;


trackdirs = {%'/home/disk/pvort/nobackup/scavallo/polar_data/track_files/original_track_files/' ...
             '/home/disk/pvort/nobackup/scavallo/polar_data/track_files/wrf_track_files/' ...
             '/home/disk/pvort/nobackup/scavallo/polar_data/track_files/wrf_track_files/' ...
};

tracknames = {%'jja_bins/track.ncep.C.jja.1990.1999.2days.60N65.dat' ...
              'djf_bins/track.wrf.C.djf.1990.1999.2days.60N65.dat' ...
              'djf_bins_norad/track.wrf.C.djf.1990.1999.2days.60N65.dat' ...
};

%legendtext = { 'NNRP' 'WRF-NNRP'};
legendtext = {'Full physics' 'No radiation'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savedir = '/home/disk/pvort/scavallo/model_default/images/';
label_fontsize = 22;
title_fontsize = 14;
ylimits = [7 16];
nhrs_min = 24*5;
smooth_plot = true;
num_passes = 5;
saveprefix = 'wrf_djf_vs_norad'
%saveprefix = 'wrf_vs_nnrp_jja';
image_format = 'eps'; % Format to save image
image_res = '-r300'; % Image resolution
cols = {'k' 'k--' 'k.-' 'k+'};
pimg = (['-d' image_format]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End user options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nhrs = 24*10;   
ntimes = nhrs / 6;
nsteps_min = nhrs_min / 6;
for ff=1:length(tracknames)
   ftopen = strcat([char(trackdirs(ff)) char(tracknames(ff))]);
   [tlat tlon tdate ttheta ttamp tvrad tpamp] = textread(ftopen,'%9.2f%8.2f%12.0f%10.2f%8.2f%8.2f%8.2f');  
   [inds] = find(tdate == 0 & ttheta > nsteps_min);      
   length(inds)
   maxtracks = max(ttheta(inds));
   amps(1:length(inds),1:maxtracks) = NaN;
   lats(1:length(inds),1:maxtracks) = NaN;
   for ii=1:length(inds)
      indnow = inds(ii);
      ntracksnow = ttheta(indnow);
      timenow = 1;
      for jj=1:ntracksnow
         amps(ii,jj) = ttamp(indnow+jj);
	 lats(ii,jj) = tlat(indnow+jj);
	 rads(ii,jj) = tvrad(indnow+jj);
      end
   end
 
   % Average for over each time step
   amps_avg(1:size(amps,2)) = NaN;
   lats_avg(1:size(amps,2)) = NaN;
   rads_avg(1:size(amps,2)) = NaN;
   for ii=1:size(amps,2)
      ampsnow = squeeze(amps(:,ii));
      tmp = find(~isnan(ampsnow));      
      amps_avg(ii) = mean(ampsnow(tmp));            
      clear tmp;
      
      latsnow = squeeze(lats(:,ii));
      tmp = find(~isnan(latsnow));      
      lats_avg(ii) = mean(latsnow(tmp)); 
      clear tmp;
      
      radsnow = squeeze(rads(:,ii));
      tmp = find(~isnan(radsnow));      
      rads_avg(ii) = mean(radsnow(tmp)); 
      clear tmp;
   end
 
   if ff == 1
      amps_control = amps_avg;      
   elseif ff == 2
      if length(amps_control) == length(amps_avg)
         %percent_change = ((amps_avg - amps_control') ./ amps_control').*100;
      else
         if length(amps_control) < length(amps_avg)
	    eind = length(amps_control);	    
	 else
	    eind = length(amps_avg);
	 end
	 %percent_change = ((amps_avg(1:eind) - amps_control(1:eind)') ./ amps_control(1:eind)').*100;
      end
   end
   
   
   hrange = 1:1:length(amps_avg);   
   hrs = (hrange-1).*6;
   [eind] = find(hrs == nhrs_min);
   
   %aa = polyfit(hrs,amps_avg,7);
   %bb = polyval(aa,hrs);
   if smooth_plot == true
      for ii=1:num_passes
	 amps_avg = smooth(amps_avg);
	 lats_avg = smooth(lats_avg);
	 rads_avg = smooth(rads_avg);
	 if ff == 2
	    %percent_change = smooth(percent_change);
	 end
      end
   end
   
   days = hrs./24;
   
   figure(1); 
   if ff == 1; clf; end;
   h = plot(days(1:eind),amps_avg(1:eind),char(cols(ff)));  hold on;
   set(h,'linewidth',3);
   xlabel('Vortex age (days)','Fontname','Times','Fontweight','bold','Fontsize',label_fontsize);
   ylabel('Amplitude','Fontname','Times','Fontweight','bold','Fontsize',label_fontsize);
   set(gca,'box','on','linewidth',3,'Fontname','Times','Fontsize',label_fontsize,'Fontweight','bold');
   if ff == length(tracknames); 
     axis tight;
     %ylims = get(gca,'Ylim'); set(gca,'Ylim',[ylims(1)-1 ylims(2)+1]);
     set(gca,'Ylim',[ylimits(1) ylimits(2)]);
     
     %L = legend(char(legendtext(1)),char(legendtext(2)),4); 
     %set(L,'Fontname','Times','Fontweight','bold','Fontsize',label_fontsize);
     %grid on;
     pfile = [savedir saveprefix '_amps_vs_time.' image_format];
     print(pimg,image_res,pfile);
   end   
   
   if 1 == 2
      figure(2); clf;
      h = plot(days(1:eind),percent_change(1:eind),char(cols(ff)));  hold on;
      set(h,'linewidth',3);
      xlabel('Vortex age (days)','Fontname','Times','Fontweight','bold','Fontsize',label_fontsize);
      ylabel('Amplitude percent change','Fontweight','bold','Fontsize',label_fontsize);
      set(gca,'box','on','linewidth',3,'Fontname','Times','Fontsize',label_fontsize,'Fontweight','bold');      
      axis tight;
      ylims = get(gca,'Ylim');
      if ( ylims(1)<0 & ylims(2)>0); set(gca,'Ylim',[ylims(1) -ylims(1)]); end; 	
      %grid on;      
      pfile = [savedir saveprefix '_amps_vs_time_percent_change.' image_format];
      print(pimg,image_res,pfile);
   end

   figure(3); 
   if ff == 1; clf; end;
   h = plot(days(1:eind),lats_avg(1:eind),char(cols(ff)));  hold on;
   set(h,'linewidth',3);
   xlabel('Vortex age (days)','Fontname','Times','Fontweight','bold','Fontsize',label_fontsize);
   ylabeltext = ['Latitude (^{\circ}N)'];
   ylabel(ylabeltext,'Fontname','Times','Fontweight','bold','Fontsize',label_fontsize);
   set(gca,'box','on','linewidth',3,'Fontname','Times','Fontsize',label_fontsize,'Fontweight','bold');
   if ff == length(tracknames); 
      axis tight;
      set(gca,'Ylim',[67 75]);
      %ylims = get(gca,'Ylim'); set(gca,'Ylim',[ylims(1)-1 ylims(2)+1]);
      %L = legend(char(legendtext(1)),char(legendtext(2)),4); 
      %set(L,'Fontweight','bold','Fontsize',label_fontsize);
      %grid on;
      pfile = [savedir saveprefix '_lats_vs_time.' image_format];
      print(pimg,image_res,pfile);
   end   

   figure(4); 
   if ff == 1; clf; end;
   h = plot(hrs(1:eind),rads_avg(1:eind),char(cols(ff)));  hold on;
   set(h,'linewidth',3);
   xlabel('Vortex age (hours)','Fontname','Times','Fontweight','bold','Fontsize',label_fontsize);
   ylabeltext = ['Radius (km)'];
   ylabel(ylabeltext,'Fontname','Times','Fontweight','bold','Fontsize',label_fontsize);
   set(gca,'box','on','linewidth',3,'Fontname','Times','Fontsize',label_fontsize,'Fontweight','bold');
   if ff == length(tracknames); 
      axis tight;
      %L = legend(char(legendtext(1)),char(legendtext(2)),4); 
      %set(L,'Fontweight','bold','Fontsize',label_fontsize);
      %grid on;
      pfile = [savedir saveprefix '_rads_vs_time.' image_format];
      print(pimg,image_res,pfile);
   end   
   
   clear amps lats;
end

