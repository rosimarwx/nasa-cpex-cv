% This script use MERRA2 specific humidity and vertical wind shear to make following plots:
% 1) maps of long-term means
% 2) cumulative probability distributions at specified locations and
% identify where CPEX-CV data fall in the distributions


%% set up parameters
% directory to save data
datadir = '/new_data/nsakaeda/CPEXCV/';

% directory to save figures
figdir = '/Users/nsakaeda/Figure/CPEXCV';

% MERRA2 data directory
datadir_merra2 = '/data/deluge/reanalysis/REANALYSIS/MERRA2/3D_PL/native/';

% range of years to collect MERRA data
yymin=2011; yymax=2021;

% month to collect MERRA data
mon_list = 9; mon_txt='Sep';

% two pressure levels to calculate shear
p_shear = [650 950];

% pressure level for specific humidity
p_qv = 925;

% lat/lon range to create long-term mean maps
lonbox_track = [-50 0];
latbox_track = [0 20];

% points to find distribution
lat_pts = [14.07 12.09];
lon_pts = [-21.78 -17.23];
% CPEX-CV vertical wind shear values at above points
shear_cpexcvobs = [4.17215208 21.18524376]; % m/s
q_cpexcvobs = [0.01507697 0.01553765]*10^3; % g/kg


%% go through MERRA2 files to calculate long-term mean and collect values at specified locations
% define the file name and location to save data
savefile = [datadir,'MERRA2comp_map_climo_',mon_txt,num2str(yymin),'-',num2str(yymax),'.mat'];
if exist(savefile,'file')==0

    % pressure levels to calculate
    p0_list = union(p_shear,p_qv);

    % load MERRA2 lat lon
    hour_merra2 = 0:3:21;
    merra2_file = [datadir_merra2,'2022/202209/MERRA2_400.inst3_3d_asm_Np.20220922.nc4'];
    lon_merra2 = ncread(merra2_file,'lon'); % -180 to +180
    lat_merra2 = ncread(merra2_file,'lat');
    lev_merra2 = ncread(merra2_file,'lev');
    lonidx_merra2 = find(lon_merra2>=min(lonbox_track)&lon_merra2<=max(lonbox_track));
    lonidx_merra2 = [lonidx_merra2(1)+(-2:-1)';lonidx_merra2;lonidx_merra2(end)+(1:2)'];
    lonidx_merra2 = lonidx_merra2(lonidx_merra2>=1&lonidx_merra2<=length(lon_merra2));
    latidx_merra2 = find(lat_merra2>=min(latbox_track)&lat_merra2<=max(latbox_track));
    latidx_merra2 = [latidx_merra2(1)+(-2:-1)';latidx_merra2;latidx_merra2(end)+(1:2)'];
    latidx_merra2 = latidx_merra2(latidx_merra2>=1&latidx_merra2<=length(lat_merra2));
    levidx_merra2 = find(ismember(lev_merra2,p0_list));
    
    % identify grids closest to the specified points
    latidx_pts = NaN(length(lat_pts),1);
    lonidx_pts = NaN(length(lon_pts),1);
    for n=1:length(lat_pts)
        latidx_pts(n) = find(abs(lat_merra2(latidx_merra2)-lat_pts(n))==min(abs(lat_merra2(latidx_merra2)-lat_pts(n))),1);
        lonidx_pts(n) = find(abs(lon_merra2(lonidx_merra2)-lon_pts(n))==min(abs(lon_merra2(lonidx_merra2)-lon_pts(n))),1);
    end
    grididx_pts = sub2ind([length(latidx_merra2),length(lonidx_merra2)],latidx_pts,lonidx_pts);
    lat_pts_merra2 = lat_merra2(latidx_merra2(latidx_pts));
    lon_pts_merra2 = lon_merra2(lonidx_merra2(lonidx_pts));

    % list of variables at the specified points
    u_list = [];
    v_list = [];
    q_list = [];

    % data points
    ndata_uv_climo = zeros(length(hour_merra2),length(levidx_merra2),length(latidx_merra2),length(lonidx_merra2));  
    ndata_q_climo = zeros(length(hour_merra2),length(levidx_merra2),length(latidx_merra2),length(lonidx_merra2));  
    % zonal wind
    usum_climo = zeros(length(hour_merra2),length(levidx_merra2),length(latidx_merra2),length(lonidx_merra2));  
    % meridional wind
    vsum_climo = zeros(length(hour_merra2),length(levidx_merra2),length(latidx_merra2),length(lonidx_merra2));  
    % water vapor
    qsum_climo = zeros(length(hour_merra2),length(levidx_merra2),length(latidx_merra2),length(lonidx_merra2)); 
    
    for yy=yymin:yymax
        for mm=mon_list
            % list of MERRA2 files for this month
            merra2_filelist = dir([datadir_merra2,sprintf('%04d/%04d%02d/',yy,yy,mm),'MERRA2_*.nc4']);
            merra2_filelist = cellfun(@(x) [datadir_merra2,sprintf('%04d/%04d%02d/',yy,yy,mm),x],{merra2_filelist(:).name},'UniformOutput',false);
            % load data
            for f=1:length(merra2_filelist)
                % date
                disp(merra2_filelist{f});
                %daysince = datenum(merra2_filelist{f}(min(strfind(merra2_filelist{f},'.nc'))+(-8:-1)),'yyyymmdd');
                daysince = ncreadatt(merra2_filelist{f},'time','units');
                daysince = datenum(daysince(length('minutes since')+2:end),'yyyy-mm-dd HH:MM:SS');
                % time
                time0 = double(ncread(merra2_filelist{f},'time'))./(24*60) + daysince;
                [~,~,~,HH0,~,~]=datevec(time0);
                % horizontal winds
                u0 = ncread(merra2_filelist{f},'U',[lonidx_merra2(1),latidx_merra2(1),1,1],[length(lonidx_merra2),length(latidx_merra2),Inf,Inf],[1,1,1,1]);
                u0 = u0(:,:,levidx_merra2,:);
                v0 = ncread(merra2_filelist{f},'V',[lonidx_merra2(1),latidx_merra2(1),1,1],[length(lonidx_merra2),length(latidx_merra2),Inf,Inf],[1,1,1,1]);
                v0 = v0(:,:,levidx_merra2,:);
                u0 = permute(u0,[4,3,2,1]); % time,lev,lat,lon
                v0 = permute(v0,[4,3,2,1]); % time,lev,lat,lon
                % water vapor
                q0 = ncread(merra2_filelist{f},'QV',[lonidx_merra2(1),latidx_merra2(1),1,1],[length(lonidx_merra2),length(latidx_merra2),Inf,Inf],[1,1,1,1]);
                q0 = q0(:,:,levidx_merra2,:);
                q0 = permute(q0,[4,3,2,1]); % time,lev,lat,lon

                for hh=1:length(hour_merra2)
                    hhidx = (HH0==hour_merra2(hh));
                    % climo
                    ndata_uv_climo(hh,:) = ndata_uv_climo(hh,:) + sum(~isnan(u0(hhidx,:)),1,'omitmissing');
                    ndata_q_climo(hh,:) = ndata_q_climo(hh,:) + sum(~isnan(q0(hhidx,:)),1,'omitmissing');
                    usum_climo(hh,:) = usum_climo(hh,:) + sum(u0(hhidx,:),1,'omitmissing');
                    vsum_climo(hh,:) = vsum_climo(hh,:) + sum(v0(hhidx,:),1,'omitmissing');
                    qsum_climo(hh,:) = qsum_climo(hh,:) + sum(q0(hhidx,:),1,'omitmissing');
                end

                % save variables at the specified points
                u_list = cat(1,u_list,u0(:,:,grididx_pts));
                v_list = cat(1,v_list,v0(:,:,grididx_pts));
                q_list = cat(1,q_list,q0(:,:,grididx_pts));

                clear u0 v0 q0 time0 HH0 daysince
            end
        end
        disp([num2str(yy),'done...']);
    end

    lat_merra2 = lat_merra2(latidx_merra2);
    lon_merra2 = lon_merra2(lonidx_merra2);
    lev_merra2 = lev_merra2(levidx_merra2);
    save(savefile,'ndata_*','usum_*','vsum_*','qsum_*','lat_merra2','lon_merra2','lev_merra2','hour_merra2',...
        'u_list','v_list','q_list','lat_pts_merra2','lon_pts_merra2','-v7.3');
    
end
load(savefile);

% find long-term mean
umean_climo = squeeze(sum(usum_climo,1)./sum(ndata_uv_climo,1));
vmean_climo = squeeze(sum(vsum_climo,1)./sum(ndata_uv_climo,1));
qmean_climo = squeeze(sum(qsum_climo,1)./sum(ndata_q_climo,1));

% grid closest to the specified points
latidx_pts = NaN(length(lat_pts),1);
lonidx_pts = NaN(length(lon_pts),1);
for n=1:length(lat_pts)
    latidx_pts(n) = find(abs(lat_merra2-lat_pts(n))==min(abs(lat_merra2-lat_pts(n))),1);
    lonidx_pts(n) = find(abs(lon_merra2-lon_pts(n))==min(abs(lon_merra2-lon_pts(n))),1);
end
grididx_pts = sub2ind([length(lat_merra2),length(lon_merra2)],latidx_pts,lonidx_pts);

% shear levels
ilev_shear = [find(lev_merra2==min(p_shear)) find(lev_merra2==max(p_shear))];
% qv levels
ilev_qv = find(lev_merra2==p_qv);



%% plot maps of long-term mean
% map of shear
ushear_climo = squeeze(umean_climo(ilev_shear(1),:,:)-umean_climo(ilev_shear(2),:,:));
vshear_climo = squeeze(vmean_climo(ilev_shear(1),:,:)-vmean_climo(ilev_shear(2),:,:));
shearmag_climo = (ushear_climo.^2 + vshear_climo.^2).^0.5;
% plotting variables
plot_u_climo = ushear_climo;
plot_v_climo = vshear_climo;
plot_fill_climo = shearmag_climo;
% colorbar
ticks_fill = 0:0.5:15;
grad_fill = flipud(pink(length(ticks_fill)-1));
cticks_fill = 0:3:15;
units_fill = 'm s^{-1}';
plot_fill_climo(plot_fill_climo<ticks_fill(1))=ticks_fill(1);
plot_fill_climo(plot_fill_climo>ticks_fill(end))=ticks_fill(end);
% wind scale
vint = 3; % plotting interval of merra2 wind vectors
vscale = 10; vdegscale = 1; % specific which speed in m/s should correspond to lat/lon deg
vlegend = 10;
lon_vlegend = lonbox(1)+[0 3]+0.1;
lat_vlegend = latbox(1)+[0 1.5]+0.1;
figure('Position',[1 1 600 425]); hold on
contourf(lon_merra2,lat_merra2,squeeze(plot_fill_climo),ticks_fill,'LineStyle','none');
colormap(gca,grad_fill); clim(ticks_fill([1 end]));
contour(lonTOPO,latTOPO,topo,[0.5 0.5],'-','LineWidth',1,'Color',[0.3 0.3 0.3]);
set(gca,'FontSize',16,'XLim',lonbox,'YLim',latbox,'XTick',[],'YTick',[],'Color',[0.7 0.7 0.7]);
text(lontick,min(latbox)*ones(size(lontick)),lonticklbl,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',14);
text(min(lonbox)*ones(size(lattick)),lattick,latticklbl,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14);
cbar=colorbar('YTick',cticks_fill,'YLim',ticks_fill([1 end]));
cbar.Label.String = units_fill;
ax_position = get(gca,'Position');
% add wind vectors
[X,Y] = meshgrid(lon_merra2(1:vint:end),lat_merra2(1:vint:end));
plot_u=plot_u_climo(1:vint:end,1:vint:end); plot_v=plot_v_climo(1:vint:end,1:vint:end);
X=X(:); Y=Y(:); plot_u=plot_u(:); plot_v=plot_v(:);
x_vlegend = interp1(lonbox,ax_position(1)+[0 ax_position(3)],lon_vlegend,'linear');
y_vlegend = interp1(latbox,ax_position(2)+[0 ax_position(4)],lat_vlegend,'linear');
for n=1:length(X)
    ltype = 'arrow';
    if X(n)<min(lonbox) || X(n)>max(lonbox) || Y(n)<min(latbox) || Y(n)>max(latbox)
        continue;
    end
    x_begin = interp1(lonbox,ax_position(1)+[0 ax_position(3)],X(n),'linear');
    x_end = interp1(lonbox,ax_position(1)+[0 ax_position(3)],X(n)+plot_u(n)*vdegscale/vscale,'linear','extrap');
    y_begin = interp1(latbox,ax_position(2)+[0 ax_position(4)],Y(n),'linear');
    y_end = interp1(latbox,ax_position(2)+[0 ax_position(4)],Y(n)+plot_v(n)*vdegscale/vscale,'linear','extrap');
    if y_end>(ax_position(2)+ax_position(4)), y_end = ax_position(2)+ax_position(4); ltype = 'line'; end
    if y_end<ax_position(2), y_end = ax_position(2); ltype = 'line'; end
    if x_end>(ax_position(1)+ax_position(3)), x_end = ax_position(1)+ax_position(3); ltype = 'line'; end
    if x_end<ax_position(1), x_end = ax_position(1); ltype = 'line'; end
    if ((x_begin>=min(x_vlegend)&&x_begin<=max(x_vlegend)) && (x_end>=min(x_vlegend)&&x_end<=max(x_vlegend))) && ...
            ((y_begin>=min(y_vlegend)&&y_begin<=max(y_vlegend)) && (y_end>=min(y_vlegend)&&y_end<=max(y_vlegend)))
        continue;
    end
    if (x_end>=min(x_vlegend)&&x_end<=max(x_vlegend)) && (y_end>=min(y_vlegend)&&y_end<=max(y_vlegend))
        ltype = 'line';
        if x_begin<min(x_vlegend)
            y_end = interp1([x_begin x_end],[y_begin y_end],min(x_vlegend));
            x_end = min(x_vlegend);
        elseif x_begin>max(x_vlegend)
            y_end = interp1([x_begin x_end],[y_begin y_end],max(x_vlegend));
            x_end = max(x_vlegend);
        elseif y_begin<min(y_vlegend)
            x_end = interp1([y_begin y_end],[x_begin x_end],min(y_vlegend));
            y_end = min(y_vlegend);
        elseif y_begin>max(y_vlegend)
            x_end = interp1([y_begin y_end],[x_begin x_end],max(y_vlegend));
            y_end = max(y_vlegend);
        end
    end
    if (x_begin>=min(x_vlegend)&&x_begin<=max(x_vlegend)) && (y_begin>=min(y_vlegend)&&y_begin<=max(y_vlegend))
        if x_end<min(x_vlegend)
            y_begin = interp1([x_begin x_end],[y_begin y_end],min(x_vlegend));
            x_begin = min(x_vlegend);
        elseif x_end>max(x_vlegend)
            y_begin = interp1([x_begin x_end],[y_begin y_end],max(x_vlegend));
            x_begin = max(x_vlegend);
        elseif y_end<min(y_vlegend)
            x_begin = interp1([y_begin y_end],[x_begin x_end],min(y_vlegend));
            y_begin = min(y_vlegend);
        elseif y_end>max(y_vlegend)
            x_begin = interp1([y_begin y_end],[x_begin x_end],max(y_vlegend));
            y_begin = max(y_vlegend);
        end
    end
    ar=annotation(ltype,[x_begin x_end],[y_begin y_end]);
    ar.Color = [1 1 1];
    ar.LineWidth = 0.5;
    ar.LineStyle = '-';
    if strcmp(ltype,'arrow')
        ar.HeadLength = 8;
        ar.HeadWidth = 8;
        ar.HeadStyle = 'vback3';
    end
end
% add vector legend
fill(lon_vlegend([1 2 2 1]),lat_vlegend([1 1 2 2]),'w','EdgeColor','k');
x_begin = interp1(lonbox,ax_position(1)+[0 ax_position(3)],lon_vlegend(1)+0.2,'linear');
x_end = interp1(lonbox,ax_position(1)+[0 ax_position(3)],lon_vlegend(1)+0.2+vlegend*vdegscale/vscale,'linear','extrap');
y_begin = interp1(latbox,ax_position(2)+[0 ax_position(4)],mean(lat_vlegend),'linear');
y_end = interp1(latbox,ax_position(2)+[0 ax_position(4)],mean(lat_vlegend),'linear','extrap');
ar=annotation('arrow',[x_begin x_end],[y_begin y_end]);
ar.Color = 'k';
ar.LineWidth = 2;
ar.LineStyle = '-';
ar.HeadLength = 10;
ar.HeadWidth = 12;
ar.HeadStyle = 'vback2';
text(lon_vlegend(1)+0.2+1.2,mean(lat_vlegend),sprintf('%d m s^{-1}',vlegend),'Color','k','FontSize',12,'HorizontalAlignment','left','VerticalAlignment','middle');
title([num2str(min(p_shear)),'-',num2str(max(p_shear)),'-hPa Shear: ',num2str(yymin),'-',num2str(yymax),' ',mon_txt]);
figname = [figdir,'/map_MERRA2_',num2str(min(p_shear)),'-',num2str(max(p_shear)),'shear_',mon_txt,'climo',num2str(yymin),'-',num2str(yymax)];
fig=gcf; fig.InvertHardcopy = 'off'; fig.Color='w';
print('-painters','-dpng',[figname,'.png']);

% map of qv
ticks_list = 6:0.2:18;
cticks_list = 6:2:18;
units_fill = 'g kg^{-1}';
vscale_list = 5;
vlegend_list = 5;
% wind scale
vint = 2; % plotting interval of merra2 wind vectors
vdegscale = 1; % specific which speed in m/s should correspond to lat/lon deg
lon_vlegend = lonbox(1)+[0 3]+0.1;
lat_vlegend = latbox(1)+[0 1.5]+0.1;
% plotting variables
plot_u_climo = squeeze(umean_climo(ilev_qv,:,:));
plot_v_climo = squeeze(vmean_climo(ilev_qv,:,:));
plot_fill_climo = 10^3*squeeze(qmean_climo(ilev_qv,:,:));
for n=1:length(grididx_pts)
    fprintf('%d hPa qv at lat %0.1f, lon %0.1f: %0.2f g/kg\n',p_qv(k),lat_pts_merra2(n),lon_pts_merra2(n),plot_fill_climo(grididx_pts(n)));
end
% colorbar
ticks_fill = ticks_list;
grad_fill = wetcolor(length(ticks_fill)-1);
cticks_fill = cticks_list{k};
plot_fill_climo(plot_fill_climo<ticks_fill(1))=ticks_fill(1);
plot_fill_climo(plot_fill_climo>ticks_fill(end))=ticks_fill(end);
% wind scale
vscale = vscale_list;
vlegend = vlegend_list;
figure('Position',[1 1 600 425]); hold on
contourf(lon_merra2,lat_merra2,squeeze(plot_fill_climo),ticks_fill,'LineStyle','none');
colormap(gca,grad_fill); clim(ticks_fill([1 end]));
contour(lonTOPO,latTOPO,topo,[0.5 0.5],'-','LineWidth',1,'Color',[0.3 0.3 0.3]);
set(gca,'FontSize',16,'XLim',lonbox,'YLim',latbox,'XTick',[],'YTick',[],'Color',[0.7 0.7 0.7]);
text(lontick,min(latbox)*ones(size(lontick)),lonticklbl,'HorizontalAlignment','center','VerticalAlignment','top','FontSize',14);
text(min(lonbox)*ones(size(lattick)),lattick,latticklbl,'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',14);
cbar=colorbar('YTick',cticks_fill,'YLim',ticks_fill([1 end]));
cbar.Label.String = units_fill;
ax_position = get(gca,'Position');
% add wind vectors
[X,Y] = meshgrid(lon_merra2(1:vint:end),lat_merra2(1:vint:end));
plot_u=plot_u_climo(1:vint:end,1:vint:end); plot_v=plot_v_climo(1:vint:end,1:vint:end);
X=X(:); Y=Y(:); plot_u=plot_u(:); plot_v=plot_v(:);
x_vlegend = interp1(lonbox,ax_position(1)+[0 ax_position(3)],lon_vlegend,'linear');
y_vlegend = interp1(latbox,ax_position(2)+[0 ax_position(4)],lat_vlegend,'linear');
for n=1:length(X)
    ltype = 'arrow';
    if X(n)<min(lonbox) || X(n)>max(lonbox) || Y(n)<min(latbox) || Y(n)>max(latbox)
        continue;
    end
    x_begin = interp1(lonbox,ax_position(1)+[0 ax_position(3)],X(n),'linear');
    x_end = interp1(lonbox,ax_position(1)+[0 ax_position(3)],X(n)+plot_u(n)*vdegscale/vscale,'linear','extrap');
    y_begin = interp1(latbox,ax_position(2)+[0 ax_position(4)],Y(n),'linear');
    y_end = interp1(latbox,ax_position(2)+[0 ax_position(4)],Y(n)+plot_v(n)*vdegscale/vscale,'linear','extrap');
    if y_end>(ax_position(2)+ax_position(4)), y_end = ax_position(2)+ax_position(4); ltype = 'line'; end
    if y_end<ax_position(2), y_end = ax_position(2); ltype = 'line'; end
    if x_end>(ax_position(1)+ax_position(3)), x_end = ax_position(1)+ax_position(3); ltype = 'line'; end
    if x_end<ax_position(1), x_end = ax_position(1); ltype = 'line'; end
    if ((x_begin>=min(x_vlegend)&&x_begin<=max(x_vlegend)) && (x_end>=min(x_vlegend)&&x_end<=max(x_vlegend))) && ...
            ((y_begin>=min(y_vlegend)&&y_begin<=max(y_vlegend)) && (y_end>=min(y_vlegend)&&y_end<=max(y_vlegend)))
        continue;
    end
    if (x_end>=min(x_vlegend)&&x_end<=max(x_vlegend)) && (y_end>=min(y_vlegend)&&y_end<=max(y_vlegend))
        ltype = 'line';
        if x_begin<min(x_vlegend)
            y_end = interp1([x_begin x_end],[y_begin y_end],min(x_vlegend));
            x_end = min(x_vlegend);
        elseif x_begin>max(x_vlegend)
            y_end = interp1([x_begin x_end],[y_begin y_end],max(x_vlegend));
            x_end = max(x_vlegend);
        elseif y_begin<min(y_vlegend)
            x_end = interp1([y_begin y_end],[x_begin x_end],min(y_vlegend));
            y_end = min(y_vlegend);
        elseif y_begin>max(y_vlegend)
            x_end = interp1([y_begin y_end],[x_begin x_end],max(y_vlegend));
            y_end = max(y_vlegend);
        end
    end
    if (x_begin>=min(x_vlegend)&&x_begin<=max(x_vlegend)) && (y_begin>=min(y_vlegend)&&y_begin<=max(y_vlegend))
        if x_end<min(x_vlegend)
            y_begin = interp1([x_begin x_end],[y_begin y_end],min(x_vlegend));
            x_begin = min(x_vlegend);
        elseif x_end>max(x_vlegend)
            y_begin = interp1([x_begin x_end],[y_begin y_end],max(x_vlegend));
            x_begin = max(x_vlegend);
        elseif y_end<min(y_vlegend)
            x_begin = interp1([y_begin y_end],[x_begin x_end],min(y_vlegend));
            y_begin = min(y_vlegend);
        elseif y_end>max(y_vlegend)
            x_begin = interp1([y_begin y_end],[x_begin x_end],max(y_vlegend));
            y_begin = max(y_vlegend);
        end
    end
    ar=annotation(ltype,[x_begin x_end],[y_begin y_end]);
    ar.Color = [0 0 0];
    ar.LineWidth = 1;
    ar.LineStyle = '-';
    if strcmp(ltype,'arrow')
        ar.HeadLength = 8;
        ar.HeadWidth = 8;
        ar.HeadStyle = 'vback3';
    end
end
% add vector legend
fill(lon_vlegend([1 2 2 1]),lat_vlegend([1 1 2 2]),'w','EdgeColor','k');
x_begin = interp1(lonbox,ax_position(1)+[0 ax_position(3)],lon_vlegend(1)+0.2,'linear');
x_end = interp1(lonbox,ax_position(1)+[0 ax_position(3)],lon_vlegend(1)+0.2+vlegend*vdegscale/vscale,'linear','extrap');
y_begin = interp1(latbox,ax_position(2)+[0 ax_position(4)],mean(lat_vlegend),'linear');
y_end = interp1(latbox,ax_position(2)+[0 ax_position(4)],mean(lat_vlegend),'linear','extrap');
ar=annotation('arrow',[x_begin x_end],[y_begin y_end]);
ar.Color = 'k';
ar.LineWidth = 2;
ar.LineStyle = '-';
ar.HeadLength = 10;
ar.HeadWidth = 12;
ar.HeadStyle = 'vback2';
text(lon_vlegend(1)+0.2+1.2,mean(lat_vlegend),sprintf('%d m s^{-1}',vlegend),'Color','k','FontSize',10,'HorizontalAlignment','left','VerticalAlignment','middle');
title([num2str(p_qv(k)),'-hPa Winds and Specific Humidity: ',num2str(yymin),'-',num2str(yymax),' ',mon_txt]);
figname = [figdir,'/map_MERRA2_',num2str(p_qv),'wind_q_',mon_txt,'climo',num2str(yymin),'-',num2str(yymax)];
fig=gcf; fig.InvertHardcopy = 'off'; fig.Color='w';
print('-painters','-dpng',[figname,'.png']);


%% plot cumulative probability distribution

% list of values at the specified points
ushear_pts_list = squeeze(u_list(:,ilev_shear(1),:)-u_list(:,ilev_shear(2),:));
vshear_pts_list = squeeze(v_list(:,ilev_shear(1),:)-v_list(:,ilev_shear(2),:));
% shear magnitude and direction (0:north, 90:east, 180:south, 270:west)
shearmag_pts_list = (ushear_pts_list.^2 + vshear_pts_list.^2).^0.5;
sheardir_pts_list = atan2(-1*ushear_pts_list,-1*vshear_pts_list)*(180/pi);
sheardir_pts_list(sheardir_pts_list<0) = sheardir_pts_list(sheardir_pts_list<0)+360;

figure('Position',[1 1 780 600]);
ytick0=0:0.2:1;
alphabet={'a','b','c','d'};
pwdt=[0.07 0.54]; phgt=[0.56 0.10]; wdt0=0.38; hgt0=0.32;
% shear magnitude
var_list = shearmag_pts_list;
var_edges = 0:2:30; xtick0=0:6:30; 
var_title = [num2str(min(p_shear)),'-',num2str(max(p_shear)),'-hPa Shear Magnitude'];
var_units = 'm s^{-1}';
var_obs = shear_cpexcvobs;
for n=1:2
    subplot('Position',[pwdt(1) phgt(n) wdt0 hgt0]); hold on; clear h
    plot(sort(var_list(:,n)),(1:size(var_list,1))./size(var_list,1),'-b','LineWidth',2);
    set(gca,'XLim',var_edges([1 end]),'XTick',xtick0,'YLim',[0 1],'YTick',0:0.2:1,'FontSize',14,'YTick',ytick0);
    xlabel(var_units);
    if n==1
        text(mean(var_edges([1 end])),1.18,var_title,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16,'FontWeight','bold');
    end
    title({sprintf('%s) Lat %0.1f%sN Lon %0.1f%sW',alphabet{2*(n-1)+1},lat_pts_merra2(n),char(176),abs(lon_pts_merra2(n)),char(176))});
    grid on
    % add obs value
    var_perc = sum(var_obs(n)>var_list(:,n))./size(var_list,1);
    var_fill = sort(var_list(:,n));
    var_fill = var_fill(var_obs(n)>var_fill);
    h=fill([var_fill;var_fill([end 1])],[(1:length(var_fill)) 0 0]./size(var_list,1),'b','FaceAlpha',0.2,'EdgeColor','none');
    plot(var_obs(n)*[1 1],[0 var_perc],'--b','LineWidth',1);
    text(var_obs(n)+0.01*(var_edges(end)-var_edges(1)),0,{sprintf('%0.2f %s',var_obs(n),var_units),sprintf('(%dth percentile)',round(var_perc*100))},'HorizontalAlignment','left','VerticalAlignment','bottom','Color','b');
end
% specific humidity
p0 = 925;
var_list = 10^3*squeeze(q_list(:,lev_merra2==p0,:));
var_edges = 4:1:20; xtick0=4:4:20;
var_title = [num2str(p0),'-hPa Specific Humidity'];
var_units = 'g kg^{-1}';
var_obs = q_cpexcvobs;
for n=1:2
    subplot('Position',[pwdt(2) phgt(n) wdt0 hgt0]); hold on
    plot(sort(var_list(:,n)),(1:size(var_list,1))./size(var_list,1),'-b','LineWidth',2);
    set(gca,'XLim',var_edges([1 end]),'XTick',xtick0,'YLim',[0 1],'YTick',0:0.2:1,'FontSize',14,'YTick',ytick0);
    xlabel(var_units);
    if n==1
        text(mean(var_edges([1 end])),1.18,var_title,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',16,'FontWeight','bold');
    end
    title({sprintf('%s) Lat %0.1f%sN Lon %0.1f%sW',alphabet{2*n},lat_pts_merra2(n),char(176),abs(lon_pts_merra2(n)),char(176))});
    grid on
    % add obs value
    var_perc = sum(var_obs(n)>var_list(:,n))./size(var_list,1);
    var_fill = sort(var_list(:,n));
    var_fill = var_fill(var_obs(n)>var_fill);
    h=fill([var_fill;var_fill([end 1])],[(1:length(var_fill)) 0 0]./size(var_list,1),'b','FaceAlpha',0.2,'EdgeColor','none');
    plot(var_obs(n)*[1 1],[0 var_perc],'--b','LineWidth',1);
    text(var_obs(n)+0.01*(var_edges(end)-var_edges(1)),0,{sprintf('%0.2f %s',var_obs(n),var_units),sprintf('(%dth percentile)',round(var_perc*100))},'HorizontalAlignment','left','VerticalAlignment','bottom','Color','b');
end
gridpts_txt=sprintf('%0.1fN%0.1fW_%0.1fN%0.1fW',lat_pts_merra2(1),abs(lon_pts_merra2(1)),lat_pts_merra2(2),abs(lon_pts_merra2(2)));
figname = [figdir,'/cumdistributions_MERRA2_',gridpts_txt,'_',mon_txt,'climo',num2str(yymin),'-',num2str(yymax)];
fig=gcf; fig.InvertHardcopy = 'off'; fig.Color='w';
print('-painters','-dpng',[figname,'.png']);
